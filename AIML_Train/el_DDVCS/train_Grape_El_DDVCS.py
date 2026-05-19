#!/usr/bin/env python3
import sys
import numpy as np
import ROOT
import torch
import json
from torch import nn
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, accuracy_score
from zmq import IN_BATCH_SIZE

# ============================================================
#  PHYSICS CONSTANTS
# ============================================================
m_e = 0.000511  # GeV


# ============================================================
#  PHYSICS UTILITIES
# ============================================================
def invariant_mass(px1, py1, pz1, px2, py2, pz2):
    """Compute invariant mass of two particles with electron mass."""
    E1 = np.sqrt(px1**2 + py1**2 + pz1**2 + m_e**2)
    E2 = np.sqrt(px2**2 + py2**2 + pz2**2 + m_e**2)
    E = E1 + E2
    px = px1 + px2
    py = py1 + py2
    pz = pz1 + pz2
    M2 = E**2 - px**2 - py**2 - pz**2
    if M2 < 0:
        M2 = 0

    return np.sqrt(M2)


# ============================================================
#  LOAD ROOT → NUMPY
# ============================================================
def load_root(filename, treename="tree"):
    f = ROOT.TFile(filename)
    t = f.Get(treename)

    X_list = []
    y_list = []

    for i in range(t.GetEntries()):
        t.GetEntry(i)

        # Extract particles:
        # em1 = index 1
        # ep  = index 2
        # em2 = index 3


        # Convert directly using Python lists
        try:
            px_list = list(t.px)
            py_list = list(t.py)
            pz_list = list(t.pz)
        except Exception as e:
            print(f"Skipping event {i}: {e}")
            continue

        # Ensure the event has enough particles
        if len(px_list) < 4:
            continue  # skip incomplete events


        # Extract particle components
        em1 = (px_list[1], py_list[1], pz_list[1])
        ep  = (px_list[2], py_list[2], pz_list[2])
        em2 = (px_list[3], py_list[3], pz_list[3])

        #print("Index = %d: Electron 1: %1.2f   %1.2f   %1.2f"%(i, px_list[1], py_list[1], pz_list[1]))

        # --- Build the two possible orderings for training ---
        # probability 0.5 to swap (data augmentation)


        # Random swap for data augmentation
        if np.random.rand() < 0.5:
            p1, p2 = em1, em2

            y = 0
        else:
            p1, p2 = em2, em1
            y = 1

        # Compute invariant masses
        m1p = invariant_mass(p1[0], p1[1], p1[2], ep[0], ep[1], ep[2])
        m2p = invariant_mass(p2[0], p2[1], p2[2], ep[0], ep[1], ep[2])

        #m1p = 1
        #m2p = 2


        # Build feature vector
        X_list.append([
            *p1, *p2, *ep,   # px, py, pz for p1, p2, ep
            m1p, m2p
        ])
        y_list.append(y)

    return np.array(X_list, dtype=float), np.array(y_list, dtype=int)


# ============================================================
#  NORMALIZATION
# ============================================================
def normalize(X):
    mean = X.mean(axis=0)
    std = X.std(axis=0)
    Xn = (X - mean) / std
    return Xn, mean, std


# ============================================================
#  MODEL DEFINITION
# ============================================================
class AssignNet(nn.Module):
    def __init__(self, n_inputs=11):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_inputs, 64), nn.GELU(), nn.Dropout(0.1),
            nn.Linear(64, 64), nn.GELU(),
            nn.Linear(64, 32), nn.GELU(),
            nn.Linear(32, 1), nn.Sigmoid()
        )

    def forward(self, x):
        return self.net(x)

# class AssignNet(nn.Module):
#     def __init__(self, n_inputs=11):
#         super().__init__()
#         self.net = nn.Sequential(
#             nn.Linear(n_inputs, 128), nn.ReLU(),
#             nn.Linear(128, 128), nn.ReLU(),
#             nn.Linear(128, 64), nn.ReLU(),
#             nn.Linear(64, 1), nn.Sigmoid()
#         )
#
#     def forward(self, x):
#         return self.net(x)


# ============================================================
#  TRAINING LOOP
# ============================================================
def train_model(X_train, y_train, X_val, y_val, epochs=30, batch=1024):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    X_train = torch.tensor(X_train, dtype=torch.float32, device=device)
    y_train = torch.tensor(y_train, dtype=torch.float32,
                           device=device).unsqueeze(1)

    X_val = torch.tensor(X_val, dtype=torch.float32, device=device)
    y_val = torch.tensor(y_val, dtype=torch.float32,
                         device=device).unsqueeze(1)

    model = AssignNet().to(device)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    loss_fn = nn.BCELoss()

    nbatches = max(1, len(X_train) // batch)

    train_losses = []
    val_losses = []

    for epoch in range(epochs):

        model.train()
        perm = torch.randperm(len(X_train))
        X_train = X_train[perm]
        y_train = y_train[perm]
        batch_losses = []

        for i in range(nbatches):
            xb = X_train[i*batch:(i+1)*batch]
            yb = y_train[i*batch:(i+1)*batch]

            pred = model(xb)
            loss = loss_fn(pred, yb)

            opt.zero_grad()
            loss.backward()
            opt.step()
            batch_losses.append(loss.item())

        # store training loss
        train_losses.append(np.mean(batch_losses))
        #train_losses.append(loss.item())

        # validation
        model.eval()
        with torch.no_grad():
            val_pred = model(X_val)
            val_loss = loss_fn(val_pred, y_val).item()
            val_losses.append(val_loss)

        print(f"Epoch {epoch+1:02d} / {epochs},  val loss = {val_loss:.5f}")

    return model, train_losses, val_losses


# ============================================================
#  EVALUATION
# ============================================================
def evaluate(model, X_test, y_test):
    device = next(model.parameters()).device

    X_test = torch.tensor(X_test, dtype=torch.float32, device=device)
    y_test = torch.tensor(y_test, dtype=torch.float32,
                          device=device).unsqueeze(1)

    model.eval()
    with torch.no_grad():
        pred = model(X_test).cpu().numpy()

    auc = roc_auc_score(y_test.cpu().numpy(), pred)
    acc = accuracy_score(y_test.cpu().numpy(), pred > 0.5)

    print("\n===== Test set performance =====")
    print("AUC      :", auc)
    print("Accuracy :", acc)

    return pred;


# ============================================================
#  MAIN
# ============================================================

# This program is intended to train on Electron DDVCS events from the GRAPE to be able to distinguish
# which electron is the beam electron and which one is the pait with the positron

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: python3 train_Grape_El_DDVCS.py $RUN")
        sys.exit(1)

    run = int(sys.argv[1])
    fname = "inpData/Train_GRAPE_Run%d.root"%run
    print("Loading:", fname)

    nEpochs = 60
    BatchSize = 2500
    print("Will do %d epochs"%(nEpochs))

    X, y = load_root(fname, "tr1")

    print("Size of X is %d"%(len(X)))
    print("Size of Y is %d"%(len(y)))

    # train/val/test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.15, random_state=42
    )

    print("Size of X_train is %d" %(len(X_train)))
    print("Size of X_test is %d" %(len(X_test)))
    print("Size of Y_train is %d" %(len(y_train)))
    print("Size of Y_test is %d" %(len(y_test)))


    X_train, X_val, y_train, y_val = train_test_split(
        X_train, y_train, test_size=0.15, random_state=42
    )

    example_input = torch.randn(1, X_train.shape[1])

    print("Example input ...")
    print(example_input)

    # normalize
    X_train, mean, std = normalize(X_train)
    X_val = (X_val - mean) / std
    X_test = (X_test - mean) / std

    model, train_losses, val_losses = train_model(X_train, y_train, X_val, y_val, nEpochs, BatchSize)

    test_pred = evaluate(model, X_test, y_test)


    traced = torch.jit.trace(model, example_input)


    # Save TorchScript model
    traced.save("assign_model_ts.pt")
    print("Saved TorchScript model → assign_model_ts.pt")

    with open("norm.json", "w") as f:
        json.dump({"mean": mean.tolist(), "std": std.tolist()}, f)

    print("Saved normalization constants → norm.json")

    # save model
    torch.save({
        "model_state": model.state_dict(),
        "mean": mean,
        "std": std
    }, "assign_model.pt")

    print("\nSaved model → assign_model.pt")


    # save training_history
    torch.save({
        "loss_train": train_losses,
        "loss_val": val_losses,
        "y_true": y_test,
        "y_pred": test_pred
    }, "training_history.pt")
