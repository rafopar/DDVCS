#!/usr/bin/env python3
import sys
import numpy as np
import ROOT
import torch
import json
from torch import nn
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, accuracy_score

# ============================================================
#  LOAD ROOT → NUMPY
# ============================================================
def load_root(filename, treename="tr1"):
    f = ROOT.TFile(filename)
    t = f.Get(treename)

    X_list = []
    y_list = []

    for i in range(t.GetEntries()):
        t.GetEntry(i)

        try:
            P_em1 = t.P_em1
            th_em1 = t.th_em1
            phi_em1 = t.phi_em1

            P_em2 = t.P_em2
            th_em2 = t.th_em2
            phi_em2 = t.phi_em2

            P_ep = t.P_ep
            th_ep = t.th_ep
            phi_ep = t.phi_ep

            Q2_em1 = t.Q2_em1
            Q2_em2 = t.Q2_em2
        except Exception as e:
            print(f"Skipping event {i}: {e}")

        em1 = (P_em1, th_em1, phi_em1)
        em2 = (P_em2, th_em2, phi_em2)
        ep = (P_ep, th_ep, phi_ep)

        # Random swap for data augmentation
        if np.random.rand() < 0.5:
            em1, em2 = em1, em2
            Q2_em1, Q2_em2 = Q2_em1, Q2_em2
            y = 0
        else:
            em1, em2 = em2, em1
            Q2_em1, Q2_em2 = Q2_em2, Q2_em1
            y = 1

        X_list.append([*em1, *em2, *ep, Q2_em1, Q2_em2])
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
            nn.Linear(n_inputs, 32), nn.GELU(), nn.Dropout(0.1),
            nn.Linear(32, 32), nn.GELU(),
            nn.Linear(32, 16), nn.GELU(),
            nn.Linear(16, 1), nn.Sigmoid()
        )

    def forward(self, x):
        return self.net(x)


# ============================================================
#  TRAINING LOOP
# ============================================================
def train_model(X_train, y_train, X_val, y_val, epochs=30, batch=1024):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    X_train = torch.tensor(X_train, dtype=torch.float32, device=device)
    y_train = torch.tensor(y_train, dtype=torch.float32, device=device).unsqueeze(1)

    X_val = torch.tensor(X_val, dtype=torch.float32, device=device)
    y_val = torch.tensor(y_val, dtype=torch.float32, device=device).unsqueeze(1)

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
    y_test = torch.tensor(y_test, dtype=torch.float32, device=device).unsqueeze(1)

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

    print("Running the %s" %(__name__));

    if len(sys.argv) < 2:
        print("Usage: python3 train_GEMC_elDDCVC.py $RUN")
        sys.exit(1)

    run = int(sys.argv[1])
    fname = "TrainingSample/TrainingSample_Run_%d.root"%(run)

    nEpochs = 1200
    BatchSize = 5000

    print("Loading %s..."%(fname))

    X, y = load_root(fname, treename="tr1")

    print("Size of X is %d" %(len(X)))
    print("Size of y is %d" %(len(y)))

    print("Example input ...")
    print("X[0], y[0]: " + str(X[0]) + " " + str(y[0]))
    print("X[1], y[1]: " + str(X[1]) + " " + str(y[1]))
    print("X[2], y[2]: " + str(X[2]) + " " + str(y[2]))
    print("X[3], y[3]: " + str(X[3]) + " " + str(y[3]))

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

    print("Size of X_train is %d" %(len(X_train)))
    print("Size of X_val is %d" %(len(X_val)))
    print("Size of Y_train is %d" %(len(y_train)))
    print("Size of Y_val is %d" %(len(y_val)))

    # Normalize input variables
    X_train, mean, std = normalize(X_train)
    X_test = (X_test - mean) / std
    X_val = (X_val - mean) / std

    print("mean: " + str(mean))
    print("std: " + str(std))

    model, train_losses, val_losses = train_model(X_train, y_train, X_val, y_val, nEpochs, BatchSize)

    test_pred = evaluate(model, X_test, y_test)

    example_input = torch.randn(1, X_train.shape[1])
    traced = torch.jit.trace(model, example_input)

    # Save TorchScript model
    traced.save("pars/assign_model_ts_Run_%d.pt"%(run))
    print("Saved TorchScript model → assign_model_ts.pt")

    with open("pars/norm_Run_%d.json"%(run), "w") as f:
        json.dump({"mean": mean.tolist(), "std": std.tolist()}, f)

    print("Saved normalization constants → norm.json")

    # save model
    torch.save({
        "model_state": model.state_dict(),
        "mean": mean,
        "std": std
    }, "pars/assign_model_Run_%d.pt"%(run))

    print("\nSaved model → assign_model.pt")


    # save training_history
    torch.save({
        "loss_train": train_losses,
        "loss_val": val_losses,
        "y_true": y_test,
        "y_pred": test_pred
    }, "pars/training_history_Run_%d.pt"%(run))

