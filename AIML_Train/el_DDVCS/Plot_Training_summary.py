import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-GUI backend
import matplotlib.pyplot as plt
matplotlib.use("Agg")   # non-GUI backend
from sklearn.metrics import roc_curve, auc, confusion_matrix
import torch

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: python3 ./bin/Plot_Training_summary.py $RUN")
        sys.exit(1)

    run = int(sys.argv[1])

    history = torch.load("pars/training_history_Run_%d.pt"%(run), weights_only=False)

    y_true = np.array(history["y_true"])
    y_pred = np.array(history["y_pred"])
    loss_train = history["loss_train"]
    loss_val   = history["loss_val"]

    # --------------------
    # ROC curve
    # --------------------
    fpr, tpr, thr = roc_curve(y_true, y_pred)
    roc_auc = auc(fpr, tpr)

    plt.figure()
    plt.plot(fpr, tpr, label="AUC = %.4f" % roc_auc)
    plt.plot([0,1],[0,1],'--')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve")
    plt.legend()
    plt.grid()
    plt.savefig("Figs/roc_curve_Run_%d.png"%(run), dpi=200)

    # --------------------
    # Confusion Matrix
    # --------------------
    cm = confusion_matrix(y_true, y_pred > 0.5)

    plt.figure()
    plt.imshow(cm, cmap="Blues")
    plt.colorbar()
    plt.title("Confusion Matrix (threshold=0.5)")
    plt.xlabel("Predicted")
    plt.ylabel("True")

    for (i,j), val in np.ndenumerate(cm):
        plt.text(j, i, f"{val}", ha="center", va="center")

    plt.savefig("Figs/confusion_matrix_Run_%d.png"%(run), dpi=200)


    # --------------------
    # Loss Curves
    # --------------------
    plt.figure()
    plt.plot(loss_train, label="Train")
    plt.plot(loss_val, label="Validation")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.title("Training History")
    plt.ylim([0.5,0.55])
    plt.legend()
    plt.grid()
    plt.savefig("Figs/loss_curve_Run_%d.png"%(run), dpi=200)

    print("Saved plots: roc_curve.png, confusion_matrix.png, loss_curve.png")
