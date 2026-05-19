import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.utils import plot_model
import os

os.makedirs("Figs", exist_ok=True)

model = keras.Sequential([
    layers.Input(shape=(11,)),
    layers.Dense(32, activation="relu"),
    layers.Dense(32, activation="relu"),
    layers.Dense(16, activation="relu"),
    layers.Dense(1, activation="sigmoid"),
])

model.summary()

plot_model(
    model,
    to_file="Figs/nn_architecture.png",
    show_shapes=True,
    show_layer_names=True,
    expand_nested=False,
)

print("Saved: Figs/nn_architecture.png")