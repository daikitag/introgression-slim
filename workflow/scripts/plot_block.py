import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

snk = snakemake  # noqa: F821

df = pd.read_csv(snk.input[0])

# Convert time into forwards in time
df["time"] = np.max(df["time"]) - df["time"]

plt.figure(figsize=(8, 6))
plt.plot(df["time"], df["num_block"], marker=".")
plt.xlabel("Time t", fontsize=13)
plt.ylabel("Number of surviving blocks", fontsize=13)
plt.savefig(snk.output.block_number)

plt.figure(figsize=(8, 6))
plt.plot(df["time"], df["block_length"], marker=".")
plt.xlabel("Time t", fontsize=13)
plt.ylabel("Average length of surviving blocks", fontsize=13)
plt.savefig(snk.output.block_length)

plt.figure(figsize=(8, 6))
plt.plot(df["time"], df["total_length"], marker=".")
plt.xlabel("Time t", fontsize=13)
plt.ylabel("Average total length of surviving blocks", fontsize=13)
plt.savefig(snk.output.total_length)
