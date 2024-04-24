import matplotlib.pyplot as plt
import pandas as pd

snk = snakemake  # noqa: F821

df = pd.read_csv(snk.input[0])

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
plt.ylabel("Average total length of surviving blocks per individual", fontsize=13)
plt.savefig(snk.output.total_length)

plt.figure(figsize=(8, 6))
plt.plot(df["time"], df["num_node"], marker=".")
plt.xlabel("Time t", fontsize=13)
plt.ylabel("Number of chromosomes inheriting the block", fontsize=13)
plt.savefig(snk.output.num_node)
