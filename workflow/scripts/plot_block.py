import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tskit

ts = tskit.load(snakemake.input.tree)
df = pd.read_csv(snakemake.input.ibd_df)
df["length"] = df["right"] - df["left"]

# All mutations are happening at the same time, so we can obtain the node for the 1st mutation
# Time is in forward simulation
mutation_node = ts.mutation(0).node
max_time = int(ts.node(mutation_node).time)

num_ind = int(snakemake.params.N)

total_length_average = []
length_average = []
num_block = []
for i in range(1, max_time+1):
    sample_df = df[df["time"]==max_time-i]
    average = np.sum(sample_df.groupby("child").sum()["length"])/num_ind
    total_length_average.append(average)
    num_block.append(len(sample_df))
    length_average.append(np.mean(sample_df["length"]))

plt.figure(figsize=(8,6))
plt.plot(np.arange(len(num_block)), num_block, marker=".")
plt.xlabel('Time t', fontsize=13)
plt.ylabel('Number of surviving blocks', fontsize=13)
plt.savefig(snakemake.output.block_number)

plt.figure(figsize=(8,6))
plt.plot(np.arange(len(length_average)), length_average, marker=".")
plt.xlabel('Time t', fontsize=13)
plt.ylabel('Average length of surviving blocks', fontsize=13)
plt.savefig(snakemake.output.block_length)

plt.figure(figsize=(8,6))
plt.plot(np.arange(len(total_length_average)), total_length_average, marker=".")
plt.xlabel('Time t', fontsize=13)
plt.ylabel('Average total length of surviving blocks', fontsize=13)
plt.savefig(snakemake.output.total_length)