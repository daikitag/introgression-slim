import matplotlib.pyplot as plt
import pandas as pd

snk = snakemake  # noqa: F821

ibd_df = pd.read_csv(snk.input.ibd_df)
trait_df = pd.read_csv(snk.input.trait_df)

figure, axs = plt.subplots(3, 2, figsize=(16, 18))
axs[0, 0].plot(ibd_df["time"], ibd_df["num_block"], marker=".")
axs[0, 0].set_xlabel("Time t", fontsize=13)
axs[0, 0].set_ylabel("Number of surviving blocks", fontsize=13)

axs[0, 1].plot(ibd_df["time"], ibd_df["block_length"], marker=".")
axs[0, 1].set_xlabel("Time t", fontsize=13)
axs[0, 1].set_ylabel("Average length of surviving blocks", fontsize=13)

axs[1, 0].plot(ibd_df["time"], ibd_df["total_length"], marker=".")
axs[1, 0].set_xlabel("Time t", fontsize=13)
axs[1, 0].set_ylabel(
    "Average total length of surviving blocks per individual", fontsize=13
)

axs[1, 1].plot(ibd_df["time"], ibd_df["num_node"], marker=".")
axs[1, 1].set_xlabel("Time t", fontsize=13)
axs[1, 1].set_ylabel("Number of chromosomes inheriting the block", fontsize=13)

axs[2, 0].plot(trait_df["time"], trait_df["trait_average"], marker=".")
axs[2, 0].set_xlabel("Time t", fontsize=13)
axs[2, 0].set_ylabel("Average trait value per individual", fontsize=13)

axs[2, 1].axis("off")

plt.savefig(snk.output[0])
