import matplotlib.pyplot as plt
import pandas as pd

snk = snakemake  # noqa: F821

trait_df = pd.read_csv(snk.input[0])

plt.figure(figsize=(8, 6))
plt.plot(trait_df["time"], trait_df["trait_average"], marker=".")
plt.xlabel("Time t", fontsize=13)
plt.ylabel("Average trait value per individual", fontsize=13)
plt.savefig(snk.output[0])
