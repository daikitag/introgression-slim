import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

snk = snakemake  # noqa: F821

trait_df = pd.read_csv(snk.input[0])

# Convert time into forwards in time
trait_df["time"] = np.max(trait_df["time"]) - trait_df["time"]

plt.figure(figsize=(8, 6))
plt.plot(trait_df["time"], trait_df["trait_average"], marker=".")
plt.xlabel("Time t", fontsize=13)
plt.ylabel("Average trait value in population", fontsize=13)
plt.savefig(snk.output[0])
