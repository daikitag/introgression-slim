import numpy as np
import pandas as pd
import tskit
from tqdm import tqdm

snk = snakemake  # noqa: F821

ts = tskit.load(snk.input[0])
table = ts.tables

num_ind = int(snk.params.N)

assert ts.num_mutations > 0, "Mutations are lost in the SLiM simulation output"

# All mutations are happening at the same time, so we can obtain the node for the 1st
# mutation.
# Time is in forward simulation
mutation_node = ts.mutation(0).node

# It would be better to put `plot_block.py` plot code here instead of storing IBD
# segments.
# Store IBD segment length information in csv and write another Python code to plot

ibd_df = []

for i in tqdm(range(int(ts.node(mutation_node).time))):
    # (ts.nodes_time==i) & (ts.nodes_individual != -1)
    # It is individual -> Guaranteed that it is the node that carries the mutation
    node_array = np.where((ts.nodes_time == i) & (ts.nodes_individual != -1))[0]

    if node_array.size > 0:
        ibd_table = table.link_ancestors(samples=node_array, ancestors=[mutation_node])
        total_length = np.sum(ibd_table.right - ibd_table.left)
        total_length_average = total_length / num_ind
        block_length_average = np.mean(ibd_table.right - ibd_table.left)
        num_block = len(ibd_table)
    else:
        total_length_average = 0
        block_length_average = 0
        num_block = 0

    ibd_df.append(
        {
            "time": i,
            "total_length": total_length_average,
            "block_length": block_length_average,
            "num_block": num_block,
            "num_node": len(np.unique(ibd_table.child)),
        }
    )

df = pd.DataFrame(ibd_df)

# Convert time into forwards in time, where (time=0) is when the simulation starts
df["time"] = np.max(df["time"]) - df["time"] + 1

df.to_csv(snk.output[0], index=False)
