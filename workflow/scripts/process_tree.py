import numpy as np
import pandas as pd
import tskit
from tqdm import tqdm

snk = snakemake  # noqa: F821

ts = tskit.load(snk.input[0])
table = ts.tables

num_ind = int(snk.params.N)

# All mutations are happening at the same time, so we can obtain the node for the 1st
# mutation.
# Time is in forward simulation
mutation_node = ts.mutation(0).node
max_time = int(ts.node(mutation_node).time)

# It would be better to put `plot_block.py` plot code here instead of storing IBD
# segments.
# Store IBD segment length information in csv and write another Python code to plot


total_length_average_array = []
block_length_average_array = []
num_block_array = []

for i in tqdm(range(int(max_time + 1))):
    # (ts.nodes_time==i) & (ts.nodes_individual != -1)
    # It is individual -> Guaranteed that it is the node that carries the mutation
    ibd_table = table.link_ancestors(
        samples=np.where(ts.nodes_time == i)[0], ancestors=[mutation_node]
    )
    total_length = np.sum(ibd_table.right - ibd_table.left)
    total_length_average = total_length / num_ind
    block_length_average = np.mean(ibd_table.right - ibd_table.left)
    num_block = len(ibd_table)

    total_length_average_array.append(total_length_average)
    block_length_average_array.append(block_length_average)
    num_block_array.append(num_block)

df = pd.DataFrame(
    {
        "time": np.arange(max_time + 1),
        "total_length": total_length_average_array,
        "block_length": block_length_average_array,
        "num_block": num_block_array,
    }
)
df.to_csv(snk.output[0])
