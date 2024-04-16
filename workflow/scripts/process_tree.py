import numpy as np
import pandas as pd
from tqdm import tqdm
import tskit

ts = tskit.load(snakemake.input[0])
table = ts.tables

# All mutations are happening at the same time, so we can obtain the node for the 1st mutation
# Time is in forward simulation
mutation_node = ts.mutation(0).node
max_time = ts.node(mutation_node).time

df = pd.DataFrame()
for i in tqdm(range(int(max_time))):
    ibd_table = table.link_ancestors(samples=np.where(ts.nodes_time==i)[0], ancestors=[mutation_node])
    ibd_df = pd.DataFrame({
        "left": ibd_table.left,
        "right": ibd_table.right,
        "parent": ibd_table.parent,
        "child": ibd_table.child,
        "time": np.ones(len(ibd_table))*i
    })
    df = pd.concat([df, ibd_df])

df.to_csv(snakemake.output[0])