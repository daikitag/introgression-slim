import numpy as np
import pandas as pd
import tskit
import tstrait

snk = snakemake  # noqa: F821

# Load saved tree sequence
ts = tskit.load(snk.input[0])

num_ind = int(snk.params.N)

assert ts.num_mutations > 0, "Mutations are lost in the SLiM simulation output"

# All mutations are happening at the same time, so we can obtain the node for the 1st
# mutation
# Time is in forward simulation
mutation_list = []
for i in range(ts.num_sites):
    mutation = ts.site(i).mutations[0]
    mutation_list.append(
        {
            "site_id": i,
            "effect_size": mutation.metadata["mutation_list"][0]["selection_coeff"],
            "causal_allele": mutation.derived_state,
            "trait_id": 0,
        }
    )

mutation_df = pd.DataFrame(mutation_list)
genetic_df = tstrait.genetic_value(ts, mutation_df)

trait_list = []
for time in range(int(ts.max_root_time)):
    individual_list = []
    node_list = np.where((ts.nodes_time == time) & (ts.nodes_individual != -1))[0]
    for j in node_list:
        individual_id = ts.node(j).individual
        if individual_id not in individual_list:
            individual_list.append(individual_id)
    genetic_select = genetic_df.iloc[individual_list]
    trait_list.append(
        {
            "time": time,
            "trait_average": np.sum(genetic_select["genetic_value"]) / num_ind,
        }
    )

trait_df = pd.DataFrame(trait_list)

# Convert time into forwards in time, where (time=0) is when the simulation starts
trait_df["time"] = np.max(trait_df["time"]) - trait_df["time"] + 1

trait_df.to_csv(snk.output[0], index=False)
