import numpy as np
import pandas as pd
import tskit
import tstrait

snk = snakemake  # noqa: F821

# Load saved tree sequence
ts = tskit.load(snk.input[0])

num_ind = int(snk.params.N)

# All mutations are happening at the same time, so we can obtain the node for the 1st
# mutation
# Time is in forward simulation
mutation_node = ts.mutation(0).node
max_time = int(ts.node(mutation_node).time)

trait_df = pd.DataFrame()
for i in range(ts.num_sites):
    mutation = ts.site(i).mutations[0]
    individual_df = pd.DataFrame(
        {
            "site_id": [i],
            "effect_size": [mutation.metadata["mutation_list"][0]["selection_coeff"]],
            "causal_allele": [mutation.derived_state],
            "trait_id": [0],
        }
    )
    trait_df = pd.concat([trait_df, individual_df], ignore_index=True)

genetic_df = tstrait.genetic_value(ts, trait_df)

trait_average = []
for time in range(max_time + 1):
    individual_list = []
    node_list = np.where(ts.nodes_time == time)[0]
    for j in node_list:
        individual_id = ts.node(j).individual
        if individual_id not in individual_list:
            individual_list.append(individual_id)
    genetic_select = genetic_df.iloc[individual_list]
    trait_average.append(np.sum(genetic_select["genetic_value"]) / num_ind)

trait_df = pd.DataFrame(
    {"time": np.arange(max_time + 1), "trait_average": trait_average}
)

trait_df.to_csv(snk.output[0])
