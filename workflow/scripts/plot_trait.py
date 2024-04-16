import matplotlib.pyplot as plt
import msprime
import numpy as np
import pandas as pd
import tskit
import tstrait

# Compile tstrait
sample_ts = msprime.sim_ancestry(
    samples=10,
    recombination_rate=1e-8,
    sequence_length=1000,
    population_size=10_000,
    random_seed=100,
)
sample_ts = msprime.sim_mutations(sample_ts, rate=1e-6, random_seed=101)
model = tstrait.trait_model(distribution="normal", mean=0, var=1)
sim_result = tstrait.sim_phenotype(
    ts=sample_ts, num_causal=2, model=model, h2=0.3, random_seed=1
)

# Load saved tree sequence
ts = tskit.load(snakemake.input[0])

num_ind = int(snakemake.params.N)

# All mutations are happening at the same time, so we can obtain the node for the 1st mutation
# Time is in forward simulation
mutation_node = ts.mutation(0).node
max_time = int(ts.node(mutation_node).time)

trait_df = pd.DataFrame()
for i in range(ts.num_sites):
    mutation = ts.site(i).mutations[0]
    individual_df = pd.DataFrame({
        "site_id": [i],
        "effect_size": [mutation.metadata["mutation_list"][0]["selection_coeff"]],
        "causal_allele": [mutation.derived_state],
        "trait_id": [0]
    })
    trait_df = pd.concat([trait_df, individual_df], ignore_index=True)
    
genetic_df = tstrait.genetic_value(ts, trait_df)

trait_average = []
for i in range(1, max_time+1):
    time = max_time - i
    individual_list = []
    node_list = np.where(ts.nodes_time==time)[0]
    for j in node_list:
        individual_id = ts.node(j).individual
        if individual_id not in individual_list:
            individual_list.append(individual_id)
    genetic_select = genetic_df.iloc[individual_list]
    trait_average.append(np.sum(genetic_select["genetic_value"])/num_ind)

plt.figure(figsize=(8,6))
plt.plot(np.arange(len(trait_average)), trait_average, marker=".")
#axs.set_xscale("log", base=10)
plt.xlabel('Time t', fontsize=13)
plt.ylabel('Average trait value in population', fontsize=13)
plt.savefig(snakemake.output[0])