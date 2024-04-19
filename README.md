# Simulating Introgression in SLiM

This repository conducts the introgression model simulation described in [Sachdeva and Barton (2018)](https://doi.org/10.1534/genetics.118.301018). This repository conducts an individual based forward-time simulation in [SLiM](https://messerlab.org/slim/).

The parameters of the simulation model can be set by manipulating the parameters in `config.yaml` file inside `config` directory, and the details of the parameters are written in `README.md` file in the directory. See the `Mathematical model` section below for details of the simulation.

## Usage:

1. Clone this repo

```
$ git clone https://github.com/daikitag/introgression-slim.git
```

2. Install `snakemake` by following the instructions written here: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html. If you have `coonda` installed on your computer, you can install `snakemake` by running:

```
$ conda install -c bioconda snakemake
```
You do not have to install other packages like `SLiM` and `tskit`, as `snakemake` will automatically install them for you when you run the command.

3. Run the simulation command by

```
$ snakemake --cores 1 --use-conda
```

Here, `--cores` command specifies the number of cores to use and `--use-conda` allows `snakemake` to be run in a `conda` environment. This `snakemake` module uses a `conda` environment to run the simulations, so please do not forget to include this argument.

## Output:

The output plots of the simulation will be generated in `output` folder. The parameters of the simulation are indicated in the file name.

## Note:

Running this simulation code for the first time will require a lot of time, as it creates a new conda environment with the required Python packages for the simulation. It will not take a lot of time for the simulation afterwards, as conda environment is already constructed in the local repository.

# Mathematical model:

The mathematical details of the simulation model is described in [Sachdeva and Barton (2018)](https://doi.org/10.1534/genetics.118.301018) as "individual-based simulation" (See section named "Individual-based simulations of long-term introgression into a finite population" in the paper for the details of the simulation). In particular, we simulate introgression based on the following assumptions:

- Discrete Wright-Fisher model
- Effect sizes and fitness are defined as the "Individual-based simulations of long-term introgression into a finite population" section
- Recombination rate is set as rate/site, instead of rate/individual in [Sachdeva and Barton (2018)](https://doi.org/10.1534/genetics.118.301018). It would be possible for the simulation to have multiple recombination breakpoints.
- Discrete loci (Recombination breakpoints will be an integer)
- Zero mutation rate, and all mutations descend from the single individual at the first generation

The output plots are obtained by calculating the summary statistics of the IBD segments that are descended from the individual carrying the mutation at the beginning of the simulation.

# References:

- Simulation model is taken from [Sachdeva and Barton (2018)](https://doi.org/10.1534/genetics.118.301018)
- The workflow of this Github repo is adapted from the structure given in the [Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows)
- The basic structure of using snakemake in SLiM is adapted from the [snakemake-tutorial](https://github.com/vsbuffalo/snakemake-tutorial/tree/master) Github repository

# Acknowledgements:

I would like to thank Jerome Kelleher for supervising and providing feedback on this Github repository, Peter Ralph for providing the initial codes for simulating introgression in SLiM and suggesting the usage of `link_ancestors` function in `tskit` to analyze the SLiM output, Ben Haller for providing suggestions on modeling recombination and effect size distributions in SLiM, and people in the tskit community for valuable suggestions regarding the usage of snakemake for combining `SLiM` and `tskit` libraries.