# Simulating Introgression in SLiM

This repository conducts the introgression model simulation described in [Sachdeva and Barton (2018)](https://doi.org/10.1534/genetics.118.301018). This repository conducts an individual based forward-time simulation in [SLiM](https://messerlab.org/slim/).

The parameters of the simulation model can be set by manipulating the parameters in `config.yaml` file inside `config` directory, and the details of the parameters are written in `README.md` file in the directory. The mathematical details of the simulation model is described in [Sachdeva and Barton (2018)](https://doi.org/10.1534/genetics.118.301018) as "individual-based simulation" (See section named "Individual-based simulations of long-term introgression into a finite population" in the paper for the details of the simulation).

1. Clone this repo

```
$ git clone https://github.com/daikitag/introgression-slim.git
```

2. Install `snakemake` by following the instructions written here: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html. If you have `coonda` installed on your computer, you can install `snakemake` by running:

```
$ conda install -c bioconda snakemake
```

3. Run the simulation command by

```
$ snakemake --cores 1 --use-conda
```

Here, `--cores` command specifies the number of cores to use and `--use-conda` allows `snakemake` to be run in a `conda` environment. This `snakemake` module uses a `conda` environment to run the simulations, so please do not forget to include this argument.

## Note:

Running this simulation code for the first time will require a lot of time, as it creates a new conda environment with the required Python packages for the simulation. However, it will not take a lot of time for the simulation afterwards, as conda environment is already constructed in the local repository.