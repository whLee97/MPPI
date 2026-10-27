# Reproduction Guide

This document provides instructions for reproducing the simulation results, figures, and tables in the MPPI paper. All commands should be run from the repository directory.

## Reproduce Figures and Tables

To generate the figures and tables from the included reference results:

```sh
Rscript --vanilla scripts/reproduce.R
```

The outputs are saved in `outputs/reference/`. This step does not rerun the simulations or require external datasets.

## Run Simulations

To perform a small test run with two repetitions per scenario:

```sh
Rscript --vanilla scripts/run_all.R --B=2 --workers=1 --results=results/smoke
```

To run all twelve simulation settings with 1,000 repetitions each:

```sh
Rscript --vanilla scripts/run_all.R --B=1000 --workers=4 --results=results/formal
```

The full simulation uses the sample sizes, auxiliary predictors, cross-fitting procedures, and other settings specified in [Simulation Designs](designs.md).

### Individual Scenarios

To run a specific scenario or a subset of repetitions, use `scripts/run_simulations.R`.

For example, to run repetitions 1–100 of scenario 9:

```sh
Rscript --vanilla scripts/run_simulations.R --scenario=9 --B=1000 --from=1 --to=100 --workers=4 --results=results/formal
```

After completing the desired simulation runs, generate summaries and figures with:

```sh
Rscript --vanilla scripts/summarize_results.R --results=results/formal
Rscript --vanilla scripts/plot_results.R --results=results/formal
```

### Resume Simulations

Completed repetitions are saved automatically and can be reused when a run is restarted.

To resume an interrupted run, use the same command with the same `B`, code, configuration, runtime, and results directory.

If the code, configuration, or runtime changes, use a new results directory.

## Simulation Outputs

Simulation outputs are organized as follows:

- `scenario_*/rep_*.rds`: individual simulation repetitions.
- `summary/`: simulation summaries, rejection curves, and diagnostic results.
- `paper/figures/`: figures generated from the simulation results.
- `paper/tables/`: numerical results for the paper.

The results directories also contain information about completed, failed, and missing repetitions.

## Linux Server and Slurm

For execution on a Linux server with an allocated eight-core job:

```sh
MPPI_AVAILABLE_CORES=8 MPPI_RESERVE_WITHIN_ALLOCATION=1 MPPI_B=1000 \
  bash scripts/run_all_server.sh . results/formal
```

For execution through a Slurm scheduler:

```sh
MPPI_CPUS=8 MPPI_B=1000 MPPI_MEMORY=16G \
  bash scripts/submit_slurm.sh . results/formal
```

Adjust the CPU and memory settings according to the available computing resources.

## DXA Application

Instructions for downloading and preparing the NHANES, NHIS, BRFSS, and HINTS datasets are provided in [Data Preparation](../data/README.md).

The DXA application can be run independently of the simulations. Raw survey datasets are not included in this repository.

## Tests

Run the following commands to check the implementation:

```sh
Rscript --vanilla tests/run_tests.R
Rscript --vanilla tests/test_framework.R
Rscript --vanilla tests/test_designs.R
Rscript --vanilla tests/test_parallel.R
python3 tests/test_launcher.py
python3 tests/test_cli.py
python3 tests/test_reference.py
```

These tests cover the numerical implementation, simulation designs, parallel execution, command-line interfaces, and reference results.