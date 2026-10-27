# Reference Results

This directory contains the reference results for the twelve simulation scenarios and the DXA application.

The main files are:

- `simulation_summary.csv`: coverage probabilities, covariance volumes, and MPPI weights.
- `mean_test_curves.csv`: rejection probabilities for mean inference.
- `regression_test_curves.csv`: rejection probabilities for regression inference.
- `normal_diagnostics.csv`: standardized errors for MPPI mean inference.
- `dxa_results.csv`: DXA application estimates and confidence intervals.
- `dxa_pvalue_curves.csv`: DXA application p-value curves.
- `figures/`: selected figures from the paper.

The simulation results are based on 1,000 replications per scenario. Random seeds and SHA-256 checksums are recorded in `provenance.json`.

To regenerate the figures and tables, run the following command from the repository directory:

```sh
Rscript --vanilla scripts/reproduce.R
```