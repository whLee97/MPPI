# Reproduction guide

## Two workflows

`scripts/reproduce.R` renders included reference results. `scripts/run_all.R` recomputes all twelve simulation settings and creates signed repetition checkpoints, CSV summaries, figures, and numeric tables. The two workflows write to distinct directories. Reference CSVs are immutable comparison outputs; they are never copied into a fresh run's checkpoints.

The full Monte Carlo workload contains 12,000 repetitions and 108,000 method fits. A `B=2` smoke run checks execution, not statistical power or coverage. Sample sizes, folds, predictors, and optimizer settings remain the formal values in smoke runs.

## Outputs

- `run_manifest.rds`: code/configuration/runtime contract and full task roster.
- `scenario_*/scenario_manifest.rds`: fixed auxiliary predictors and scenario signature.
- `scenario_*/rep_*.rds`: atomic, signed repetition checkpoints.
- `summary/replication_results.csv`: method-level estimates, covariances, weights, and statuses.
- `summary/simulation_summary.csv`: coverage, Monte Carlo standard error, Vol, weights, and expected/valid/failed/missing counts.
- `summary/mean_test_curves.csv`, `regression_test_curves.csv`: rejection probabilities from valid repetitions.
- `summary/normal_diagnostics.csv`: standardized errors for all available methods and coordinates.
- `paper/figures/`: paired linear/nonlinear panels for mean and regression in each shift regime, and MPPI mean normal diagnostics. A DXA p-value figure is added when DXA results exist.
- `paper/tables/`: numeric ACP/Vol and MPPI-weight tables by target/shift; DXA main and sensitivity tables when available.

Included reference normal diagnostics retain only the 6,000 MPPI mean standardized errors needed for the paper figure. Fresh runs retain all method/coordinate diagnostics. Reference summaries record successful results for all 108 setting-method combinations.

## Resume and failures

Use the same output directory, `B`, code, and runtime to resume. Existing completed checkpoints are reused. The runner does not silently merge results from different configurations or software versions. A changed run must use a new directory. Incomplete tasks remain missing, and failures retain their method rows and summary denominators. Investigate `summary/failures.csv` before interpreting a partial run. Output-directory locks prevent overlapping managers from writing the same scenario/run; inspect lock ownership before removing a stale lock after a terminated job.

## Linux server and Slurm

For direct execution inside an explicit 8-core allocation, reserving one core within that allocation:

```sh
MPPI_AVAILABLE_CORES=8 MPPI_RESERVE_WITHIN_ALLOCATION=1 MPPI_B=1000 \
  bash scripts/run_all_server.sh . results/formal
```

The launcher checks the visible and Slurm-allocated CPU limits, sets BLAS/OpenMP threads to one, runs all twelve scenarios, and writes per-step logs. Simulations do not require survey data. Do not launch another manager against the same results directory.

For Slurm, put `Rscript` on PATH and specify the resources appropriate for your cluster:

```sh
MPPI_CPUS=8 MPPI_B=1000 MPPI_PARTITION=your_partition MPPI_MEMORY=16G \
  bash scripts/submit_slurm.sh . results/formal
```

`MPPI_NODE` is optional. The submission script uses the scheduler's default partition if `MPPI_PARTITION` is omitted. The DXA analysis is a separate command after data preparation.

## Validation

The tests check covariance construction, analytic/finite-difference gradients, score acceptance, affine transport, deterministic seeds, resume signatures, missing/failed denominators, complementary DGP moments, Gaussian Hermite orthogonality, all twelve designs at full sample sizes, serial/parallel equivalence, and launcher failure/lock handling. GitHub Actions runs checks on Linux and macOS. Real-data fitting is run separately because survey files are external inputs.

R 4.5.3 generated the reference results; R 4.5.1 with `foreign` 0.8-90 is used for the local publication check. The Python preprocessing dependencies are pinned in `requirements.txt`. Floating-point differences across operating systems, BLAS implementations, and R versions can affect the last digits and optimizer ties; compare with tolerances and report the runtime recorded in each manifest.
