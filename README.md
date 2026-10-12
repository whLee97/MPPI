# Multi-Source Prediction-Powered Inference (MPPI)

Reproduction code for **Multi-Source Prediction-Powered Inference**, by Wenhui Li, Fen Jiang, and Xinyu Zhang. The repository implements MPPI and its comparison methods, the paper's twelve simulation settings, and the DXA application.

## Requirements

- R 4.5.x; simulations and figures use base/recommended R packages only.
- R package `foreign` for the DXA application (tested with 0.8-90).
- Python 3.11 or later and `requirements.txt` only for HINTS preparation.
- Linux/macOS for parallel workers; use `--workers=1` on Windows.

The reference runs used R 4.5.3 on Linux. The public implementation is also checked locally with R 4.5.1. Numerical tolerances may be needed across R/BLAS/platform versions.

## Quick start: reproduce the reported figures and tables

From the repository directory:

```sh
Rscript --vanilla scripts/reproduce.R
```

This reads the included `reference/` CSVs and writes eight figures and fourteen numeric tables under `outputs/reference/`. It requires no survey downloads and does not rerun the Monte Carlo experiments. The reference results contain 1,000 repetitions in each of the twelve formal scenarios; `reference/provenance.json` records the scenario groups and checksums.

## Rerun the simulations

```sh
# Smoke run, using the paper's sample sizes in every setting:
Rscript --vanilla scripts/run_all.R --B=2 --workers=1 --results=results/smoke

# Formal run: all 12 settings, 1,000 repetitions per setting:
Rscript --vanilla scripts/run_all.R --B=1000 --workers=4 --results=results/formal
```

The default is all twelve scenarios, with target size 5,000, source sizes 10,000/20,000/30,000, five cross-fitting folds, and nominal level 0.05. Each scenario compares MPPI, Classic, EW, three source-specific PPI estimators, and three source-specific PPI++ estimators.

For a single scenario or a range of repetitions:

```sh
Rscript --vanilla scripts/run_simulations.R --scenario=9 --B=1000 --from=1 --to=100 --workers=4 --results=results/formal
Rscript --vanilla scripts/summarize_results.R --results=results/formal
Rscript --vanilla scripts/plot_results.R --results=results/formal
```

Use the same `B` when resuming. Completed repetitions are reused only when the code, configuration, fitted auxiliary predictors, and runtime signatures agree. Changing the code/configuration/runtime requires a new results directory. Missing and failed repetitions remain visible in the summary denominators; exit status 2 indicates statistical failures and status 3 from the single-scenario runner indicates infrastructure failure. See [the reproduction guide](docs/reproduction.md).

## DXA application

Obtain the ten source files listed in [the data guide](data/README.md), then run:

```sh
python3 -m pip install -r requirements.txt
python3 scripts/prepare_hints_2017_2018.py
Rscript --vanilla scripts/prepare_dxa.R
Rscript --vanilla scripts/run_dxa.R --results=results/formal/dxa
Rscript --vanilla scripts/reproduce.R --source=run --results=results/formal --output=outputs/formal
```

The data validation step checks all accepted male/female sample counts. Raw and processed survey records are local inputs and are excluded from Git. Simulation runs do not depend on these files.

## Repository layout

| Path | Contents |
| --- | --- |
| `R/` | DGPs, fixed auxiliary predictors, nuisance estimation, inference, optimization, summaries, and paper outputs |
| `config/` | Formal simulation configuration |
| `scripts/` | Simulation, data preparation, analysis, plotting, and optional Slurm entry points |
| `tests/` | Numerical, framework, design, parallel execution, and launcher checks |
| `reference/` | Compact reference summaries, rejection curves, MPPI normal diagnostics, and DXA results |
| `docs/` | Design specification and reproduction instructions |
| `data/` | Download instructions and input checksums |

The nonlinear regression settings use complementary quadratic features. Nonlinear mean inference under domain shift uses complementary cubic Hermite features. The affine alignment model and inference corrections are shared with the other domain-shift settings. The full scenario roster and predictor definitions are documented in [the design specification](docs/designs.md).

## Tests

```sh
Rscript --vanilla tests/run_tests.R
Rscript --vanilla tests/test_framework.R
Rscript --vanilla tests/test_designs.R
Rscript --vanilla tests/test_parallel.R
python3 tests/test_launcher.py
python3 tests/test_cli.py
python3 tests/test_reference.py
```

GitHub Actions is configured to run these checks and rebuild reference outputs after pushes and pull requests.

## License and citation

Code and accompanying repository documentation are released under the [MIT License](LICENSE). External survey data remain subject to the terms of their original providers. See [CITATION.cff](CITATION.cff) for authorship and citation metadata. Please cite the accompanying paper when using this implementation or the reference results in research.
