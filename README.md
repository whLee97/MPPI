# Multi-Source Prediction-Powered Inference (MPPI)

This repository contains the code for **Multi-Source Prediction-Powered Inference**, by Wenhui Li, Fen Jiang, and Xinyu Zhang. It includes the proposed method, competing methods, twelve simulation settings, and the DXA application.

## Requirements

- R 4.5.x with the packages `foreign`, `ggplot2`, and `cowplot`.
- Python 3.11 or later for HINTS data preparation.
- Linux or macOS is recommended for parallel execution. Use `--workers=1` on Windows.

Install the required R packages with:

```sh
Rscript -e 'install.packages(c("foreign", "ggplot2", "cowplot"), repos="https://cloud.r-project.org")'
```

## Reproduce Figures and Tables

From the repository directory, run:

```sh
Rscript --vanilla scripts/reproduce.R
```

This generates the paper's figures and tables using the results provided in `reference/`. The outputs are saved in `outputs/reference/`. No simulation runs or external datasets are required.

## Run Simulations

To run all twelve simulation settings with 1,000 repetitions each:

```sh
Rscript --vanilla scripts/run_all.R --B=1000 --workers=4 --results=results/formal
```

The simulation results, figures, and tables are saved under `results/formal/`.

For individual scenarios and additional options, see [Reproduction Guide](docs/reproduction.md).

## DXA Application

Download the NHANES, NHIS, BRFSS, and HINTS datasets following the instructions in [Data Preparation](data/README.md), and place the required files in `data/raw_dxa/`.

Then run:

```sh
python3 -m pip install -r requirements.txt
python3 scripts/prepare_hints_2017_2018.py
Rscript --vanilla scripts/prepare_dxa.R
Rscript --vanilla scripts/run_dxa.R --results=results/formal/dxa
```

The DXA results are saved in `results/formal/dxa/`.

After completing both the simulations and the DXA analysis, generate figures and tables from the newly computed results with:

```sh
Rscript --vanilla scripts/reproduce.R --source=run --results=results/formal --output=outputs/formal
```

Raw survey data are not included in this repository.

## Documentation

- [Simulation Designs](docs/designs.md)
- [Reproduction Guide](docs/reproduction.md)
- [Data Preparation](data/README.md)

## License and Citation

This repository is released under the [MIT License](LICENSE). Please cite the accompanying paper when using the code or results. Citation information is available in [CITATION.cff](CITATION.cff).