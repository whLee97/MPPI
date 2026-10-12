# DXA application data

The application uses 2017-2018 NHANES target data and pooled 2017/2018 NHIS, BRFSS, and HINTS sources. Obtain the public files directly from their providers and place them in `data/raw_dxa/` with the following local names. Keep the NHIS, BRFSS, and HINTS archives compressed.

| Local filename | Provider and file |
| --- | --- |
| `DEMO_J.xpt` | [NHANES 2017-2018 demographics](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.htm) |
| `DXX_J.xpt` | [NHANES 2017-2018 DXA](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DXX_J.htm) |
| `WHQ_J.xpt` | [NHANES 2017-2018 weight history](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/WHQ_J.htm) |
| `BMX_J.xpt` | [NHANES 2017-2018 body measures](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BMX_J.htm) |
| `nhis_2017_samadultcsv.zip` | [NHIS 2017](https://ftp.cdc.gov/pub/health_statistics/nchs/datasets/NHIS/2017/), `samadultcsv.zip` |
| `nhis_2018_samadultcsv.zip` | [NHIS 2018](https://ftp.cdc.gov/pub/health_statistics/nchs/datasets/NHIS/2018/), `samadultcsv.zip` |
| `brfss_2017_xpt.zip` | [BRFSS 2017](https://www.cdc.gov/brfss/annual_data/annual_2017.html), SAS transport ZIP |
| `brfss_2018_xpt.zip` | [BRFSS 2018](https://www.cdc.gov/brfss/annual_data/annual_2018.html), SAS transport ZIP |
| `hints5_cycle1_2017_stata.zip` | [HINTS data downloads](https://hints.cancer.gov/data/download-data.aspx), HINTS 5 Cycle 1, Stata ZIP |
| `hints5_cycle2_2018_stata.zip` | [HINTS data downloads](https://hints.cancer.gov/data/download-data.aspx), HINTS 5 Cycle 2, Stata ZIP |

Follow each provider's download and use conditions, including any HINTS access acknowledgment. `input_checksums.json` records SHA-256 hashes of the exact raw files used for the reference analysis; revised provider files may have different hashes. `scripts/check_data.py` checks the supplied files against these hashes without downloading data.

```sh
python3 scripts/check_data.py
python3 -m pip install -r requirements.txt
python3 scripts/prepare_hints_2017_2018.py
Rscript --vanilla scripts/prepare_dxa.R
Rscript --vanilla scripts/run_dxa.R
```

The HINTS script writes `data/processed_dxa/source_hints_2017_2018.csv`. The R loader prepares NHANES/NHIS/BRFSS directly from the raw files and validates the following complete-case counts, restricting all groups to ages 40-59:

| Group | Target | NHIS | BRFSS | HINTS |
| --- | ---: | ---: | ---: | ---: |
| Male | 529 | 7,487 | 119,843 | 886 |
| Female | 584 | 8,320 | 133,665 | 1,337 |

The target outcome is total body-fat percentage. The target BMI uses self-reported height/weight (`WHD010`/`WHD020`) in these files; the loader retains its documented measured-BMI fallback if those variables are absent. The main thresholds are 30 for men and 40 for women; sensitivity thresholds are 28/32 and 38/42. Source assignments are NHIS/Deurenberg, BRFSS/CUN-BAE, and HINTS/Gallagher. Cross-fitting uses five folds and seed 20261005.

Survey records and processed records are ignored by Git and are not part of the public code archive. Their provider terms are separate from the repository's MIT license.
