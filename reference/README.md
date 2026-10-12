# Reference results

These compact CSVs contain the paper's approved results: 1,000 successful repetitions per setting, twelve settings, and nine methods. They combine the formal runs for the unchanged mean settings, the complementary regression settings, and domain-shift mean inference. `provenance.json` records source-summary hashes, scenario groups, sample sizes, and reference-file SHA-256 hashes.

`simulation_summary.csv` retains unrounded coverage, estimated covariance determinant (Vol), weights, Monte Carlo standard errors, and success/failure denominators. Rejection curves are computed from the same repetition estimates and covariance matrices. Domain-mean curves include additional local nulls within the original display range. `normal_diagnostics.csv` retains the 6,000 MPPI mean standardized errors needed for the six-panel figure; fresh simulation summaries include all method/coordinate diagnostics. DXA results include the main and sensitivity thresholds; its p-value curves use the main thresholds.

Run `Rscript --vanilla scripts/reproduce.R` to render figures and export numeric tables. Reference results are comparison data, not resumable checkpoints. Fresh results are stored in a separate signed run directory. Small platform/runtime differences can affect the last digits of a rerun.
