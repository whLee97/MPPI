# Validation of the public implementation

Checked on 2026-10-12 with R 4.5.1 and `foreign` 0.8-90 on macOS. Reference simulations were produced with R 4.5.3 on Linux.

- Core numerical tests: passed (nine groups).
- Framework tests: passed, including signatures, resume, failure and missing-task denominators.
- Formal design tests: passed for all twelve scenarios at the full sample sizes.
- Serial/parallel and recovery tests: passed, including deliberately interrupted workers and serialization failures.
- Linux launcher control-flow tests: passed using local fake Rscript/nproc executables.
- CLI validation and checkout discovery from another working directory, including paths with spaces and parentheses: passed.
- Clean source archive without survey data: reference rendering and simulation preflight passed.
- Reference data hashes, rosters, weights, and true-null rejection/coverage consistency: passed.
- End-to-end smoke run: two repetitions in each of twelve scenarios, 216 successful method records, using two workers.
- The same 216 method records agree with the approved formal runs at tolerance 1e-7 for estimates, covariances, weights, coverage, Vol, and PPI++ tuning parameters. This tolerance allows floating-point differences between R/platform versions.
- All ten raw survey input hashes match the recorded files. Reconstructed HINTS data agree with the accepted processed data, and all male/female complete-case counts match.
- All 54 DXA method/threshold/group results agree with the reference table at tolerance 1e-7.
- Reference rendering: eight single-page PDF figures and fourteen numeric CSV tables generated; figure pages rendered and visually checked.

This publication check verifies the consolidated implementation and reference exports. It does not claim a new 1,000-repetition rerun on macOS. The included references retain the approved 1,000-repetition results. The supplied GitHub Actions workflow has not yet run on GitHub.
