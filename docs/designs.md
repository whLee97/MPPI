# Simulation specification

The fixed formal configuration is defined in `R/config.R`. There are three sources, two covariates, target size 5,000, source sizes 10,000/20,000/30,000, five folds, 1,000 repetitions, and nominal level 0.05. Target covariates are standard bivariate Gaussian with independent coordinates. Every source predictor is fitted once using an independent labeled auxiliary sample of the corresponding source size, then held fixed throughout inference.

## Scenario roster

| ID | Shift | Target | DGP |
| --- | --- | --- | --- |
| 1 | Homogeneous | Mean | Linear |
| 2 | Covariate | Mean | Linear |
| 3 | Domain | Mean | Linear |
| 4 | Homogeneous | Regression | Linear |
| 5 | Covariate | Regression | Linear |
| 6 | Domain | Regression | Linear |
| 7 | Homogeneous | Mean | Nonlinear |
| 8 | Covariate | Mean | Nonlinear |
| 9 | Domain | Mean | Nonlinear |
| 10 | Homogeneous | Regression | Nonlinear |
| 11 | Covariate | Regression | Nonlinear |
| 12 | Domain | Regression | Nonlinear |

Homogeneous sources have the target covariate law. Covariate-shift sources have identity covariance and means (-0.5,-0.5), (-1,-1), and (1.5,1.5). For domain shift, source covariates are `X_s = A_s Z + b_s`, with `Z ~ N(0,I_2)`, `A_s = diag(1 + 0.05 s (0.5,1))`, and `b_s = s (0.05,-0.025)`.

## Mean inference

In homogeneous and covariate-shift settings, the linear signal is `1.2 X_1 - 0.8 X_2`; the nonlinear signal is `2 X_1 + sin(2 pi X_1) - 0.8 X_2`. Independent noise is standard normal. The target mean is zero. Auxiliary predictors are linear OLS without an intercept using both source covariates.

In the linear domain-shift setting, the signal is `1 + 1.2 Z_1 + 1.2 Z_2`, with independent normal noise of standard deviation 7. Source 1 uses `(1,u_1)`, source 2 uses `(1,u_2)`, and source 3 uses an intercept. Here each `u` is standardized using the corresponding independent auxiliary sample; preprocessing is frozen before inference. The target mean is one.

In the nonlinear domain-shift setting, the signal is `1 + h_3(Z_1) + h_3(Z_2)`, where `h_3(t) = (t^3 - 3t)/sqrt(6)`, with independent standard normal noise. Each source auxiliary response follows the full signal. Source 1 uses `(1,u_1,h_3(u_1))`, source 2 uses `(1,u_2,h_3(u_2))`, and source 3 uses an intercept. The target mean remains one. Full Gaussian affine alignment, five-fold nuisance estimation, and the calibration contributions in the estimated covariance are retained.

## Regression

The inferential target is the two-dimensional population least-squares coefficient without an intercept. Linear responses have signal `1.2 X_1 - 0.8 X_2` and standard normal noise. For domain shift this response is specified in latent target coordinates `Z`; linear source auxiliary responses use `1.2 Z_1`, `-0.8 Z_2`, and `1.2 Z_1 - 0.8 Z_2`, respectively. Auxiliary linear predictors use both observed covariates, with an intercept in the domain-shift setting.

All nonlinear regression settings have full response signal `1.2 Z_1 - 0.8 Z_2 + (Z_1^2-1) + (Z_2^2-1)` and standard normal noise. For homogeneous/covariate-shift settings, replace `Z` by `X` in the auxiliary response law. Source 1 uses `(1,X_1,X_2,X_1^2-1)`, source 2 uses `(1,X_1,X_2,X_2^2-1)`, and source 3 uses `(1,X_1,X_2)` on its observed scale. The target working model remains linear without an intercept. Gaussian orthogonality gives the same target coefficient `(1.2,-0.8)` in all regression settings.

## Comparisons and diagnostics

Every method uses the same auxiliary predictors and inference samples in a given repetition. MPPI estimates simplex weights; Classic uses only target labels; EW uses equal weights over the target and three sources; source-specific PPI uses a source endpoint. PPI++ uses the existing variance/trace tuning implementation, with heterogeneous covariance corrections.

Summary `mean_volume` is the average estimated asymptotic covariance determinant (the scalar covariance in mean inference), following the paper's Vol convention. It is not the geometric length/area of the confidence region. Coverage uses the corresponding Wald interval or joint region. Mean rejection curves evaluate candidate scalar nulls; regression rejection curves average over 100 evenly spaced directions at each positive radius and use the single true null at radius zero. The level is 0.05. The domain-mean display uses the original [0.7,1.3] range, with additional local null values in [0.9,1.1] for line resolution.

Base seeds are auxiliary 123, inference 12001, folds 27011, and optimization candidates 41017. Scenario/repetition/stage seed derivation is deterministic and independent of worker scheduling. Auxiliary seeds are intentionally shared in the unchanged designs where the predictor training law is shared; nonlinear regression has its own fixed auxiliary keys. See `mppi_seed` and `mppi_fit_auxiliary_predictors` for the exact mapping.
