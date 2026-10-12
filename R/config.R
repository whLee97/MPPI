mppi_simulation_config <- function(B = 1000L, K = 5L) {
  scenarios <- expand.grid(
    shift = c("homogeneous", "covariate", "domain"),
    target = c("mean", "regression"),
    dgp = c("linear", "nonlinear"),
    stringsAsFactors = FALSE
  )
  scenarios$id <- seq_len(nrow(scenarios))
  scenarios <- scenarios[, c("id", "shift", "target", "dgp")]
  list(
    version = "mppi-paper-v1",
    B = as.integer(B), K = as.integer(K), alpha = 0.05,
    regression_quadratic = c(a = 1, b = 1),
    mean_nonlinear = list(noise_sd = 1),
    S = 3L, d = 2L, N0 = 5000L, Ns = c(10000L, 20000L, 30000L),
    beta = c(1.2, -0.8),
    covariate_means = list(rep(-0.5, 2), rep(-1, 2), rep(1.5, 2)),
    domain_A = lapply(seq_len(3L), function(s) diag(1 + 0.05 * s * c(0.5, 1))),
    domain_b = lapply(seq_len(3L), function(s) s * c(0.05, -0.025)),
    scenarios = scenarios,
    seeds = list(auxiliary = 123L, inference = 12001L, folds = 27011L, candidates = 41017L),
    optimizer = list(M_n = NULL, max_iter = 5000L, eps_w = 1e-7,
                     eta0 = 1, rho = 0.5, armijo = 0.1, max_backtrack = 60L),
    finite_difference = list(relative_step = 2e-5, min_step = 1e-7, spd_rel_tol = 1e-10)
  )
}

mppi_validate_config <- function(cfg) {
  mppi_assert(length(cfg$B) == 1L && is.finite(cfg$B) && cfg$B >= 1L, "B must be a positive integer.")
  mppi_assert(cfg$S == 3L && cfg$d == 2L, "The confirmed simulations require S=3 and d=2.")
  mppi_assert(cfg$N0 == 5000L && identical(as.integer(cfg$Ns), c(10000L, 20000L, 30000L)),
              "The confirmed sample sizes have been changed.")
  mppi_assert(cfg$K == 5L, "The confirmed simulations use five-fold cross-fitting.")
  invisible(cfg)
}

mppi_g_standard <- function(x, dgp) {
  x <- as.matrix(x)
  if (dgp == "linear") return(1.2 * x[, 1L] - 0.8 * x[, 2L])
  if (dgp == "nonlinear") return(2 * x[, 1L] + sin(2 * pi * x[, 1L]) - 0.8 * x[, 2L])
  stop("Unknown standard DGP.", call. = FALSE)
}

mppi_h3 <- function(z) (z^3 - 3*z) / sqrt(6)

mppi_g_domain_mean <- function(z, dgp) {
  z <- as.matrix(z)
  if (dgp == "linear") return(1 + 1.2 * z[, 1L] + 1.2 * z[, 2L])
  if (dgp == "nonlinear") return(1 + mppi_h3(z[, 1L]) + mppi_h3(z[, 2L]))
  stop("Unknown domain-shift mean DGP.", call. = FALSE)
}

mppi_g_complementary_regression <- function(x, dgp, coefficients = c(a = 1, b = 1)) {
  x <- as.matrix(x)
  linear <- 1.2 * x[, 1L] - 0.8 * x[, 2L]
  if (dgp == "linear") return(linear)
  if (dgp == "nonlinear")
    return(linear + coefficients[1L] * (x[, 1L]^2 - 1) + coefficients[2L] * (x[, 2L]^2 - 1))
  stop("Unknown regression DGP.", call. = FALSE)
}


mppi_domain_partial_signal <- function(z, source, dgp) {
  z <- as.matrix(z)
  if (dgp == "linear") {
    if (source == 1L) return(1.2 * z[, 1L])
    if (source == 2L) return(-0.8 * z[, 2L])
    return(1.2 * z[, 1L] - 0.8 * z[, 2L])
  }
  stop("Partial-response auxiliary models are defined only for linear regression.", call. = FALSE)
}

mppi_truth <- function(target, shift, dgp) {
  if (target == "mean") return(if (shift == "domain") 1 else 0)
  c(1.2, -0.8)
}
