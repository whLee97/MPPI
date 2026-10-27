# Source-specific nonlinear features are evaluated on source-scale covariates.
# OLS includes an intercept to represent source-specific mean shifts.
mppi_feature_complementary <- function(x, source) {
  x <- as.matrix(x)
  linear <- cbind(1, x[, 1L], x[, 2L])
  if (source == 1L) return(cbind(linear, x[, 1L]^2 - 1))
  if (source == 2L) return(cbind(linear, x[, 2L]^2 - 1))
  if (source == 3L) return(linear)
  stop("Unknown complementary predictor source.", call. = FALSE)
}

mppi_feature_domain <- function(x, target, dgp, source, preprocessing = NULL, training = FALSE) {
  x <- as.matrix(x)
  if (target == "mean") {
    if (is.null(preprocessing)) {
      mppi_assert(training, "Mean-predictor preprocessing must be learned from its auxiliary training sample.")
      preprocessing <- list(center = colMeans(x), scale = pmax(apply(x, 2L, sd), 1e-8))
    }
    u <- sweep(sweep(x, 2L, preprocessing$center, "-"), 2L, preprocessing$scale, "/")
    if (source == 1L) {
      out <- if (dgp == "linear") {
        cbind(1, u[, 1L]^3)
      } else {
        cbind(1, u[, 1L], mppi_h3(u[, 1L]))
      }
    } else if (source == 2L) {
      out <- if (dgp == "linear") {
        cbind(1, u[, 2L]^3)
      } else {
        cbind(1, u[, 2L], mppi_h3(u[, 2L]))
      }
    } else {
      out <- matrix(1, nrow(x), 1L)
    }
    return(list(x = out, preprocessing = preprocessing))
  }
  mppi_assert(dgp == "linear", "Nonlinear regression uses complementary predictor features.")
  out <- cbind(1, x[, 1L], x[, 2L])
  list(x = out, preprocessing = list())
}

mppi_fit_predictor <- function(x, y, kind = c("standard", "domain", "complementary"), target = NULL,
                               dgp = NULL, source = NULL) {
  kind <- match.arg(kind)
  x <- as.matrix(x); y <- as.numeric(y)
  if (kind == "complementary") {
    design <- mppi_feature_complementary(x, source)
    prep <- list()
  } else if (kind == "standard") {
    design <- x
    prep <- list()
  } else {
    f <- mppi_feature_domain(x, target, dgp, source, training = TRUE)
    design <- f$x; prep <- f$preprocessing
  }
  fit <- lm.fit(design, y)
  mppi_assert(fit$rank == ncol(design) && all(is.finite(fit$coefficients)), "Auxiliary OLS fit is rank deficient.")
  structure(list(coefficients = as.numeric(fit$coefficients), kind = kind, target = target,
                 dgp = dgp, source = source, preprocessing = prep), class = "mppi_predictor")
}

mppi_predict <- function(object, x) {
  mppi_assert(inherits(object, "mppi_predictor"), "Invalid MPPI predictor object.")
  x <- as.matrix(x)
  if (object$kind == "dxa_threshold") {
    age <- x[, 1L]; bmi <- x[, 2L]; male <- object$male
    continuous <- switch(object$formula,
      Deurenberg = 1.20 * bmi + 0.23 * age - 10.8 * male - 5.4,
      `CUN-BAE` = -44.988 + 0.503 * age + 10.689 * (1 - male) + 3.172 * bmi -
        0.026 * bmi^2 + 0.181 * bmi * (1 - male) - 0.020 * bmi * age -
        0.005 * bmi^2 * (1 - male) + 0.00021 * bmi^2 * age,
      Gallagher = 64.5 - 848 / bmi + 0.079 * age - 16.4 * male +
        0.05 * male * age + 39 * male / bmi,
      stop("Unknown DXA prediction equation.", call. = FALSE))
    return(as.numeric(continuous > object$threshold))
  }
  design <- if (object$kind == "complementary") mppi_feature_complementary(x, object$source) else
    if (object$kind == "standard") x else
    mppi_feature_domain(x, object$target, object$dgp, object$source,
                        preprocessing = object$preprocessing, training = FALSE)$x
  as.numeric(design %*% object$coefficients)
}

mppi_dxa_predictor <- function(formula, threshold, male) {
  structure(list(kind = "dxa_threshold", formula = formula, threshold = threshold,
                 male = male), class = "mppi_predictor")
}

mppi_fit_auxiliary_predictors <- function(cfg, shift, target, dgp, scenario_id) {
  auxiliary_key <- if (shift == "domain") {
    match(paste(target, dgp), c("mean linear", "regression linear", "mean nonlinear", "regression nonlinear")) + 4L
  } else {
    match(paste(shift, dgp), c("homogeneous linear", "covariate linear",
                              "homogeneous nonlinear", "covariate nonlinear"))
  }
  mppi_assert(!is.na(auxiliary_key), "Unknown auxiliary-predictor configuration.")
  if (target == "regression" && dgp == "nonlinear" && shift != "domain")
    auxiliary_key <- 10L + match(shift, c("homogeneous", "covariate"))
  seed <- mppi_seed(cfg$seeds$auxiliary, auxiliary_key, 0L, 1L)
  mppi_with_seed(seed, {
    models <- vector("list", cfg$S)
    for (s in seq_len(cfg$S)) {
      if (shift == "domain") {
        z <- mppi_rstd_logistic(cfg$Ns[s], cfg$d)
        x <- sweep(z %*% t(cfg$domain_A[[s]]), 2L, cfg$domain_b[[s]], "+")
        signal <- if (target == "mean") mppi_g_domain_mean(z, dgp) else
          if (dgp == "nonlinear") mppi_g_domain_regression(z, dgp, cfg$regression_quadratic) else
          mppi_domain_partial_signal(z, s, dgp)
        noise_sd <- if (target != "mean") 1 else if (dgp == "nonlinear")
          cfg$mean_nonlinear$noise_sd else cfg$domain_mean_linear_noise_sd
        y <- signal + rnorm(cfg$Ns[s], sd = noise_sd)
        kind <- if (target == "regression" && dgp == "nonlinear") "complementary" else "domain"
        models[[s]] <- mppi_fit_predictor(x, y, kind, target, dgp, s)
      } else {
        mu <- if (shift == "covariate") cfg$covariate_means[[s]] else rep(0, cfg$d)
        x <- mppi_rmvnorm(cfg$Ns[s], mu, diag(cfg$d))
        signal <- if (target == "regression") mppi_g_complementary_regression(x, dgp, cfg$regression_quadratic) else mppi_g_standard(x, dgp)
        y <- signal + rnorm(cfg$Ns[s])
        kind <- if (target == "regression" && dgp == "nonlinear") "complementary" else "standard"
        models[[s]] <- mppi_fit_predictor(x, y, kind, target, dgp, s)
      }
    }
    models
  })
}
