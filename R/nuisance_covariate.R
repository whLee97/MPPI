mppi_design <- function(target, x) {
  x <- as.matrix(x)
  if (target == "mean") matrix(1, nrow(x), 1L) else x
}

mppi_sensitivity <- function(component, theta) {
  p <- length(theta); q <- ncol(component$dc)
  out <- matrix(0, p, q)
  for (j in seq_len(q)) out[, j] <- component$dH[, , j, drop = TRUE] %*% theta - component$dc[, j]
  out
}

mppi_none_nuisance <- function(x0, xs, predictors, target) {
  x0 <- as.matrix(x0)
  list(
    type = "known",
    aligned_sources = lapply(xs, as.matrix),
    pred_target = lapply(predictors, mppi_predict, x = x0),
    pred_sources = Map(mppi_predict, predictors, xs),
    ratios = lapply(xs, function(x) rep(1, nrow(x))),
    influence_source = lapply(xs, function(x) matrix(0, nrow(x), 0L)),
    influence_target = lapply(xs, function(x) matrix(0, nrow(x0), 0L)),
    sensitivity_components = lapply(xs, function(x) NULL),
    sensitivity = function(s, theta) matrix(0, length(theta), 0L),
    diagnostics = list(type = "known", nuisance_uncertainty = "none")
  )
}

mppi_gaussian_dr_nuisance <- function(x0, xs, predictors, target, folds0, folds_sources, K = 5L) {
  x0 <- as.matrix(x0); xs <- lapply(xs, as.matrix)
  S <- length(xs); p <- ncol(mppi_design(target, x0)); d <- ncol(x0)
  ratios <- influence_source <- influence_target <- vector("list", S)
  components <- vector("list", S); diagnostics <- vector("list", S)
  pred0 <- lapply(predictors, mppi_predict, x = x0)
  preds <- Map(mppi_predict, predictors, xs)
  for (s in seq_len(S)) {
    n_s <- nrow(xs[[s]])
    ratio <- numeric(n_s); a <- matrix(NA_real_, n_s, d)
    grad_log_r <- matrix(NA_real_, n_s, d)
    mu_by_fold <- vector("list", K)
    for (k in seq_len(K)) {
      train <- folds_sources[[s]] != k; test <- !train
      mppi_assert(any(train) && any(test), "Empty Gaussian density-ratio fold.")
      mu <- colMeans(xs[[s]][train, , drop = FALSE])
      mu_by_fold[[k]] <- mu
      log_r <- as.numeric(-xs[[s]][test, , drop = FALSE] %*% mu + sum(mu^2) / 2)
      mppi_assert(all(is.finite(log_r)) && max(log_r) < log(.Machine$double.xmax),
                  "Nonfinite Gaussian density ratio; no clipping was applied.")
      ratio[test] <- exp(log_r)
      a[test, ] <- sweep(xs[[s]][test, , drop = FALSE], 2L, mu, "-")
      grad_log_r[test, ] <- matrix(mu, sum(test), d, byrow = TRUE) -
        xs[[s]][test, , drop = FALSE]
    }
    z <- mppi_design(target, xs[[s]])
    qgrad <- grad_log_r * ratio
    dH <- array(0, c(p, p, d)); dc <- matrix(0, p, d)
    for (j in seq_len(d)) {
      dH[, , j] <- crossprod(z, z * qgrad[, j]) / n_s
      dc[, j] <- as.numeric(crossprod(z, preds[[s]] * qgrad[, j]) / n_s)
    }
    ratios[[s]] <- ratio; influence_source[[s]] <- a
    influence_target[[s]] <- matrix(0, nrow(x0), d)
    components[[s]] <- list(dH = dH, dc = dc)
    diagnostics[[s]] <- list(
      source = s, fold_sizes = tabulate(folds_sources[[s]], K), target_fold_sizes = tabulate(folds0, K),
      fitted_means = mu_by_fold, ratio_min = min(ratio), ratio_max = max(ratio),
      ratio_mean = mean(ratio), effective_sample_size = sum(ratio)^2 / sum(ratio^2), clipped = FALSE
    )
  }
  list(type = "gaussian_density_ratio", aligned_sources = xs, pred_target = pred0,
       pred_sources = preds, ratios = ratios, influence_source = influence_source,
       influence_target = influence_target, sensitivity_components = components,
       sensitivity = function(s, theta) mppi_sensitivity(components[[s]], theta),
       folds_target = folds0, folds_sources = folds_sources, diagnostics = diagnostics)
}

mppi_logistic_features <- function(x) {
  x <- as.matrix(x)
  cbind(1, age = (x[, 1L] - 50) / 10, bmi = (x[, 2L] - 30) / 10)
}

mppi_fit_logistic_ratio <- function(x_source, x_target) {
  zs <- mppi_logistic_features(x_source); z0 <- mppi_logistic_features(x_target)
  z <- rbind(zs, z0); label <- c(rep(1, nrow(zs)), rep(0, nrow(z0)))
  fit <- glm.fit(z, label, family = binomial())
  mppi_assert(isTRUE(fit$converged) && fit$rank == ncol(z) && all(is.finite(fit$coefficients)),
              "Logistic density-ratio fit failed or is rank deficient.")
  gamma <- as.numeric(fit$coefficients)
  e <- plogis(as.numeric(z %*% gamma))
  H <- crossprod(z, z * (e * (1 - e))) / nrow(z)
  Hinv <- mppi_spd_solve(H, diag(ncol(z)), "logistic information")
  es <- e[seq_len(nrow(zs))]
  e0 <- e[nrow(zs) + seq_len(nrow(z0))]
  us <- zs * (1 - es)
  u0 <- -z0 * e0
  list(gamma = gamma, Hinv = Hinv, mean_us = colMeans(us), mean_u0 = colMeans(u0),
       pi_s = nrow(zs) / nrow(z), pi_0 = nrow(z0) / nrow(z),
       prior_ratio = nrow(zs) / nrow(z0))
}

mppi_logistic_dr_base <- function(x0, xs, folds0, folds_sources, K = 5L) {
  x0 <- as.matrix(x0); xs <- lapply(xs, as.matrix)
  S <- length(xs); q <- ncol(mppi_logistic_features(x0))
  ratios <- influence_source <- influence_target <- vector("list", S)
  grad_log_ratio <- vector("list", S); diagnostics <- vector("list", S)
  for (s in seq_len(S)) {
    n_s <- nrow(xs[[s]])
    ratio <- numeric(n_s); a <- matrix(NA_real_, n_s, q); b <- matrix(NA_real_, nrow(x0), q)
    grad_log_r <- matrix(NA_real_, n_s, q); fits <- vector("list", K)
    for (k in seq_len(K)) {
      tr_s <- folds_sources[[s]] != k; te_s <- !tr_s
      tr_0 <- folds0 != k; te_0 <- !tr_0
      fit <- mppi_fit_logistic_ratio(xs[[s]][tr_s, , drop = FALSE], x0[tr_0, , drop = FALSE])
      fits[[k]] <- fit
      zs <- mppi_logistic_features(xs[[s]][te_s, , drop = FALSE])
      z0 <- mppi_logistic_features(x0[te_0, , drop = FALSE])
      log_r <- log(fit$prior_ratio) - as.numeric(zs %*% fit$gamma)
      mppi_assert(all(is.finite(log_r)) && max(log_r) < log(.Machine$double.xmax),
                  "Nonfinite logistic density ratio; no clipping was applied.")
      ratio[te_s] <- exp(log_r)
      grad_log_r[te_s, ] <- -zs
      es <- plogis(as.numeric(zs %*% fit$gamma))
      e0 <- plogis(as.numeric(z0 %*% fit$gamma))
      us <- sweep(zs * (1 - es), 2L, fit$mean_us, "-")
      u0 <- sweep(-z0 * e0, 2L, fit$mean_u0, "-")
      a[te_s, ] <- fit$pi_s * t(fit$Hinv %*% t(us))
      b[te_0, ] <- fit$pi_0 * t(fit$Hinv %*% t(u0))
    }
    ratios[[s]] <- ratio; influence_source[[s]] <- a; influence_target[[s]] <- b
    grad_log_ratio[[s]] <- grad_log_r
    diagnostics[[s]] <- list(source = s, fold_sizes = tabulate(folds_sources[[s]], K),
      target_fold_sizes = tabulate(folds0, K), fitted = fits, ratio_min = min(ratio),
      ratio_max = max(ratio), ratio_mean = mean(ratio),
      effective_sample_size = sum(ratio)^2 / sum(ratio^2), clipped = FALSE,
      labels = c(source = 1L, target = 0L))
  }
  list(x0 = x0, xs = xs, ratios = ratios, influence_source = influence_source,
       influence_target = influence_target, grad_log_ratio = grad_log_ratio,
       folds_target = folds0, folds_sources = folds_sources, diagnostics = diagnostics)
}

mppi_logistic_dr_nuisance <- function(x0, xs, predictors, target, folds0 = NULL,
                                      folds_sources = NULL, K = 5L, base = NULL) {
  if (is.null(base)) base <- mppi_logistic_dr_base(x0, xs, folds0, folds_sources, K)
  x0 <- base$x0; xs <- base$xs; S <- length(xs)
  p <- ncol(mppi_design(target, x0)); q <- ncol(base$influence_source[[1L]])
  pred0 <- lapply(predictors, mppi_predict, x = x0)
  preds <- Map(mppi_predict, predictors, xs)
  components <- vector("list", S)
  for (s in seq_len(S)) {
    z <- mppi_design(target, xs[[s]])
    qgrad <- base$grad_log_ratio[[s]] * base$ratios[[s]]
    dH <- array(0, c(p, p, q)); dc <- matrix(0, p, q)
    for (j in seq_len(q)) {
      dH[, , j] <- crossprod(z, z * qgrad[, j]) / nrow(z)
      dc[, j] <- as.numeric(crossprod(z, preds[[s]] * qgrad[, j]) / nrow(z))
    }
    components[[s]] <- list(dH = dH, dc = dc)
  }
  list(type = "logistic_density_ratio", aligned_sources = xs, pred_target = pred0,
       pred_sources = preds, ratios = base$ratios, influence_source = base$influence_source,
       influence_target = base$influence_target, sensitivity_components = components,
       sensitivity = function(s, theta) mppi_sensitivity(components[[s]], theta),
       folds_target = base$folds_target, folds_sources = base$folds_sources,
       diagnostics = base$diagnostics)
}
