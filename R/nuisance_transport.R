mppi_pack_transport <- function(mu_s, sigma_s, mu_0, sigma_0) {
  c(as.numeric(mu_s), mppi_vech(sigma_s), as.numeric(mu_0), mppi_vech(sigma_0))
}

mppi_unpack_transport <- function(eta, d) {
  qv <- d * (d + 1L) / 2L
  i1 <- seq_len(d); i2 <- d + seq_len(qv); i3 <- d + qv + seq_len(d); i4 <- 2L * d + qv + seq_len(qv)
  list(mu_s = eta[i1], sigma_s = mppi_unvech(eta[i2], d),
       mu_0 = eta[i3], sigma_0 = mppi_unvech(eta[i4], d))
}

mppi_fit_transport_eta <- function(x_source, x_target) {
  x_source <- as.matrix(x_source); x_target <- as.matrix(x_target)
  mppi_pack_transport(colMeans(x_source), mppi_cov_n(x_source),
                      colMeans(x_target), mppi_cov_n(x_target))
}

mppi_transport_matrices <- function(eta, d, rel_tol = 1e-10) {
  z <- mppi_unpack_transport(eta, d)
  s0h <- mppi_matrix_power(z$sigma_0, 0.5, "target training covariance", rel_tol)
  s0i <- mppi_matrix_power(z$sigma_0, -0.5, "target training covariance", rel_tol)
  ssh <- mppi_matrix_power(z$sigma_s, 0.5, "source training covariance", rel_tol)
  ssi <- mppi_matrix_power(z$sigma_s, -0.5, "source training covariance", rel_tol)
  list(par = z, T = s0h %*% ssi, S = ssh %*% s0i)
}

mppi_transport_T <- function(x, eta, d, rel_tol = 1e-10) {
  mm <- mppi_transport_matrices(eta, d, rel_tol)
  sweep(sweep(as.matrix(x), 2L, mm$par$mu_s, "-") %*% t(mm$T), 2L, mm$par$mu_0, "+")
}

mppi_transport_S <- function(z, eta, d, rel_tol = 1e-10) {
  mm <- mppi_transport_matrices(eta, d, rel_tol)
  sweep(sweep(as.matrix(z), 2L, mm$par$mu_0, "-") %*% t(mm$S), 2L, mm$par$mu_s, "+")
}

mppi_transport_if <- function(x, eta, d, side = c("source", "target")) {
  side <- match.arg(side); x <- as.matrix(x); par <- mppi_unpack_transport(eta, d)
  qv <- d * (d + 1L) / 2L; out <- matrix(0, nrow(x), 2L * (d + qv))
  if (side == "source") {
    xc <- sweep(x, 2L, par$mu_s, "-")
    out[, seq_len(d)] <- xc
    for (i in seq_len(nrow(x))) out[i, d + seq_len(qv)] <- mppi_vech(tcrossprod(xc[i, ]) - par$sigma_s)
  } else {
    xc <- sweep(x, 2L, par$mu_0, "-")
    off <- d + qv
    out[, off + seq_len(d)] <- xc
    for (i in seq_len(nrow(x))) out[i, off + d + seq_len(qv)] <- mppi_vech(tcrossprod(xc[i, ]) - par$sigma_0)
  }
  out
}

mppi_transport_fold_moments <- function(eta, x_source, x_target, predictor, target,
                                        Ns_total, N0_total, rel_tol) {
  d <- ncol(x_target)
  t_source <- mppi_transport_T(x_source, eta, d, rel_tol)
  reconstructed <- mppi_transport_S(t_source, eta, d, rel_tol)
  source_scale_target <- mppi_transport_S(x_target, eta, d, rel_tol)
  fs <- mppi_predict(predictor, reconstructed)
  f0 <- mppi_predict(predictor, source_scale_target)
  zs <- mppi_design(target, t_source); z0 <- mppi_design(target, x_target)
  list(H = crossprod(zs) / Ns_total - crossprod(z0) / N0_total,
       c = as.numeric(crossprod(zs, fs) / Ns_total - crossprod(z0, f0) / N0_total),
       T = t_source, reconstructed = reconstructed, source_scale_target = source_scale_target,
       fs = fs, f0 = f0)
}

mppi_transport_derivative <- function(eta, j, fold_args, fd) {
  h <- max(fd$min_step, fd$relative_step * max(1, abs(eta[j])))
  for (attempt in seq_len(12L)) {
    ep <- em <- eta; ep[j] <- ep[j] + h; em[j] <- em[j] - h
    plus <- tryCatch(do.call(mppi_transport_fold_moments, c(list(eta = ep), fold_args)), error = identity)
    minus <- tryCatch(do.call(mppi_transport_fold_moments, c(list(eta = em), fold_args)), error = identity)
    if (!inherits(plus, "error") && !inherits(minus, "error")) {
      return(list(dH = (plus$H - minus$H) / (2 * h), dc = (plus$c - minus$c) / (2 * h),
                  step = h, attempts = attempt))
    }
    h <- h / 2
  }
  stop("Could not construct a positive-definite central difference for transport sensitivity.", call. = FALSE)
}

mppi_transport_nuisance <- function(x0, xs, predictors, target, folds0, folds_sources,
                                    K = 5L, finite_difference = list(relative_step = 2e-5,
                                    min_step = 1e-7, spd_rel_tol = 1e-10)) {
  x0 <- as.matrix(x0); xs <- lapply(xs, as.matrix)
  S <- length(xs); d <- ncol(x0); p <- ncol(mppi_design(target, x0)); qv <- d * (d + 1L) / 2L
  q <- 2L * (d + qv)
  aligned <- pred_sources <- pred_target <- influence_source <- influence_target <- vector("list", S)
  ratios <- lapply(xs, function(x) rep(1, nrow(x)))
  components <- vector("list", S); diagnostics <- vector("list", S)
  for (s in seq_len(S)) {
    n_s <- nrow(xs[[s]])
    aligned_s <- matrix(NA_real_, n_s, d); fs_s <- numeric(n_s); f0_s <- numeric(nrow(x0))
    a_s <- matrix(NA_real_, n_s, q); b_s <- matrix(NA_real_, nrow(x0), q)
    dH <- array(0, c(p, p, q)); dc <- matrix(0, p, q)
    eta_folds <- vector("list", K); step_folds <- vector("list", K)
    inverse_error <- numeric(K)
    for (k in seq_len(K)) {
      tr_s <- folds_sources[[s]] != k; te_s <- !tr_s
      tr_0 <- folds0 != k; te_0 <- !tr_0
      mppi_assert(any(tr_s) && any(te_s) && any(tr_0) && any(te_0), "Empty transport cross-fitting fold.")
      eta <- mppi_fit_transport_eta(xs[[s]][tr_s, , drop = FALSE], x0[tr_0, , drop = FALSE])
      eta_folds[[k]] <- eta
      args <- list(x_source = xs[[s]][te_s, , drop = FALSE], x_target = x0[te_0, , drop = FALSE],
                   predictor = predictors[[s]], target = target, Ns_total = n_s,
                   N0_total = nrow(x0), rel_tol = finite_difference$spd_rel_tol)
      base <- do.call(mppi_transport_fold_moments, c(list(eta = eta), args))
      aligned_s[te_s, ] <- base$T; fs_s[te_s] <- base$fs; f0_s[te_0] <- base$f0
      a_s[te_s, ] <- mppi_transport_if(xs[[s]][te_s, , drop = FALSE], eta, d, "source")
      b_s[te_0, ] <- mppi_transport_if(x0[te_0, , drop = FALSE], eta, d, "target")
      inverse_error[k] <- sqrt(mean(rowSums((base$reconstructed - xs[[s]][te_s, , drop = FALSE])^2)))
      step_k <- numeric(q)
      for (j in seq_len(q)) {
        der <- mppi_transport_derivative(eta, j, args, finite_difference)
        dH[, , j] <- dH[, , j] + der$dH
        dc[, j] <- dc[, j] + der$dc
        step_k[j] <- der$step
      }
      step_folds[[k]] <- step_k
    }
    mppi_assert(all(is.finite(aligned_s)) && all(is.finite(fs_s)) && all(is.finite(f0_s)) &&
                  all(is.finite(a_s)) && all(is.finite(b_s)), "Incomplete transport cross-fitting output.")
    aligned[[s]] <- aligned_s; pred_sources[[s]] <- fs_s; pred_target[[s]] <- f0_s
    influence_source[[s]] <- a_s; influence_target[[s]] <- b_s
    components[[s]] <- list(dH = dH, dc = dc)
    diagnostics[[s]] <- list(source = s, fold_sizes = tabulate(folds_sources[[s]], K),
      target_fold_sizes = tabulate(folds0, K), eta = eta_folds,
      finite_difference_steps = step_folds, inverse_rmse = inverse_error,
      covariance_divisor = "training sample size", ridge = FALSE, shrinkage = FALSE)
  }
  list(type = "affine_gaussian_transport", aligned_sources = aligned,
       pred_target = pred_target, pred_sources = pred_sources, ratios = ratios,
       influence_source = influence_source, influence_target = influence_target,
       sensitivity_components = components,
       sensitivity = function(s, theta) mppi_sensitivity(components[[s]], theta),
       folds_target = folds0, folds_sources = folds_sources, diagnostics = diagnostics)
}
