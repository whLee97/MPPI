mppi_parameter_components <- function(data, nuisance, target) {
  x0 <- as.matrix(data$x0); y0 <- as.numeric(data$y0); xs <- lapply(data$xs, as.matrix)
  z0 <- mppi_design(target, x0)
  zs <- lapply(nuisance$aligned_sources, function(x) mppi_design(target, x))
  n0 <- nrow(z0); Ns <- vapply(zs, nrow, integer(1)); S <- length(zs)
  H0 <- crossprod(z0) / n0
  Hs <- lapply(seq_len(S), function(s) crossprod(zs[[s]], zs[[s]] * nuisance$ratios[[s]]) / Ns[s])
  h0 <- as.numeric(crossprod(z0, y0) / n0)
  hs <- lapply(seq_len(S), function(s) as.numeric(crossprod(zs[[s]], nuisance$ratios[[s]] * nuisance$pred_sources[[s]]) / Ns[s]))
  h0s <- lapply(seq_len(S), function(s) as.numeric(crossprod(z0, nuisance$pred_target[[s]]) / n0))
  list(z0 = z0, zs = zs, H0 = H0, Hs = Hs, h0 = h0, hs = hs, h0s = h0s,
       n0 = n0, Ns = Ns, p = ncol(z0), S = S)
}

mppi_solve_parameter <- function(components, w) {
  w <- as.numeric(w); K <- components$S + 1L
  mppi_assert(length(w) == K && all(w >= -1e-12) && abs(sum(w) - 1) <= 1e-8, "Invalid simplex weights.")
  H <- w[1L] * components$H0
  rhs <- components$h0
  for (s in seq_len(components$S)) {
    H <- H + w[s + 1L] * components$Hs[[s]]
    rhs <- rhs + w[s + 1L] * (components$hs[[s]] - components$h0s[[s]])
  }
  theta <- as.numeric(mppi_spd_solve(H, rhs, "weighted squared-loss Hessian"))
  score <- as.numeric(H %*% theta) - rhs
  score_norm <- sqrt(sum(score^2))
  score_tolerance <- 1 / components$n0
  mppi_assert(all(is.finite(theta)) && is.finite(score_norm) && score_norm <= score_tolerance,
              sprintf("Parameter score %.4e exceeds tolerance %.4e.", score_norm, score_tolerance))
  dg <- mppi_eigen_diagnostics(H)
  list(theta = theta, score_norm = score_norm, score_tolerance = score_tolerance,
       score_passed = TRUE, H = H,
       min_eigen_H = dg$min, condition_H = dg$condition)
}

mppi_build_Q <- function(theta, data, nuisance, components) {
  p <- components$p; S <- components$S; n0 <- components$n0
  G0 <- mppi_score(theta, components$z0, data$y0)
  delta <- source_scores <- vector("list", S)
  Dlist <- vector("list", S)
  for (s in seq_len(S)) {
    D <- nuisance$sensitivity(s, theta); Dlist[[s]] <- D
    pred0_score <- mppi_score(theta, components$z0, nuisance$pred_target[[s]])
    correction0 <- if (ncol(nuisance$influence_target[[s]]) == 0L) matrix(0, n0, p) else
      nuisance$influence_target[[s]] %*% t(D)
    delta[[s]] <- G0 - pred0_score + correction0
    source_base <- mppi_score(theta, components$zs[[s]], nuisance$pred_sources[[s]]) * nuisance$ratios[[s]]
    correction_s <- if (ncol(nuisance$influence_source[[s]]) == 0L) matrix(0, components$Ns[s], p) else
      nuisance$influence_source[[s]] %*% t(D)
    source_scores[[s]] <- source_base + correction_s
  }
  target_joint <- do.call(cbind, c(list(G0), delta))
  Q <- mppi_emp_cov(target_joint)
  for (s in seq_len(S)) {
    idx <- s * p + seq_len(p)
    Q[idx, idx] <- Q[idx, idx, drop = FALSE] + (n0 / components$Ns[s]) * mppi_emp_cov(source_scores[[s]])
  }
  Q <- mppi_sym(Q)
  qdiag <- mppi_eigen_diagnostics(Q)
  mppi_assert(qdiag$min >= -1e-8 * max(abs(qdiag$max), 1),
              sprintf("Full Q has a material negative eigenvalue %.4e.", qdiag$min))
  list(Q = Q, A = components$H0, target_joint_scores = target_joint,
       source_scores = source_scores, sensitivities = Dlist,
       diagnostics = list(min_eigen_Q = qdiag$min, max_eigen_Q = qdiag$max,
                          condition_Q = qdiag$condition, covariance_divisor = "N"))
}

mppi_covariance <- function(w, Q, A) {
  p <- nrow(A); K <- length(w)
  Ai <- mppi_spd_solve(A, diag(p), "target Hessian")
  W <- kronecker(matrix(w, ncol = 1L), diag(p))
  Sigma <- mppi_sym(Ai %*% crossprod(W, Q %*% W) %*% t(Ai))
  dg <- mppi_eigen_diagnostics(Sigma)
  mppi_assert(dg$min > 1e-11 * max(abs(dg$max), 1),
              sprintf("Final covariance is singular or indefinite (minimum eigenvalue %.4e).", dg$min))
  list(Sigma = Sigma, min_eigen_Sigma = dg$min, condition_Sigma = dg$condition)
}

mppi_inference <- function(theta, Sigma, n0, alpha = 0.05, truth = NULL) {
  p <- length(theta)
  if (p == 1L) {
    half <- qnorm(1 - alpha / 2) * sqrt(Sigma[1L, 1L] / n0)
    interval <- c(lower = theta - half, upper = theta + half)
    covered <- if (is.null(truth)) NA else as.integer(truth >= interval[1L] && truth <= interval[2L])
  } else {
    half <- sqrt(qchisq(1 - alpha, p) * diag(Sigma) / n0)
    interval <- cbind(lower = theta - half, upper = theta + half)
    if (is.null(truth)) covered <- NA else {
      diff <- theta - truth
      stat <- as.numeric(n0 * crossprod(diff, solve(Sigma, diff)))
      covered <- as.integer(stat <= qchisq(1 - alpha, p))
    }
  }
  list(interval = interval, covered = covered, volume = det(Sigma), se = sqrt(diag(Sigma) / n0))
}

mppi_fit_weight <- function(w, data, nuisance, components, alpha, truth, label) {
  point <- mppi_solve_parameter(components, w)
  qfinal <- mppi_build_Q(point$theta, data, nuisance, components)
  covariance <- mppi_covariance(w, qfinal$Q, qfinal$A)
  inference <- mppi_inference(point$theta, covariance$Sigma, components$n0, alpha, truth)
  list(method = label, w = as.numeric(w), theta = point$theta, Sigma = covariance$Sigma,
       interval = inference$interval, covered = inference$covered, volume = inference$volume,
       se = inference$se, diagnostics = c(point[c("score_norm", "score_tolerance", "score_passed", "min_eigen_H", "condition_H")],
         covariance[c("min_eigen_Sigma", "condition_Sigma")], qfinal$diagnostics))
}
