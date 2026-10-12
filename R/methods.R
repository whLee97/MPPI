mppi_ppi_plus <- function(source, data, nuisance, components, alpha, truth) {
  S <- components$S; n0 <- components$n0; Ns <- components$Ns
  lambda_pre <- Ns[source] / (n0 + Ns[source])
  w_pre <- numeric(S + 1L); w_pre[1L] <- 1 - lambda_pre; w_pre[source + 1L] <- lambda_pre
  preliminary <- mppi_solve_parameter(components, w_pre)
  qpre <- mppi_build_Q(preliminary$theta, data, nuisance, components)
  p <- components$p; Ai <- mppi_spd_solve(qpre$A, diag(p), "PPI++ target Hessian")
  block <- function(j, k) {
    ij <- (j - 1L) * p + seq_len(p); ik <- (k - 1L) * p + seq_len(p)
    Ai %*% qpre$Q[ij, ik, drop = FALSE] %*% t(Ai)
  }
  V00 <- block(1L, 1L); Vss <- block(source + 1L, source + 1L)
  V0s <- block(1L, source + 1L); Vs0 <- block(source + 1L, 1L)
  a <- sum(diag(V00)); b <- sum(diag(Vss)); c0 <- 0.5 * sum(diag(V0s + Vs0))
  denominator <- a + b - 2 * c0
  flat_tol <- 1e-12 * max(1, abs(a), abs(b), abs(c0))
  if (denominator <= flat_tol) {
    lambda <- 0; status <- "flat_or_nonconvex_trace_target_target_tie_rule"
  } else {
    lambda <- min(1, max(0, (a - c0) / denominator)); status <- "analytic_trace_minimizer"
  }
  w <- numeric(S + 1L); w[1L] <- 1 - lambda; w[source + 1L] <- lambda
  fit <- mppi_fit_weight(w, data, nuisance, components, alpha, truth,
                         paste0("PPI++ source ", source))
  fit$lambda <- lambda
  fit$diagnostics$ppi_plus <- list(source = source, lambda_pre = lambda_pre, lambda = lambda,
                                  criterion = if (p == 1L) "variance" else "trace",
                                  a = a, b = b, c = c0, denominator = denominator, status = status)
  fit
}


mppi_method_labels <- function(S = 3L) {
  c("MPPI", "Classic", paste0("PPI source ", seq_len(S)), "EW", paste0("PPI++ source ", seq_len(S)))
}

mppi_attempt_method <- function(label, code) {
  tryCatch({ fit <- force(code); fit$method <- label; fit$status <- "success"; fit$message <- ""; fit },
    error = function(e) {
      mppi_catch_statistical_error(e)
      list(method = label, status = "failed", message = conditionMessage(e))
    })
}

# Each source can fail independently during nuisance construction.
mppi_combine_nuisance <- function(parts) {
  fields <- c("aligned_sources", "pred_target", "pred_sources", "ratios", "influence_source", "influence_target")
  out <- setNames(lapply(fields, function(f) lapply(parts, function(x) {
    if (inherits(x, "error")) NULL else x[[f]][[1L]]
  })), fields)
  out$source_errors <- lapply(parts, function(x) if (inherits(x, "error")) conditionMessage(x) else NULL)
  out$diagnostics <- lapply(parts, function(x) if (inherits(x, "error")) list(error = conditionMessage(x)) else if (identical(x$type, "known")) x$diagnostics else x$diagnostics[[1L]])
  out$sensitivity <- function(s, theta) parts[[s]]$sensitivity(1L, theta)
  out
}

mppi_method_context <- function(data, nuisance, target, sources) {
  if (inherits(nuisance, "error")) stop(conditionMessage(nuisance), call. = FALSE)
  for (s in sources) if (!is.null(nuisance$source_errors[[s]]))
    stop(sprintf("Source %d nuisance failed: %s", s, nuisance$source_errors[[s]]), call. = FALSE)
  dat <- data; dat$xs <- data$xs[sources]; nu <- nuisance
  for (f in c("aligned_sources", "pred_target", "pred_sources", "ratios", "influence_source", "influence_target"))
    nu[[f]] <- nuisance[[f]][sources]
  nu$sensitivity <- function(s, theta) nuisance$sensitivity(sources[s], theta)
  list(data = dat, nuisance = nu, components = mppi_parameter_components(dat, nu, target))
}

mppi_classic <- function(data, target, alpha, truth) {
  z <- mppi_design(target, data$x0); n <- nrow(z); H <- crossprod(z) / n
  comp <- list(S = 0L, n0 = n, H0 = H, h0 = as.numeric(crossprod(z, data$y0) / n))
  point <- mppi_solve_parameter(comp, 1)
  Q <- mppi_emp_cov(mppi_score(point$theta, z, data$y0))
  cv <- mppi_covariance(1, Q, H)
  inf <- mppi_inference(point$theta, cv$Sigma, n, alpha, truth)
  list(method = "Classic", theta = point$theta, w = c(1, rep(0, length(data$xs))),
       Sigma = cv$Sigma, interval = inf$interval, covered = inf$covered, volume = inf$volume, se = inf$se,
       diagnostics = c(point[c("score_norm", "score_tolerance", "score_passed", "min_eigen_H", "condition_H")],
                       cv[c("min_eigen_Sigma", "condition_Sigma")]))
}

mppi_fit_all_methods <- function(data, nuisance, target, truth, alpha, optimizer_seed, optimizer_control) {
  S <- length(data$xs); K <- S + 1L
  Classic <- mppi_attempt_method("Classic", mppi_classic(data, target, alpha, truth))
  MPPI <- mppi_attempt_method("MPPI", {
    ctx <- mppi_method_context(data, nuisance, target, seq_len(S))
    pilot <- mppi_solve_parameter(ctx$components, c(1, rep(0, S)))
    qpilot <- mppi_build_Q(pilot$theta, ctx$data, ctx$nuisance, ctx$components)
    wf <- mppi_optimize_weights(qpilot$A, qpilot$Q, ctx$components$n0, optimizer_seed, optimizer_control)
    fit <- mppi_fit_weight(wf$w, ctx$data, ctx$nuisance, ctx$components, alpha, truth, "MPPI")
    fit$pilot <- pilot$theta; fit$pilot_Q <- qpilot$Q; fit$optimization <- wf$diagnostics; fit
  })
  if (MPPI$status == "failed" && Classic$status == "success") {
    reason <- MPPI$message; MPPI <- Classic; MPPI$method <- "MPPI"
    MPPI$status <- "fallback"; MPPI$message <- reason
    MPPI$optimization <- list(status = "target_only_fallback", converged = FALSE)
  }
  fit_source <- function(s, plus) {
    label <- paste0(if (plus) "PPI++ source " else "PPI source ", s)
    mppi_attempt_method(label, {
      ctx <- mppi_method_context(data, nuisance, target, s)
      fit <- if (plus) mppi_ppi_plus(1L, ctx$data, ctx$nuisance, ctx$components, alpha, truth) else
        mppi_fit_weight(c(0, 1), ctx$data, ctx$nuisance, ctx$components, alpha, truth, label)
      w <- numeric(K); w[c(1L, s + 1L)] <- fit$w; fit$w <- w
      if (plus) fit$diagnostics$ppi_plus$source <- s
      fit
    })
  }
  PPI <- lapply(seq_len(S), fit_source, plus = FALSE)
  PPI_plus <- lapply(seq_len(S), fit_source, plus = TRUE)
  EW <- mppi_attempt_method("EW", {
    ctx <- mppi_method_context(data, nuisance, target, seq_len(S))
    mppi_fit_weight(rep(1 / K, K), ctx$data, ctx$nuisance, ctx$components, alpha, truth, "EW")
  })
  list(MPPI = MPPI, Classic = Classic, PPI = PPI, EW = EW, PPI_plus = PPI_plus,
       pilot = Classic$theta, nuisance_diagnostics = if (inherits(nuisance, "error"))
         list(error = conditionMessage(nuisance)) else nuisance$diagnostics)
}

mppi_method_list <- function(methods) {
  c(list(methods$MPPI, methods$Classic), methods$PPI, list(methods$EW), methods$PPI_plus)
}
