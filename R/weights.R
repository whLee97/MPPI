mppi_project_simplex <- function(v) {
  v <- as.numeric(v)
  mppi_assert(length(v) > 0L && all(is.finite(v)), "Invalid simplex projection input.")
  v <- v - max(v)
  u <- sort(v, decreasing = TRUE)
  rho <- max(which(u > (cumsum(u) - 1) / seq_along(u)))
  w <- pmax(v - (sum(u[seq_len(rho)]) - 1) / rho, 0)
  w[which.max(w)] <- w[which.max(w)] + 1 - sum(w)
  w
}

mppi_dirichlet_candidates <- function(M, K, seed) {
  mppi_with_seed(seed, {
    u <- matrix(rexp(M * K), M, K, byrow = TRUE)
    u / rowSums(u)
  })
}

mppi_weight_objective <- function(A, Q) {
  A <- mppi_sym(as.matrix(A)); Q <- mppi_sym(as.matrix(Q)); p <- nrow(A)
  mppi_assert(ncol(A) == p && nrow(Q) %% p == 0L && nrow(Q) == ncol(Q), "Incompatible A and Q dimensions.")
  K <- nrow(Q) %/% p
  qe <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  qscale <- max(abs(qe), 1)
  mppi_assert(min(qe) >= -1e-9 * qscale, "Pilot Q is not positive semidefinite within numerical tolerance.")
  Ai <- mppi_spd_solve(A, diag(p), "target Hessian")
  transform <- kronecker(diag(K), Ai)
  R <- mppi_sym(transform %*% Q %*% t(transform))
  scale <- max(abs(R), 1e-300); Rn <- R / scale
  blocks <- lapply(seq_len(K), function(j) ((j - 1L) * p + 1L):(j * p))
  covariance <- function(w, normalized = TRUE) {
    out <- matrix(0, p, p)
    for (j in seq_len(K)) for (k in seq_len(K))
      out <- out + w[j] * w[k] * Rn[blocks[[j]], blocks[[k]], drop = FALSE]
    if (normalized) mppi_sym(out) else mppi_sym(out * scale)
  }
  evaluate <- function(w, gradient = TRUE) {
    C <- covariance(w, TRUE)
    value <- if (p == 1L) C[1L, 1L] else if (p == 2L) C[1L, 1L] * C[2L, 2L] - C[1L, 2L]^2 else det(C)
    if (!gradient) return(as.numeric(value))
    if (p == 1L) {
      cof <- matrix(1, 1, 1)
    } else if (p == 2L) {
      cof <- matrix(c(C[2L, 2L], -C[1L, 2L], -C[2L, 1L], C[1L, 1L]), 2L, 2L)
    } else {
      cof <- matrix(0, p, p)
      for (i in seq_len(p)) for (j in seq_len(p)) cof[i, j] <- (-1)^(i + j) * det(C[-i, -j, drop = FALSE])
    }
    grad <- numeric(K)
    for (j in seq_len(K)) {
      dC <- matrix(0, p, p)
      for (k in seq_len(K)) {
        dC <- dC + w[k] * (Rn[blocks[[j]], blocks[[k]], drop = FALSE] +
                             Rn[blocks[[k]], blocks[[j]], drop = FALSE])
      }
      grad[j] <- sum(cof * dC)
    }
    list(value = as.numeric(value), gradient = grad, Sigma_normalized = C)
  }
  batch <- function(W) apply(as.matrix(W), 1L, evaluate, gradient = FALSE)
  list(p = p, K = K, R = Rn, scale = scale, covariance = covariance,
       evaluate = evaluate, batch = batch, min_eigen_Q = min(qe))
}

mppi_scalar_qp <- function(Q, tol = 1e-10) {
  K <- nrow(Q); best_w <- diag(K)[which.min(diag(Q)), ]; best_value <- min(diag(Q))
  faces <- skipped <- 0L
  for (m in seq_len(K)) {
    supports <- combn(seq_len(K), m, simplify = FALSE)
    for (J in supports) {
      faces <- faces + 1L
      H <- rbind(cbind(2 * Q[J, J, drop = FALSE], rep(1, m)), c(rep(1, m), 0))
      sol <- tryCatch(as.numeric(solve(H, c(rep(0, m), 1))), error = function(e) NULL)
      if (is.null(sol)) { skipped <- skipped + 1L; next }
      if (any(!is.finite(sol)) || any(sol[seq_len(m)] < -1e-10)) next
      w <- numeric(K); w[J] <- pmax(sol[seq_len(m)], 0); w <- w / sum(w)
      value <- as.numeric(crossprod(w, Q %*% w))
      if (value < best_value) { best_value <- value; best_w <- w }
    }
  }
  gradient <- as.numeric(2 * Q %*% best_w)
  gap <- max(0, sum(best_w * gradient) - min(gradient))
  list(w = best_w, value = best_value, global_gap = gap, certified = gap <= tol,
       faces = faces, singular_faces = skipped)
}

mppi_projected_gradient <- function(objective, start, control) {
  w <- start; value <- objective$evaluate(w, FALSE); initial <- value
  status <- "budget_stop"; steps <- 0L; evaluations <- 1L; eta_last <- NA_real_
  residual <- Inf
  for (iteration in seq_len(control$max_iter)) {
    vg <- objective$evaluate(w); evaluations <- evaluations + 1L
    eta <- control$eta0; accepted <- FALSE
    for (bt in seq_len(control$max_backtrack)) {
      candidate <- mppi_project_simplex(w - eta * vg$gradient)
      step <- candidate - w
      candidate_value <- objective$evaluate(candidate, FALSE); evaluations <- evaluations + 1L
      if (is.finite(candidate_value) && candidate_value <= vg$value - control$armijo * sum(step^2) / eta) {
        accepted <- TRUE; break
      }
      eta <- eta * control$rho
    }
    if (!accepted) { status <- "line_search_or_roundoff_stop"; break }
    residual <- sqrt(sum(step^2)) / eta
    w <- candidate; value <- candidate_value; steps <- steps + 1L; eta_last <- eta
    if (residual <= control$eps_w) { status <- "stationary_tolerance"; break }
  }
  vg <- objective$evaluate(w)
  unit_residual <- sqrt(sum((mppi_project_simplex(w - vg$gradient) - w)^2))
  list(w = w, value = vg$value, initial_value = initial, residual = residual,
       unit_residual = unit_residual, residual_definition = "accepted_step_norm / accepted_step_size",
       status = status, converged = status == "stationary_tolerance", iterations = steps,
       evaluations = evaluations, eta_last = eta_last)
}

mppi_optimize_weights <- function(A, Q, n, seed, control = list()) {
  defaults <- list(M_n = max(256L, ceiling(log(n + 1)^2)), max_iter = 5000L, eps_w = 1e-7,
                   qp_tol = min(1e-10, 1 / n), eta0 = 1, rho = 0.5, armijo = 0.1,
                   max_backtrack = 60L)
  for (nm in names(control)) defaults[[nm]] <- control[[nm]]
  ctl <- defaults
  if (is.null(ctl$M_n)) ctl$M_n <- max(256L, ceiling(log(n + 1)^2))
  obj <- mppi_weight_objective(A, Q)
  if (obj$p == 1L) {
    qp <- mppi_scalar_qp(obj$R, ctl$qp_tol)
    mppi_assert(qp$certified, "Scalar QP did not achieve its convex global-gap tolerance.")
    return(list(w = qp$w, objective_scaled = qp$value,
      diagnostics = list(method = "scalar_QP_faces", status = "global_gap_tolerance",
        converged = TRUE, global_gap_scaled = qp$global_gap,
        global_gap = qp$global_gap * obj$scale, faces = qp$faces,
        singular_faces = qp$singular_faces, iterations = 0L, M_n = 0L,
        objective_scaled = qp$value, log_objective_scale = log(obj$scale),
        min_eigen_Q = obj$min_eigen_Q)))
  }
  candidates <- mppi_dirichlet_candidates(as.integer(ctl$M_n), obj$K, seed)
  candidates <- rbind(candidates, rep(1 / obj$K, obj$K), diag(obj$K))
  values <- obj$batch(candidates)
  mppi_assert(all(is.finite(values)), "Nonfinite determinant among candidate weights.")
  start <- candidates[which.min(values), ]
  pg <- mppi_projected_gradient(obj, start, ctl)
  mppi_assert(pg$value <= min(values) + 1e-12 * max(1, abs(min(values))),
              "Projected gradient increased the best candidate objective.")
  list(w = pg$w, objective_scaled = pg$value,
       diagnostics = list(method = "dirichlet_candidates_projected_gradient", status = pg$status,
         converged = pg$converged, residual = pg$residual, unit_residual = pg$unit_residual,
         residual_definition = pg$residual_definition, iterations = pg$iterations,
         evaluations = pg$evaluations, eta_last = pg$eta_last, M_n = ctl$M_n,
         initial_value_scaled = min(values), objective_scaled = pg$value,
         log_objective_scale = obj$p * log(obj$scale), min_eigen_Q = obj$min_eigen_Q))
}
