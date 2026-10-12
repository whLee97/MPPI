mppi_flatten_method <- function(fit, scenario, replicate, truth) {
  status <- mppi_null(fit$status, "success")
  if (!status %in% c("success", "fallback")) {
    p <- length(truth)
    fit$theta <- rep(NA_real_, p); fit$Sigma <- matrix(NA_real_, p, p)
    fit$volume <- NA_real_; fit$covered <- NA_integer_; fit$w <- rep(NA_real_, 4L)
  }
  data.frame(status = status, message = mppi_null(fit$message, ""),
    optimization_status = mppi_null(fit$optimization$status, "not_applicable"),
    score_norm = mppi_null(fit$diagnostics$score_norm, NA_real_),
    score_tolerance = mppi_null(fit$diagnostics$score_tolerance, NA_real_),
    scenario_id = scenario$id, shift = scenario$shift, target = scenario$target,
    dgp = scenario$dgp, replicate = replicate, method = fit$method,
    theta1 = fit$theta[1L], theta2 = if (length(fit$theta) > 1L) fit$theta[2L] else NA_real_,
    truth1 = truth[1L], truth2 = if (length(truth) > 1L) truth[2L] else NA_real_,
    Sigma11 = fit$Sigma[1L, 1L], Sigma12 = if (nrow(fit$Sigma) > 1L) fit$Sigma[1L, 2L] else NA_real_,
    Sigma22 = if (nrow(fit$Sigma) > 1L) fit$Sigma[2L, 2L] else NA_real_,
    volume = fit$volume, covered = fit$covered,
    w0 = fit$w[1L], w1 = fit$w[2L], w2 = fit$w[3L], w3 = fit$w[4L],
    lambda = mppi_null(fit$lambda, NA_real_), stringsAsFactors = FALSE)
}


mppi_collect_results <- function(results_root) {
  manifest_path <- file.path(results_root, "run_manifest.rds")
  mppi_assert(file.exists(manifest_path), "Missing run manifest; refusing to aggregate unversioned results.")
  manifest <- readRDS(manifest_path); cfg <- manifest$config
  mppi_assert(identical(manifest$schema, 2L), "Unsupported result schema.")
  contract <- manifest[c("schema", "config", "code_manifest", "rng_kind", "R_version")]
  mppi_assert(identical(mppi_object_md5(contract), manifest$run_signature), "Run manifest checksum mismatch.")
  mppi_assert(identical(manifest$expected_repetitions, seq_len(cfg$B)), "Manifest repetition roster disagrees with its configuration.")
  total <- nrow(cfg$scenarios) * length(manifest$expected_repetitions)
  rows <- vector("list", total); timings <- list(); expected_paths <- character(); j <- 0L
  for (i in seq_len(nrow(cfg$scenarios))) {
    scenario <- as.list(cfg$scenarios[i, ]); truth <- mppi_truth(scenario$target, scenario$shift, scenario$dgp)
    directory <- mppi_scenario_dir(results_root, scenario)
    sc_path <- file.path(directory, "scenario_manifest.rds")
    scm <- if (file.exists(sc_path)) readRDS(sc_path) else NULL
    if (!is.null(scm)) mppi_assert(identical(scm$run_signature, manifest$run_signature) &&
      identical(scm$scenario, scenario), "Scenario manifest does not match this run.")
    if (!is.null(scm)) mppi_assert(identical(scm$scenario_signature,
      mppi_object_md5(list(run=scm$run_signature, scenario=scm$scenario, predictors=scm$predictors))),
      "Scenario manifest checksum mismatch.")
    scenario_paths <- file.path(directory, sprintf("rep_%04d.rds", manifest$expected_repetitions))
    expected_paths <- c(expected_paths, scenario_paths)
    for (rep in manifest$expected_repetitions) {
      path <- scenario_paths[rep]
      x <- NULL
      if (file.exists(path)) {
        x <- readRDS(path)
        mppi_assert(!is.null(scm) && identical(x$run_signature, manifest$run_signature) &&
          identical(x$scenario_signature, scm$scenario_signature) && identical(x$version, manifest$version) &&
          identical(x$scenario, scenario) && identical(x$replicate, rep),
          paste("Result signature or task identity mismatch:", path))
      }
      if (is.null(x) || is.null(x$methods)) {
        fits <- lapply(mppi_method_labels(cfg$S), function(label) list(method = label,
          status = if (is.null(x)) "missing" else "failed",
          message = if (is.null(x)) "Not yet computed" else x$message))
      } else {
        fits <- mppi_method_list(x$methods)
        mppi_assert(identical(vapply(fits, function(z) z$method, character(1)), mppi_method_labels(cfg$S)), "Method roster mismatch.")
        if (!is.null(x$timing_seconds)) timings[[length(timings) + 1L]] <- data.frame(
          scenario_id = scenario$id, replicate = rep, t(x$timing_seconds), check.names = FALSE)
      }
      j <- j + 1L
      rows[[j]] <- do.call(rbind, lapply(fits, mppi_flatten_method, scenario=scenario, replicate=rep, truth=truth))
    }
  }
  found <- list.files(results_root, pattern = "^rep_[0-9]+.*[.]rds$", recursive = TRUE, full.names = TRUE)
  mppi_assert(!length(setdiff(found, expected_paths)), "Unexpected replication files found outside the run task roster.")
  rows <- do.call(rbind, rows)
  mppi_assert(all(rows$status %in% c("success", "fallback", "failed", "missing")), "Unknown method status.")
  list(rows = rows, timings = if (length(timings)) do.call(rbind, timings) else data.frame(),
       failures = rows[rows$status %in% c("failed", "fallback"), c("scenario_id", "replicate", "method", "status", "message")],
       successful_files = expected_paths[file.exists(expected_paths)], manifest = manifest)
}

mppi_summary_table <- function(rows, N0 = 5000L) {
  if (!nrow(rows)) return(data.frame())
  keys <- interaction(rows$scenario_id, rows$method, drop = TRUE)
  do.call(rbind, lapply(split(rows, keys), function(x) {
    valid <- x$status %in% c("success", "fallback"); z <- x[valid, , drop = FALSE]
    coverage <- if (nrow(z)) mean(z$covered) else NA_real_
    data.frame(scenario_id = x$scenario_id[1L], shift = x$shift[1L], target = x$target[1L], dgp = x$dgp[1L],
      method = x$method[1L], B_expected = nrow(x), B_success = sum(x$status == "success"),
      B_fallback = sum(x$status == "fallback"), B_failed = sum(x$status == "failed"),
      B_missing = sum(x$status == "missing"), B_valid = nrow(z),
      B_budget_stop = sum(x$optimization_status == "budget_stop"),
      coverage = coverage, coverage_mcse = if (nrow(z)) sqrt(coverage * (1 - coverage) / nrow(z)) else NA_real_,
      mean_volume = mean(z$volume), mean_theta1 = mean(z$theta1), mean_theta2 = mean(z$theta2, na.rm = TRUE),
      bias1 = mean(z$theta1 - z$truth1), empirical_sd1 = sd(z$theta1), mean_se1 = mean(sqrt(z$Sigma11 / N0)),
      bias2 = mean(z$theta2 - z$truth2), empirical_sd2 = sd(z$theta2), mean_se2 = mean(sqrt(z$Sigma22 / N0)),
      mean_w0 = mean(z$w0), mean_w1 = mean(z$w1), mean_w2 = mean(z$w2), mean_w3 = mean(z$w3),
      stringsAsFactors = FALSE)
  }))
}

mppi_mean_test_curves <- function(rows, N0 = 5000L) {
  scalar <- rows[rows$target == "mean", , drop = FALSE]
  if (!nrow(scalar)) return(data.frame(scenario_id = integer(), method = character(), theta0 = numeric(),
                                       distance = numeric(), rejection_frequency = numeric()))
  output <- list(); j <- 0L
  for (key in unique(interaction(scalar$scenario_id, scalar$method, drop = TRUE))) {
    x <- scalar[interaction(scalar$scenario_id, scalar$method, drop = TRUE) == key, , drop = FALSE]
    radius <- if (x$shift[1L] == "domain") 0.3 else 0.1
    grid <- seq(x$truth1[1L] - radius, x$truth1[1L] + radius, length.out = 41L)
    if (x$shift[1L] == "domain")
      grid <- sort(unique(round(c(grid, seq(0.9, 1.1, length.out = 41L)), 12L)))
    for (theta0 in grid) {
      reject <- N0 * (x$theta1 - theta0)^2 / x$Sigma11 > qchisq(0.95, 1L)
      j <- j + 1L
      output[[j]] <- data.frame(scenario_id = x$scenario_id[1L], method = x$method[1L],
        theta0 = theta0, distance = abs(theta0 - x$truth1[1L]), rejection_frequency = mean(reject))
    }
  }
  do.call(rbind, output)
}

mppi_regression_test_curves <- function(rows, N0 = 5000L, directions = 100L) {
  vector <- rows[rows$target == "regression", , drop = FALSE]
  if (!nrow(vector)) return(data.frame(scenario_id = integer(), method = character(), radius = numeric(),
                                       rejection_frequency = numeric()))
  angles <- seq(0, 2 * pi, length.out = directions + 1L)[-(directions + 1L)]
  unit <- cbind(cos(angles), sin(angles)); radii <- c(0, seq(0.01, 0.1, by = 0.01))
  output <- list(); j <- 0L
  for (key in unique(interaction(vector$scenario_id, vector$method, drop = TRUE))) {
    x <- vector[interaction(vector$scenario_id, vector$method, drop = TRUE) == key, , drop = FALSE]
    truth <- c(x$truth1[1L], x$truth2[1L]); rejections <- numeric(length(radii))
    for (r in seq_along(radii)) {
      nulls <- if (radii[r] == 0) matrix(truth, 1L, 2L) else sweep(unit * radii[r], 2L, truth, "+")
      d1 <- outer(x$theta1, nulls[, 1L], "-")
      d2 <- outer(x$theta2, nulls[, 2L], "-")
      determinant <- x$Sigma11 * x$Sigma22 - x$Sigma12^2
      mppi_assert(all(is.finite(determinant) & determinant > 0), "Invalid covariance in test-curve input.")
      statistic <- N0 * (x$Sigma22 * d1^2 - 2 * x$Sigma12 * d1 * d2 + x$Sigma11 * d2^2) / determinant
      reject <- statistic > qchisq(0.95, 2L)
      rejections[r] <- mean(reject)
    }
    for (r in seq_along(radii)) {
      j <- j + 1L
      output[[j]] <- data.frame(scenario_id = x$scenario_id[1L], method = x$method[1L],
                                radius = radii[r], rejection_frequency = rejections[r])
    }
  }
  do.call(rbind, output)
}

mppi_normal_diagnostics <- function(rows, N0 = 5000L) {
  z1 <- sqrt(N0) * (rows$theta1 - rows$truth1) / sqrt(rows$Sigma11)
  out1 <- data.frame(rows[, c("scenario_id", "shift", "target", "dgp", "replicate", "method")],
                     coordinate = rep(1L, nrow(rows)), z = z1)
  use2 <- is.finite(rows$theta2)
  out2 <- data.frame(rows[use2, c("scenario_id", "shift", "target", "dgp", "replicate", "method")],
                     coordinate = rep(2L, sum(use2)),
                     z = sqrt(N0) * (rows$theta2[use2] - rows$truth2[use2]) / sqrt(rows$Sigma22[use2]))
  rbind(out1, out2)
}

mppi_write_summaries <- function(results_root, output_dir = file.path(results_root, "summary")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(output_dir, ".summary_lock"); mppi_lock(lock, "write summaries")
  on.exit(unlink(lock, recursive=TRUE), add=TRUE)
  collected <- mppi_collect_results(results_root)
  n0 <- collected$manifest$config$N0
  table <- mppi_summary_table(collected$rows, n0)
  valid <- collected$rows[collected$rows$status %in% c("success", "fallback"), , drop = FALSE]
  mean_curves <- mppi_mean_test_curves(valid, n0)
  regression_curves <- mppi_regression_test_curves(valid, n0)
  normal <- mppi_normal_diagnostics(valid, n0)
  mppi_atomic_csv(collected$rows, file.path(output_dir, "replication_results.csv"))
  mppi_atomic_csv(table, file.path(output_dir, "simulation_summary.csv"))
  mppi_atomic_csv(mean_curves, file.path(output_dir, "mean_test_curves.csv"))
  mppi_atomic_csv(regression_curves, file.path(output_dir, "regression_test_curves.csv"))
  mppi_atomic_csv(normal, file.path(output_dir, "normal_diagnostics.csv"))
  mppi_atomic_csv(collected$timings, file.path(output_dir, "timings.csv"))
  failures <- collected$failures
  mppi_atomic_csv(failures, file.path(output_dir, "failures.csv"))
  invisible(list(table = table, mean_curves = mean_curves, regression_curves = regression_curves,
                 normal = normal, failures = failures))
}
