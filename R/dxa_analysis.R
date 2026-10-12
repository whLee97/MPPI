mppi_dxa_group <- function(dxa, male) {
  target <- dxa$target[dxa$target$x_male == male, , drop = FALSE]
  sources <- lapply(dxa$sources, function(x) x[x$x_male == male, , drop = FALSE])
  list(target = target, sources = sources)
}

mppi_dxa_folds <- function(group, male, K = 5L, seed = 20261005L) {
  fold0 <- mppi_near_equal_folds(nrow(group$target), K, seed + 101L + male)
  folds_s <- lapply(seq_along(group$sources), function(s) {
    mppi_near_equal_folds(nrow(group$sources[[s]]), K, seed + 1009L * s + male)
  })
  list(target = fold0, sources = folds_s)
}

mppi_dxa_method_rows <- function(result, group, threshold) {
  flatten <- function(fit) {
    if (fit$status == "failed") {
      fit$theta <- fit$se <- NA_real_; fit$Sigma <- matrix(NA_real_, 1L, 1L)
      fit$interval <- c(NA_real_, NA_real_); fit$w <- rep(NA_real_, 4L)
    }
    data.frame(Status = fit$status, Message = fit$message, Group = group, Threshold = threshold, Method = fit$method,
    theta_hat = fit$theta, Sigma = fit$Sigma[1L, 1L], SE = fit$se,
    CI_lower = fit$interval[1L], CI_upper = fit$interval[2L],
    w_Target = fit$w[1L], w_NHIS = fit$w[2L], w_BRFSS = fit$w[3L], w_HINTS = fit$w[4L],
    lambda = mppi_null(fit$lambda, NA_real_), stringsAsFactors = FALSE)
  }
  rows <- c(list(result$MPPI, result$Classic), result$PPI, list(result$EW), result$PPI_plus)
  do.call(rbind, lapply(rows, flatten))
}

mppi_dxa_pvalue_curve <- function(rows, points = 401L) {
  rows <- rows[rows$Status %in% c("success", "fallback"), , drop = FALSE]
  if (!nrow(rows)) return(data.frame(Group=character(), Threshold=numeric(), Method=character(), theta0=numeric(), p_value=numeric()))
  grid <- seq(0, 1, length.out = points)
  do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
    z <- abs(rows$theta_hat[i] - grid) / rows$SE[i]
    data.frame(Group = rows$Group[i], Threshold = rows$Threshold[i], Method = rows$Method[i],
               theta0 = grid, p_value = 2 * pnorm(z, lower.tail = FALSE))
  }))
}

mppi_run_dxa <- function(root = ".", results_root = file.path(root, "results", mppi_simulation_config()$version, "dxa"),
                         K = 5L, alpha = 0.05, seed = 20261005L,
                         optimizer_control = list(M_n = NULL, max_iter = 5000L, eps_w = 1e-7,
                                                  eta0 = 1, rho = 0.5, armijo = 0.1, max_backtrack = 60L)) {
  dir.create(results_root, recursive=TRUE, showWarnings=FALSE)
  lock <- file.path(results_root, ".dxa_lock"); mppi_lock(lock, "DXA analysis")
  on.exit(unlink(lock, recursive=TRUE), add=TRUE)
  mppi_assert(requireNamespace("foreign", quietly=TRUE), "DXA requires the recommended R package foreign.")
  inputs <- sort(c(file.path(root, "data", "raw_dxa", mppi_dxa_input_files()),
                   file.path(root, "data", "processed_dxa", "source_hints_2017_2018.csv")))
  input_checksums <- tools::md5sum(inputs)
  mppi_assert(all(!is.na(input_checksums)), "Cannot checksum DXA input files.")
  settings <- list(K=K, alpha=alpha, seed=seed, optimizer=optimizer_control, rng_kind=RNGkind(),
                   R_version=as.character(getRversion()), foreign_version=as.character(packageVersion("foreign")))
  code_manifest <- mppi_code_manifest(root)
  run_signature <- mppi_object_md5(list(code=code_manifest, inputs=unname(input_checksums), settings=settings))
  result_path <- file.path(results_root, "dxa_results.rds")
  if (file.exists(result_path)) {
    out <- readRDS(result_path)
    mppi_assert(identical(out$run_signature, run_signature), "Existing DXA output belongs to a different run; choose a new results directory.")
    # Repair missing/interrupted CSV exports from the authoritative checkpoint without fitting again.
    mppi_atomic_csv(out$table, file.path(results_root, "dxa_results.csv"))
    mppi_atomic_csv(out$pvalue_curves, file.path(results_root, "dxa_pvalue_curves.csv"))
    return(out)
  }
  dxa <- mppi_load_dxa(root, validate = TRUE)
  specs <- list(
    list(label = "Male 40-59", male = 1L, thresholds = c(28, 30, 32), main = 30),
    list(label = "Female 40-59", male = 0L, thresholds = c(38, 40, 42), main = 40))
  all_results <- list(); all_rows <- list(); idx <- 0L
  for (spec in specs) {
    group <- mppi_dxa_group(dxa, spec$male)
    folds <- mppi_dxa_folds(group, spec$male, K, seed)
    x0 <- as.matrix(group$target[, c("x_age", "x_bmi")])
    xs <- lapply(group$sources, function(x) as.matrix(x[, c("x_age", "x_bmi")]))
    dr_base <- lapply(seq_along(xs), function(s) tryCatch(
      mppi_logistic_dr_base(x0, xs[s], folds$target, folds$sources[s], K), error = mppi_catch_statistical_error))
    for (threshold in spec$thresholds) {
      predictors <- list(mppi_dxa_predictor("Deurenberg", threshold, spec$male),
                         mppi_dxa_predictor("CUN-BAE", threshold, spec$male),
                         mppi_dxa_predictor("Gallagher", threshold, spec$male))
      parts <- lapply(seq_along(xs), function(s) tryCatch({
        if (inherits(dr_base[[s]], "error")) stop(dr_base[[s]])
        mppi_logistic_dr_nuisance(x0, xs[s], predictors[s], "mean", base = dr_base[[s]])
      }, error = mppi_catch_statistical_error))
      nuisance <- mppi_combine_nuisance(parts)
      data <- list(x0 = x0, y0 = as.numeric(group$target$Y > threshold), xs = xs)
      result <- mppi_fit_all_methods(data, nuisance, "mean", truth = NULL, alpha = alpha,
        optimizer_seed = seed + 7919L * threshold + spec$male, optimizer_control = optimizer_control)
      idx <- idx + 1L
      all_results[[idx]] <- list(group = spec$label, male = spec$male, threshold = threshold,
                                 main = threshold == spec$main, result = result)
      all_rows[[idx]] <- mppi_dxa_method_rows(result, spec$label, threshold)
    }
  }
  rows <- do.call(rbind, all_rows)
  curves <- mppi_dxa_pvalue_curve(rows[rows$Threshold %in% c(30, 40), , drop = FALSE])
  out <- list(version = mppi_simulation_config()$version, code_manifest = code_manifest,
              data_diagnostics = dxa[c("counts", "expected_counts", "target_bmi_source", "years")],
              analyses = all_results, table = rows, pvalue_curves = curves)
  out$input_checksums <- input_checksums
  out$settings <- settings; out$session_info <- sessionInfo(); out$run_signature <- run_signature
  mppi_atomic_save_rds(out, result_path)
  mppi_atomic_csv(rows, file.path(results_root, "dxa_results.csv"))
  mppi_atomic_csv(curves, file.path(results_root, "dxa_pvalue_curves.csv"))
  out
}
