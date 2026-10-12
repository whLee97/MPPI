mppi_generate_inference_data <- function(cfg, scenario, replicate) {
  seed <- mppi_seed(cfg$seeds$inference, scenario$id, replicate, 1L)
  mppi_with_seed(seed, {
    x0 <- mppi_rmvnorm(cfg$N0, rep(0, cfg$d), diag(cfg$d))
    if (scenario$shift == "homogeneous") {
      xs <- lapply(cfg$Ns, mppi_rmvnorm, mean = rep(0, cfg$d), sigma = diag(cfg$d))
    } else if (scenario$shift == "covariate") {
      xs <- lapply(seq_len(cfg$S), function(s) {
        mppi_rmvnorm(cfg$Ns[s], cfg$covariate_means[[s]], diag(cfg$d))
      })
    } else {
      latent <- lapply(cfg$Ns, mppi_rmvnorm, mean = rep(0, cfg$d), sigma = diag(cfg$d))
      xs <- lapply(seq_len(cfg$S), function(s) {
        sweep(latent[[s]] %*% t(cfg$domain_A[[s]]), 2L, cfg$domain_b[[s]], "+")
      })
    }
    signal <- if (scenario$shift == "domain") {
      if (scenario$target == "mean") mppi_g_domain_mean(x0, scenario$dgp) else
        mppi_g_complementary_regression(x0, scenario$dgp, cfg$regression_quadratic)
    } else if (scenario$target == "regression") {
      mppi_g_complementary_regression(x0, scenario$dgp, cfg$regression_quadratic)
    } else mppi_g_standard(x0, scenario$dgp)
    noise_sd <- if (scenario$shift == "domain" && scenario$target == "mean") {
      if (scenario$dgp == "nonlinear") cfg$mean_nonlinear$noise_sd else 7
    } else 1
    list(x0 = x0, y0 = as.numeric(signal + rnorm(cfg$N0, sd = noise_sd)), xs = xs)
  })
}

mppi_make_folds <- function(cfg, scenario, replicate, data) {
  fold0 <- mppi_near_equal_folds(cfg$N0, cfg$K,
    mppi_seed(cfg$seeds$folds, scenario$id, replicate, 1L))
  folds_s <- lapply(seq_len(cfg$S), function(s) {
    mppi_near_equal_folds(nrow(data$xs[[s]]), cfg$K,
      mppi_seed(cfg$seeds$folds, scenario$id, replicate, s + 1L))
  })
  list(target = fold0, sources = folds_s)
}

mppi_fit_simulation_nuisance <- function(cfg, scenario, data, predictors, folds) {
  parts <- lapply(seq_along(data$xs), function(s) tryCatch({
    if (scenario$shift == "homogeneous")
      return(mppi_none_nuisance(data$x0, data$xs[s], predictors[s], scenario$target))
    if (scenario$shift == "covariate")
      mppi_gaussian_dr_nuisance(data$x0, data$xs[s], predictors[s], scenario$target,
                               folds$target, folds$sources[s], cfg$K)
    else mppi_transport_nuisance(data$x0, data$xs[s], predictors[s], scenario$target,
                                 folds$target, folds$sources[s], cfg$K, cfg$finite_difference)
  }, error = mppi_catch_statistical_error))
  mppi_combine_nuisance(parts)
}

mppi_run_replication <- function(cfg, scenario, replicate, predictors) {
  started <- proc.time()[[3L]]
  data <- mppi_generate_inference_data(cfg, scenario, replicate)
  t_data <- proc.time()[[3L]]
  folds <- mppi_make_folds(cfg, scenario, replicate, data)
  nuisance <- tryCatch(mppi_fit_simulation_nuisance(cfg, scenario, data, predictors, folds), error = mppi_catch_statistical_error)
  t_nuisance <- proc.time()[[3L]]
  truth <- mppi_truth(scenario$target, scenario$shift, scenario$dgp)
  candidate_seed <- mppi_seed(cfg$seeds$candidates, scenario$id, replicate, 1L)
  methods <- mppi_fit_all_methods(data, nuisance, scenario$target, truth, cfg$alpha,
                                  candidate_seed, cfg$optimizer)
  finished <- proc.time()[[3L]]
  list(version = cfg$version, scenario = scenario, replicate = replicate,
       seeds = list(inference = mppi_seed(cfg$seeds$inference, scenario$id, replicate, 1L),
                    folds = mppi_seed(cfg$seeds$folds, scenario$id, replicate, 1L),
                    candidates = candidate_seed),
       truth = truth, methods = methods, config = cfg, folds = folds,
       rng_kind = RNGkind(),
       timing_seconds = c(data = t_data - started, nuisance = t_nuisance - t_data,
                          methods = finished - t_nuisance, total = finished - started))
}

mppi_scenario_dir <- function(results_root, scenario) {
  file.path(results_root, sprintf("scenario_%02d_%s_%s_%s", scenario$id,
                                  scenario$shift, scenario$target, scenario$dgp))
}


mppi_prepare_run <- function(cfg, results_root, code_manifest) {
  dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
  contract <- list(schema = 2L, config = cfg, code_manifest = code_manifest,
                   rng_kind = RNGkind(), R_version = as.character(getRversion()))
  signature <- mppi_object_md5(contract)
  path <- file.path(results_root, "run_manifest.rds")
  if (file.exists(path)) {
    existing <- readRDS(path)
    mppi_assert(identical(existing$run_signature, signature),
      "Result directory belongs to a different code/configuration/runtime signature. Use a new results directory.")
    return(existing)
  }
  lock <- file.path(results_root, ".manifest_lock")
  mppi_lock(lock, "initialize run manifest")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (file.exists(path)) {
    existing <- readRDS(path)
    mppi_assert(identical(existing$run_signature, signature), "Concurrent initialization signature mismatch.")
    return(existing)
  }
  previous <- list.files(results_root, pattern = "^rep_[0-9]+.*[.]rds$", recursive = TRUE)
  mppi_assert(!length(previous), "Unmanifested replication results already exist. Use a new results directory.")
  manifest <- c(contract, list(run_signature = signature, version = cfg$version,
    expected_repetitions = seq_len(cfg$B), session_info = sessionInfo(), created_at = as.character(Sys.time())))
  mppi_atomic_save_rds(manifest, path)
  manifest
}

mppi_run_scenario <- function(cfg, scenario_id, repetitions = seq_len(cfg$B),
                              results_root = file.path("results", cfg$version), workers = 1L) {
  mppi_validate_config(cfg)
  mppi_assert(length(workers) == 1L && is.finite(workers) && workers >= 1 && workers == as.integer(workers), "Invalid worker count.")
  mppi_assert(length(repetitions) > 0L && all(is.finite(repetitions)) &&
    all(repetitions == as.integer(repetitions)) && !anyDuplicated(repetitions) &&
    all(repetitions >= 1 & repetitions <= cfg$B), "Repetition IDs must be distinct integers in 1:B.")
  scenario <- cfg$scenarios[cfg$scenarios$id == scenario_id, , drop = FALSE]
  mppi_assert(nrow(scenario) == 1L, "Unknown scenario id.")
  scenario <- as.list(scenario[1L, ])
  code_manifest <- mppi_null(cfg$code_manifest, mppi_code_manifest("."))
  manifest <- mppi_prepare_run(cfg, results_root, code_manifest)
  out_dir <- mppi_scenario_dir(results_root, scenario)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # A manager owns a scenario while dispatching its workers; disjoint scenarios may run concurrently.
  lock <- file.path(out_dir, ".scenario_lock")
  mppi_lock(lock, paste("scenario", scenario_id))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  predictors <- mppi_fit_auxiliary_predictors(cfg, scenario$shift, scenario$target, scenario$dgp, scenario$id)
  scenario_signature <- mppi_object_md5(list(run = manifest$run_signature, scenario = scenario, predictors = predictors))
  sc_manifest <- list(run_signature = manifest$run_signature, scenario_signature = scenario_signature,
                      scenario = scenario, predictors = predictors)
  sc_path <- file.path(out_dir, "scenario_manifest.rds")
  if (file.exists(sc_path)) mppi_assert(identical(readRDS(sc_path), sc_manifest), "Scenario predictor/signature mismatch.")
  else mppi_atomic_save_rds(sc_manifest, sc_path)
  run_one <- function(rep) {
    path <- file.path(out_dir, sprintf("rep_%04d.rds", rep))
    if (file.exists(path)) {
      old <- readRDS(path)
      mppi_assert(identical(old$version, cfg$version) && identical(old$replicate, rep) &&
        identical(old$scenario, scenario) && identical(old$run_signature, manifest$run_signature) &&
        identical(old$scenario_signature, scenario_signature), "Existing replication signature mismatch.")
      return(list(replicate = rep, status = old$status, reused = TRUE, path = path))
    }
    ans <- tryCatch(mppi_run_replication(cfg, scenario, rep, predictors), error = mppi_catch_statistical_error)
    if (inherits(ans, "error")) {
      ans <- list(version = cfg$version, scenario = scenario, replicate = rep,
                  status = "failed", message = conditionMessage(ans), config = cfg,
                  truth = mppi_truth(scenario$target, scenario$shift, scenario$dgp))
    } else {
      statuses <- vapply(mppi_method_list(ans$methods), function(x) x$status, character(1))
      ans$status <- if (all(statuses == "success")) "completed" else "partial"
    }
    ans$run_signature <- manifest$run_signature; ans$scenario_signature <- scenario_signature
    ans$code_manifest <- code_manifest
    mppi_atomic_save_rds(ans, path, compress = "xz")
    list(replicate = rep, status = ans$status, reused = FALSE, path = path)
  }
  repetitions <- as.integer(repetitions)
  safe_run_one <- function(rep) tryCatch(run_one(rep), error = function(e)
    list(replicate = rep, status = "infrastructure_error", message = conditionMessage(e)))
  verify <- function(value, rep) {
    path <- file.path(out_dir, sprintf("rep_%04d.rds", rep))
    # A worker can disappear after committing its file. The committed record is authoritative.
    saved <- if (file.exists(path)) tryCatch(readRDS(path), error = function(e) NULL) else NULL
    if (!is.null(saved) && identical(saved$replicate, rep) && identical(saved$scenario, scenario) &&
        identical(saved$version, cfg$version) && identical(saved$run_signature, manifest$run_signature) &&
        identical(saved$scenario_signature, scenario_signature) &&
        saved$status %in% c("completed", "partial", "failed")) {
      return(list(replicate=rep, status=saved$status, path=path,
        reused=is.list(value) && isTRUE(value$reused),
        recovered_worker_output=!is.list(value) || is.null(value$status)))
    }
    reason <- if (is.list(value) && !is.null(value$message)) value$message else
      "Worker returned without a valid committed result (possible termination, memory or I/O failure)."
    list(replicate=rep, status="infrastructure_error", message=reason, path=path)
  }
  workers <- min(as.integer(workers), length(repetitions))
  if (.Platform$OS.type != "unix" && workers > 1L)
    stop("This parallel runner requires Linux/macOS; use workers=1 on other platforms.", call.=FALSE)
  batch_size <- min(length(repetitions), max(1L, 2L * workers))
  batches <- split(repetitions, ceiling(seq_along(repetitions) / batch_size))
  results <- vector("list", length(repetitions)); done <- 0L
  for (batch in batches) {
    raw <- if (workers <= 1L) lapply(batch, safe_run_one) else
      parallel::mclapply(batch, safe_run_one, mc.cores=workers, mc.preschedule=FALSE,
                        mc.set.seed=FALSE, mc.cleanup=TRUE)
    checked <- Map(verify, raw, batch)
    results[done + seq_along(batch)] <- checked; done <- done + length(batch)
    counts <- table(vapply(results[seq_len(done)], function(x) x$status, character(1)))
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "scenario", scenario_id,
        "checked", done, "/", length(repetitions),
        paste(names(counts), as.integer(counts), sep="=", collapse=" "), "\n")
    flush.console()
    if (any(vapply(checked, function(x) x$status == "infrastructure_error", logical(1)))) {
      errors <- Filter(function(x) x$status == "infrastructure_error", checked)
      for (e in errors) cat("INFRASTRUCTURE_ERROR replicate", e$replicate, e$message, "\n", file=stderr())
      # Leave unfinished tasks missing so a restart retries them with their original seed.
      return(results[seq_len(done)])
    }
  }
  results
}
