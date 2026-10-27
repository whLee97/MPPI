mppi_null <- function(x, y) if (is.null(x)) y else x

mppi_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

mppi_as_matrix <- function(x, name = deparse(substitute(x))) {
  x <- as.matrix(x)
  mppi_assert(is.numeric(x) && length(x) > 0L && all(is.finite(x)),
              paste0(name, " must be a finite numeric matrix."))
  x
}

mppi_sym <- function(x) (x + t(x)) / 2

mppi_emp_cov <- function(x) {
  x <- mppi_as_matrix(x)
  n <- nrow(x)
  mppi_assert(n >= 2L, "At least two rows are required for covariance estimation.")
  xc <- sweep(x, 2L, colMeans(x), "-")
  mppi_sym(crossprod(xc) / n)
}

mppi_eigen_diagnostics <- function(x) {
  values <- eigen(mppi_sym(x), symmetric = TRUE, only.values = TRUE)$values
  scale <- max(abs(values), 1)
  list(min = min(values), max = max(values), condition = if (min(values) > 0) max(values) / min(values) else Inf,
       tolerance = 100 * .Machine$double.eps * scale)
}

mppi_spd_solve <- function(a, b, name = "matrix", rel_tol = 1e-11) {
  a <- mppi_sym(mppi_as_matrix(a, name))
  ev <- eigen(a, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(abs(ev), 1)
  mppi_assert(min(ev) > rel_tol * scale,
              sprintf("%s is not numerically positive definite (minimum eigenvalue %.4e).", name, min(ev)))
  r <- chol(a)
  backsolve(r, forwardsolve(t(r), b))
}

mppi_matrix_power <- function(a, power, name = "covariance", rel_tol = 1e-11) {
  a <- mppi_sym(mppi_as_matrix(a, name))
  ee <- eigen(a, symmetric = TRUE)
  scale <- max(abs(ee$values), 1)
  mppi_assert(min(ee$values) > rel_tol * scale,
              sprintf("%s is not positive definite (minimum eigenvalue %.4e).", name, min(ee$values)))
  mppi_sym(ee$vectors %*% diag(ee$values^power, nrow(a)) %*% t(ee$vectors))
}

mppi_cov_n <- function(x) {
  x <- mppi_as_matrix(x)
  mu <- colMeans(x)
  xc <- sweep(x, 2L, mu, "-")
  mppi_sym(crossprod(xc) / nrow(x))
}

mppi_near_equal_folds <- function(n, k, seed) {
  mppi_assert(n >= k && k >= 2L, "Cross-fitting requires n >= K >= 2.")
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))
  sample(rep(seq_len(k), length.out = n))
}

mppi_with_seed <- function(seed, code) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old, envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

mppi_rmvnorm <- function(n, mean, sigma) {
  mean <- as.numeric(mean)
  sigma <- mppi_sym(as.matrix(sigma))
  z <- matrix(rnorm(n * length(mean)), n, length(mean))
  sweep(z %*% chol(sigma), 2L, mean, "+")
}

mppi_seed <- function(base, scenario, replicate, stage) {
  z <- (as.double(base) + 1000003 * as.double(scenario) + 1009 * as.double(replicate) +
          97 * as.double(stage)) %% 2147483646
  as.integer(z + 1)
}

mppi_score <- function(theta, z, y) {
  z <- as.matrix(z)
  sweep(z, 1L, as.numeric(z %*% theta) - as.numeric(y), "*")
}

mppi_vech <- function(a) {
  a <- as.matrix(a)
  a[lower.tri(a, diag = TRUE)]
}

mppi_unvech <- function(v, d) {
  a <- matrix(0, d, d)
  a[lower.tri(a, diag = TRUE)] <- v
  a <- a + t(a) - diag(diag(a))
  a
}

mppi_atomic_save_rds <- function(object, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(paste0(".", basename(path), ".tmp-"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  mppi_assert(file.rename(tmp, path), paste0("Could not atomically move result to ", path))
  invisible(path)
}

mppi_object_md5 <- function(object) {
  path <- tempfile("mppi-signature-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 3)
  unname(tools::md5sum(path))
}

mppi_code_manifest <- function(root = ".") {
  files <- c(list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE),
             list.files(file.path(root, "scripts"), pattern = "[.](R|py|sh)$", full.names = TRUE),
             list.files(file.path(root, "config"), pattern = "[.]R$", full.names = TRUE),
             list.files(file.path(root, "tests"), pattern = "[.](R|py)$", full.names = TRUE),
             file.path(root, "README.md"))
  files <- sort(files[file.exists(files)])
  values <- unname(tools::md5sum(files))
  names(values) <- substring(normalizePath(files), nchar(normalizePath(root)) + 2L)
  values
}

mppi_source_files <- function(root = ".") {
  files <- c("utils.R", "config.R", "predictors.R", "nuisance_covariate.R",
             "nuisance_transport.R", "weights.R", "estimation.R", "methods.R",
             "simulation.R", "dxa_data.R", "dxa_analysis.R", "summaries.R", "cli.R", "plot_domain_raw_zoom.R", "paper_outputs.R")
  for (f in files) source(file.path(root, "R", f), local = .GlobalEnv)
  invisible(files)
}


mppi_lock <- function(path, purpose) {
  mppi_assert(dir.create(path, showWarnings = FALSE),
    paste("Resource is locked:", path, "Inspect owner.rds and verify the owner has stopped before clearing a stale lock."))
  owner <- list(host = unname(Sys.info()["nodename"]), pid = Sys.getpid(),
                started = as.character(Sys.time()), purpose = purpose)
  tryCatch(mppi_atomic_save_rds(owner, file.path(path, "owner.rds")), error = function(e) {
    unlink(path, recursive = TRUE); stop(e)
  })
  invisible(path)
}

mppi_integer_argument <- function(value, name, minimum = 1L) {
  mppi_assert(length(value) == 1L && !is.na(value) && grepl("^[0-9]+$", as.character(value)),
              paste(name, "must be an integer."))
  n <- suppressWarnings(as.numeric(value))
  mppi_assert(is.finite(n) && n >= minimum && n <= .Machine$integer.max,
              paste(name, "is outside its allowed range."))
  as.integer(n)
}

mppi_atomic_csv <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(paste0(".", basename(path), ".tmp-"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  write.csv(object, tmp, row.names = FALSE)
  mppi_assert(file.rename(tmp, path), paste("Cannot commit CSV:", path))
}


mppi_catch_statistical_error <- function(e) {
  # Resource exhaustion must stop dispatch; it is not a statistical method failure.
  if (grepl("cannot allocate|memory exhausted|out of memory|std::bad_alloc|no space left on device|too many open files",
            conditionMessage(e), ignore.case=TRUE)) stop(e)
  e
}
