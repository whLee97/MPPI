mppi_numeric_clean <- function(x, valid_min = NULL, valid_max = NULL, missing_values = numeric()) {
  out <- suppressWarnings(as.numeric(x))
  if (length(missing_values)) out[out %in% missing_values] <- NA_real_
  if (!is.null(valid_min)) out[out < valid_min] <- NA_real_
  if (!is.null(valid_max)) out[out > valid_max] <- NA_real_
  out
}

mppi_bmi_imperial <- function(height_in, weight_lb) {
  height <- mppi_numeric_clean(height_in, 36, 90, c(777, 999, 7777, 9999))
  weight <- mppi_numeric_clean(weight_lb, 50, 700, c(777, 999, 7777, 9999))
  703 * weight / height^2
}

mppi_add_dxa_predictions <- function(dat) {
  age <- dat$x_age; male <- dat$x_male; bmi <- dat$x_bmi
  dat$f_deurenberg <- 1.20 * bmi + 0.23 * age - 10.8 * male - 5.4
  female <- 1 - male
  dat$f_cun_bae <- -44.988 + 0.503 * age + 10.689 * female + 3.172 * bmi -
    0.026 * bmi^2 + 0.181 * bmi * female - 0.020 * bmi * age -
    0.005 * bmi^2 * female + 0.00021 * bmi^2 * age
  dat$f_gallagher <- 64.5 - 848 / bmi + 0.079 * age - 16.4 * male +
    0.05 * male * age + 39 * male / bmi
  dat
}

mppi_first_existing <- function(dat, candidates) {
  lookup <- setNames(names(dat), toupper(names(dat)))
  for (candidate in candidates) if (toupper(candidate) %in% names(lookup)) return(lookup[[toupper(candidate)]])
  NULL
}

mppi_read_nhis_zip <- function(path) {
  csv_name <- utils::unzip(path, list = TRUE)$Name[1L]
  read.csv(unz(path, csv_name), check.names = FALSE)
}

mppi_read_xpt_zip <- function(path) {
  mppi_assert(requireNamespace("foreign", quietly = TRUE), "The foreign package is required to read DXA XPT files.")
  exdir <- tempfile("mppi-xpt-"); dir.create(exdir); on.exit(unlink(exdir, recursive = TRUE), add = TRUE)
  files <- utils::unzip(path, exdir = exdir)
  foreign::read.xport(files[1L])
}

mppi_build_dxa_target <- function(raw_dir) {
  mppi_assert(requireNamespace("foreign", quietly = TRUE), "The foreign package is required to read NHANES XPT files.")
  demo <- foreign::read.xport(file.path(raw_dir, "DEMO_J.xpt"))
  dxx <- foreign::read.xport(file.path(raw_dir, "DXX_J.xpt"))
  whq <- foreign::read.xport(file.path(raw_dir, "WHQ_J.xpt"))
  bmx <- foreign::read.xport(file.path(raw_dir, "BMX_J.xpt"))
  dxx <- dxx[, intersect(c("SEQN", "DXDTOPF", "DXAEXSTS"), names(dxx)), drop = FALSE]
  demo <- demo[, intersect(c("SEQN", "RIDAGEYR", "RIAGENDR"), names(demo)), drop = FALSE]
  bmx <- bmx[, intersect(c("SEQN", "BMXBMI"), names(bmx)), drop = FALSE]
  dat <- merge(merge(merge(demo, dxx, by = "SEQN", all = FALSE), whq, by = "SEQN", all.x = TRUE),
               bmx, by = "SEQN", all.x = TRUE)
  sex <- mppi_numeric_clean(dat$RIAGENDR, missing_values = c(7, 9))
  out <- data.frame(
    Y = mppi_numeric_clean(dat$DXDTOPF, 1, 70, c(7777, 9999)),
    x_age = mppi_numeric_clean(dat$RIDAGEYR, 40, 59, c(777, 999)),
    x_male = ifelse(sex == 1, 1, ifelse(sex == 2, 0, NA_real_)))
  if (all(c("WHD010", "WHD020") %in% names(dat))) {
    out$x_bmi <- mppi_bmi_imperial(dat$WHD010, dat$WHD020)
    bmi_source <- "NHANES self-reported WHD010/WHD020"
  } else {
    out$x_bmi <- mppi_numeric_clean(dat$BMXBMI, 10, 80, c(7777, 9999))
    bmi_source <- "NHANES measured BMXBMI"
  }
  out <- out[out$x_age >= 40 & out$x_age <= 59, , drop = FALSE]
  out <- out[complete.cases(out[, c("Y", "x_age", "x_male", "x_bmi")]), , drop = FALSE]
  out <- mppi_add_dxa_predictions(out)
  attr(out, "bmi_source") <- bmi_source
  out
}

mppi_standard_source_frame <- function(dat, source, year) {
  names(dat) <- toupper(names(dat))
  age_col <- mppi_first_existing(dat, c("_AGE80", "X_AGE80", "AGE_P", "AGEP_A", "AGE", "AGE_R"))
  sex_col <- mppi_first_existing(dat, c("SEX1", "SEXVAR", "_SEX", "X_SEX", "SEX", "SEX_A"))
  bmi_col <- mppi_first_existing(dat, c("_BMI5", "X_BMI5", "BMI", "BMINDX", "BMXBMI"))
  h_col <- mppi_first_existing(dat, c("AHEIGHT", "HEIGHT", "HEIGHT_I", "HTIN"))
  w_col <- mppi_first_existing(dat, c("AWEIGHTP", "WEIGHT", "WEIGHT2", "WTLB"))
  mppi_assert(!is.null(age_col) && !is.null(sex_col), paste(source, year, "is missing age or sex."))
  out <- data.frame(Y = NA_real_,
    x_age = mppi_numeric_clean(dat[[age_col]], 18, 100, c(97, 98, 99, 997, 998, 999)))
  sex <- mppi_numeric_clean(dat[[sex_col]], missing_values = c(7, 8, 9))
  out$x_male <- ifelse(sex == 1, 1, ifelse(sex == 2, 0, NA_real_))
  if (!is.null(bmi_col)) {
    bmi <- mppi_numeric_clean(dat[[bmi_col]], 10, 8000, c(7777, 9999, 77777, 99999))
    if (toupper(bmi_col) %in% c("_BMI5", "X_BMI5") || median(bmi, na.rm = TRUE) > 100) bmi <- bmi / 100
    out$x_bmi <- bmi
  } else {
    mppi_assert(!is.null(h_col) && !is.null(w_col), paste(source, year, "is missing BMI and height/weight."))
    out$x_bmi <- mppi_bmi_imperial(dat[[h_col]], dat[[w_col]])
  }
  out$source <- source; out$year <- year
  out <- out[out$x_age >= 40 & out$x_age <= 59, , drop = FALSE]
  out <- out[complete.cases(out[, c("x_age", "x_male", "x_bmi")]), , drop = FALSE]
  mppi_add_dxa_predictions(out)
}

mppi_build_dxa_sources <- function(raw_dir, processed_dir) {
  nhis <- rbind(
    mppi_standard_source_frame(mppi_read_nhis_zip(file.path(raw_dir, "nhis_2017_samadultcsv.zip")), "NHIS", 2017),
    mppi_standard_source_frame(mppi_read_nhis_zip(file.path(raw_dir, "nhis_2018_samadultcsv.zip")), "NHIS", 2018))
  brfss <- rbind(
    mppi_standard_source_frame(mppi_read_xpt_zip(file.path(raw_dir, "brfss_2017_xpt.zip")), "BRFSS", 2017),
    mppi_standard_source_frame(mppi_read_xpt_zip(file.path(raw_dir, "brfss_2018_xpt.zip")), "BRFSS", 2018))
  hints_path <- file.path(processed_dir, "source_hints_2017_2018.csv")
  mppi_assert(file.exists(hints_path), "Missing processed HINTS 2017-2018 source file.")
  hints <- read.csv(hints_path, check.names = FALSE)
  hints$Y <- NA_real_
  hints$x_age <- mppi_numeric_clean(hints$x_age, 40, 59)
  hints$x_male <- mppi_numeric_clean(hints$x_male, 0, 1)
  hints$x_bmi <- mppi_numeric_clean(hints$x_bmi, 10, 80)
  hints <- hints[complete.cases(hints[, c("x_age", "x_male", "x_bmi")]), , drop = FALSE]
  hints <- mppi_add_dxa_predictions(hints)
  list(NHIS = nhis, BRFSS = brfss, HINTS = hints)
}

mppi_load_dxa <- function(root = ".", validate = TRUE) {
  raw_dir <- file.path(root, "data", "raw_dxa")
  processed_dir <- file.path(root, "data", "processed_dxa")
  target <- mppi_build_dxa_target(raw_dir)
  sources <- mppi_build_dxa_sources(raw_dir, processed_dir)
  counts <- rbind(
    male = c(Target = sum(target$x_male == 1), vapply(sources, function(x) sum(x$x_male == 1), integer(1))),
    female = c(Target = sum(target$x_male == 0), vapply(sources, function(x) sum(x$x_male == 0), integer(1))))
  expected <- rbind(male = c(Target = 529, NHIS = 7487, BRFSS = 119843, HINTS = 886),
                    female = c(Target = 584, NHIS = 8320, BRFSS = 133665, HINTS = 1337))
  if (validate) mppi_assert(identical(dim(counts), dim(expected)) && all(counts == expected),
    paste0("DXA sample sizes differ from the accepted values. Observed: ", paste(counts, collapse = ", ")))
  list(target = target, sources = sources, counts = counts, expected_counts = expected,
       target_bmi_source = attr(target, "bmi_source"), years = "2017-2018")
}
