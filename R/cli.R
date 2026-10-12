# Shared validation for command-line entry points. Options use --name=value.
mppi_cli <- function(allowed, defaults = list()) {
  args <- commandArgs(trailingOnly = TRUE)
  mppi_assert(all(grepl("^--[A-Za-z-]+=.+$", args)), "Use --name=value for every option.")
  keys <- sub("^--([^=]+)=.*$", "\\1", args)
  mppi_assert(!anyDuplicated(keys), "Duplicate command-line option.")
  mppi_assert(all(keys %in% c("root", allowed)), paste("Unknown command-line option:", paste(setdiff(keys, c("root", allowed)), collapse=", ")))
  values <- setNames(sub("^--[^=]+=", "", args), keys)
  for (key in keys) defaults[[key]] <- values[[key]]
  defaults
}

mppi_dxa_input_files <- function() c(
  "DEMO_J.xpt", "DXX_J.xpt", "WHQ_J.xpt", "BMX_J.xpt",
  "nhis_2017_samadultcsv.zip", "nhis_2018_samadultcsv.zip",
  "brfss_2017_xpt.zip", "brfss_2018_xpt.zip",
  "hints5_cycle1_2017_stata.zip", "hints5_cycle2_2018_stata.zip")
