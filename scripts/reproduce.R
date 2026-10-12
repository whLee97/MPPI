#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
script <- gsub("~+~", " ", script, fixed=TRUE) # Rscript encodes spaces in --file.
root_option <- grep("^--root=", args, value=TRUE)
if (length(root_option) > 1L) stop("Duplicate root option.", call.=FALSE)
root <- normalizePath(if (length(root_option)) sub("^--root=", "", root_option) else
                      file.path(dirname(script), ".."), mustWork=TRUE)
source(file.path(root, "R", "utils.R")); mppi_source_files(root)
opt <- mppi_cli(c("source", "results", "output"), list(source="reference"))
mppi_assert(opt$source %in% c("reference", "run"), "source must be reference or run.")
if (opt$source == "reference") {
  summary_dir <- dxa_dir <- file.path(root,"reference")
} else {
  results <- mppi_null(opt$results, file.path(root,"results",mppi_simulation_config()$version))
  summary_dir <- file.path(results,"summary"); dxa_dir <- file.path(results,"dxa")
}
output <- mppi_null(opt$output, file.path(root,"outputs",opt$source))
mppi_export_paper(summary_dir, dxa_dir, output)
cat("Paper outputs written to", output, "\n")
