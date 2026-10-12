#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
script <- gsub("~+~", " ", script, fixed=TRUE) # Rscript encodes spaces in --file.
root_option <- grep("^--root=", args, value=TRUE)
if (length(root_option) > 1L) stop("Duplicate root option.", call.=FALSE)
root <- normalizePath(if (length(root_option)) sub("^--root=", "", root_option) else
                      file.path(dirname(script), ".."), mustWork=TRUE)
source(file.path(root, "R", "utils.R")); mppi_source_files(root)
opt <- mppi_cli(c("scenario", "B", "from", "to", "workers", "results"),
                list(scenario="1", B="1000", from="1", workers="1"))
B <- mppi_integer_argument(opt$B, "B")
scenario <- mppi_integer_argument(opt$scenario, "scenario")
from <- mppi_integer_argument(opt$from, "from")
to <- mppi_integer_argument(mppi_null(opt$to, as.character(B)), "to")
workers <- mppi_integer_argument(opt$workers, "workers")
mppi_assert(from <= to && to <= B, "Require 1 <= from <= to <= B.")
cfg <- mppi_simulation_config(B=B); cfg$code_manifest <- mppi_code_manifest(root)
results <- mppi_null(opt$results, file.path(root, "results", cfg$version))
status <- mppi_run_scenario(cfg, scenario, seq.int(from,to), results, workers)
print(table(vapply(status, function(x) x$status, character(1))))
if (any(vapply(status, function(x) x$status == "infrastructure_error", logical(1)))) quit(status=3L)
if (any(vapply(status, function(x) x$status != "completed", logical(1)))) quit(status=2L)
