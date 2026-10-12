#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
script <- gsub("~+~", " ", script, fixed=TRUE) # Rscript encodes spaces in --file.
root_option <- grep("^--root=", args, value=TRUE)
if (length(root_option) > 1L) stop("Duplicate root option.", call.=FALSE)
root <- normalizePath(if (length(root_option)) sub("^--root=", "", root_option) else
                      file.path(dirname(script), ".."), mustWork=TRUE)
source(file.path(root, "R", "utils.R")); mppi_source_files(root)
opt <- mppi_cli("results")
results <- mppi_null(opt$results,file.path(root,"results",mppi_simulation_config()$version,"dxa"))
out <- mppi_run_dxa(root,results)
print(out$table)
if (any(out$table$Status != "success")) quit(status=2L)
