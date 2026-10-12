#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
script <- gsub("~+~", " ", script, fixed=TRUE) # Rscript encodes spaces in --file.
root_option <- grep("^--root=", args, value=TRUE)
if (length(root_option) > 1L) stop("Duplicate root option.", call.=FALSE)
root <- normalizePath(if (length(root_option)) sub("^--root=", "", root_option) else
                      file.path(dirname(script), ".."), mustWork=TRUE)
source(file.path(root, "R", "utils.R")); mppi_source_files(root)
opt <- mppi_cli(c("B", "workers", "results"), list(B="1000",workers="1"))
cfg <- mppi_simulation_config(B=mppi_integer_argument(opt$B,"B"))
cfg$code_manifest <- mppi_code_manifest(root)
workers <- mppi_integer_argument(opt$workers,"workers")
results <- mppi_null(opt$results,file.path(root,"results",cfg$version))
dir.create(results,recursive=TRUE,showWarnings=FALSE)
lock <- file.path(results,".launcher_lock"); mppi_lock(lock,"all simulation scenarios")
# Top-level scripts use tryCatch(finally=...) to release the run lock.
partial <- FALSE
tryCatch({
  for (id in cfg$scenarios$id) {
    status <- mppi_run_scenario(cfg,id,seq_len(cfg$B),results,workers)
    if (any(vapply(status,function(x)x$status=="infrastructure_error",logical(1))))
      stop("Infrastructure failure; resume after resolving the reported problem.",call.=FALSE)
    partial <- partial || any(vapply(status,function(x)x$status!="completed",logical(1)))
  }
  mppi_write_summaries(results)
  mppi_export_paper(file.path(results,"summary"),file.path(results,"dxa"),file.path(results,"paper"))
},finally=unlink(lock,recursive=TRUE))
if (partial) quit(status=2L)
