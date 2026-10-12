#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
script <- gsub("~+~", " ", script, fixed=TRUE) # Rscript encodes spaces in --file.
root_option <- grep("^--root=", args, value=TRUE)
if (length(root_option) > 1L) stop("Duplicate root option.", call.=FALSE)
root <- normalizePath(if (length(root_option)) sub("^--root=", "", root_option) else
                      file.path(dirname(script), ".."), mustWork=TRUE)
source(file.path(root, "R", "utils.R")); mppi_source_files(root)
opt <- mppi_cli(c("B", "results", "dxa"), list(B="1000",dxa="false"))
mppi_assert(opt$dxa %in% c("true","false"), "dxa must be true or false.")
cfg <- mppi_simulation_config(B=mppi_integer_argument(opt$B,"B"))
cfg$code_manifest <- mppi_code_manifest(root)
results <- mppi_null(opt$results,file.path(root,"results",cfg$version))
print(sessionInfo())
if (opt$dxa == "true") {
  mppi_assert(requireNamespace("foreign",quietly=TRUE), "DXA requires R package foreign.")
  print(mppi_load_dxa(root,validate=TRUE)$counts)
}
dir.create(results,recursive=TRUE,showWarnings=FALSE)
probe <- tempfile(".write-probe-",tmpdir=results)
mppi_atomic_save_rds(list(ok=TRUE),probe)
mppi_assert(isTRUE(readRDS(probe)$ok), "Output filesystem read/write verification failed.")
unlink(probe)
invisible(mppi_prepare_run(cfg,results,cfg$code_manifest))
cat("PREFLIGHT_PASSED\n")
