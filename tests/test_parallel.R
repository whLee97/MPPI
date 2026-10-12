#!/usr/bin/env Rscript
root<-normalizePath(if(file.exists('R/utils.R'))'.'else'..',mustWork=TRUE)
source(file.path(root,'R/utils.R'));mppi_source_files(root)
check<-function(ok,msg)if(!isTRUE(ok))stop(msg,call.=FALSE)
fails<-function(code)inherits(tryCatch({force(code);NULL},error=identity),'error')
if(.Platform$OS.type!='unix')stop('Parallel tests require Linux/macOS.')
cfg<-mppi_simulation_config(B=6L);cfg$code_manifest<-mppi_code_manifest(root)
serial<-tempfile('mppi-serial-');parallel_root<-tempfile('mppi-parallel-')
invisible(mppi_run_scenario(cfg,1L,1:6,serial,1L))
invisible(mppi_run_scenario(cfg,1L,1:6,parallel_root,2L))
sc<-as.list(cfg$scenarios[1,]);dir1<-mppi_scenario_dir(serial,sc);dir2<-mppi_scenario_dir(parallel_root,sc)
for(i in 1:6){a<-readRDS(file.path(dir1,sprintf('rep_%04d.rds',i)));b<-readRDS(file.path(dir2,sprintf('rep_%04d.rds',i)));check(identical(a$methods,b$methods)&&identical(a$folds,b$folds),'Serial/fork result mismatch')}
# Interrupt-like loss of an uncommitted task: existing records remain byte-for-byte unchanged on resume.
saved<-tools::md5sum(file.path(dir2,sprintf('rep_%04d.rds',1:5)))
unlink(file.path(dir2,'rep_0006.rds'))
x<-mppi_run_scenario(cfg,1L,1:6,parallel_root,2L)
check(identical(saved,tools::md5sum(names(saved)))&&sum(vapply(x,function(z)z$reused,logical(1)))==5,'Resume rewrote committed records')
# A source manager lock must exclude even a second nonoverlapping launcher.
lock<-file.path(dir2,'.scenario_lock');mppi_lock(lock,'test active manager')
check(file.exists(file.path(lock,'owner.rds'))&&fails(mppi_run_scenario(cfg,1L,1L,parallel_root,1L)),'Scenario lock did not exclude another manager');unlink(lock,recursive=TRUE)
# Corrupt output is an infrastructure error, not a statistical failure or an overwritten result.
path<-file.path(dir2,'rep_0006.rds');writeLines('corrupt',path)
x<-mppi_run_scenario(cfg,1L,6L,parallel_root,2L)
check(x[[1]]$status=='infrastructure_error'&&identical(readLines(path),'corrupt'),'Corrupt file was silently accepted or overwritten')
# Deliberately terminate only our own forked test worker. Preserve completed peers, leave missing task resumable.
original_runner<-mppi_run_replication
mppi_run_replication<-function(cfg,scenario,replicate,predictors){
 if(replicate==2L)tools::pskill(Sys.getpid(),signal=9L)
 original_runner(cfg,scenario,replicate,predictors)
}
killed<-tempfile('mppi-killed-worker-')
x<-suppressWarnings(mppi_run_scenario(cfg,1L,1:6,killed,2L))
mppi_run_replication<-original_runner
check(length(x)==4 && x[[2]]$status=='infrastructure_error','Killed worker was not detected at the batch boundary')
kdir<-mppi_scenario_dir(killed,sc)
check(file.exists(file.path(kdir,'rep_0001.rds'))&&!file.exists(file.path(kdir,'rep_0002.rds'))&&!file.exists(file.path(kdir,'rep_0005.rds')),'Killed task or undispatched tasks were misrecorded')
y<-mppi_run_scenario(cfg,1L,1:6,killed,2L)
check(all(vapply(y,function(z)z$status=='completed',logical(1))),'Killed-worker restart did not recover')
for(i in 1:6){a<-readRDS(file.path(dir1,sprintf('rep_%04d.rds',i)));b<-readRDS(file.path(kdir,sprintf('rep_%04d.rds',i)));check(identical(a$methods,b$methods),'Restart altered the seeded experiment')}
check(fails(mppi_integer_argument('2.5','B'))&&fails(mppi_integer_argument('-1','B')),'Invalid numeric arguments accepted')
# Failure during serialization does not leave a visible or stale temporary result.
bad<-tempfile('mppi-save-fail-');dir.create(bad)
check(fails(mppi_atomic_save_rds(list(a=1),file.path(bad,'x.rds'),compress='invalid')),'Injected serializer failure unexpectedly succeeded')
check(!length(list.files(bad,all.files=TRUE,no..=TRUE)),'Temporary file leaked after serialization failure')
# R-level memory errors also stop dispatch instead of becoming statistical failures/fallback.
check(fails(mppi_attempt_method('MPPI', stop('cannot allocate vector of size 1 Gb'))),'Memory error became a method failure')
mppi_run_replication<-function(...)stop('cannot allocate vector of size 1 Gb')
mem<-tempfile('mppi-memory-error-');z<-mppi_run_scenario(cfg,1L,1:6,mem,1L)
mppi_run_replication<-original_runner
check(all(vapply(z,function(x)x$status=='infrastructure_error',logical(1))) &&
      !length(list.files(mem,pattern='^rep_[0-9]+[.]rds$',recursive=TRUE)),'Memory error was committed as a statistical result')
cat('PARALLEL_AND_RECOVERY_TESTS_PASSED\n')
