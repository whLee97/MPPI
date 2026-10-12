#!/usr/bin/env Rscript
root <- normalizePath(if (file.exists('R/utils.R')) '.' else '..', mustWork=TRUE)
source(file.path(root,'R/utils.R'));mppi_source_files(root)
check <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)
fails <- function(code) inherits(tryCatch({force(code);NULL},error=identity),'error')
set.seed(64); dat<-list(x0=matrix(rnorm(160),80,2),y0=rnorm(80),xs=lapply(c(90,100,110),function(n)matrix(rnorm(2*n),n,2)))
preds<-lapply(1:3,function(s)mppi_fit_predictor(matrix(rnorm(200),100,2),rnorm(100),'standard'))
nu<-mppi_none_nuisance(dat$x0,dat$xs,preds,'mean')
fitall<-function(nuisance=nu)mppi_fit_all_methods(dat,nuisance,'mean',0,.05,132,list())
normal<-fitall()
check(all(vapply(mppi_method_list(normal),function(x)x$status=='success',logical(1))),'Normal methods did not succeed')
# Method-level error isolation.
original_fit<-mppi_fit_weight
mppi_fit_weight<-function(w,data,nuisance,components,alpha,truth,label){if(label=='PPI source 2')stop('injected method error');original_fit(w,data,nuisance,components,alpha,truth,label)}
bad<-fitall();mppi_fit_weight<-original_fit
check(bad$PPI[[2]]$status=='failed' && all(vapply(mppi_method_list(bad)[-4],function(x)x$status=='success',logical(1))),'One method failure affected other methods')
# A single source nuisance fails without poisoning other sources or Classic.
badnu<-nu;badnu$source_errors<-list(NULL,'injected source error',NULL)
x<-fitall(badnu)
check(x$MPPI$status=='fallback' && x$Classic$status=='success' && x$EW$status=='failed','Incorrect common-source failure states')
check(x$PPI[[1]]$status=='success' && x$PPI[[2]]$status=='failed' && x$PPI[[3]]$status=='success','Source PPI isolation failed')
check(x$PPI_plus[[1]]$status=='success' && x$PPI_plus[[2]]$status=='failed' && x$PPI_plus[[3]]$status=='success','Source PPI++ isolation failed')
check(identical(x$MPPI$theta,x$Classic$theta),'MPPI fallback is not target-only')
x<-fitall(simpleError('shared nuisance error'))
check(x$Classic$status=='success' && x$MPPI$status=='fallback','Classic depends on shared nuisance success')
# Real per-source nuisance assembly also preserves healthy sources.
cf<-mppi_simulation_config(B=1L);sc<-as.list(cf$scenarios[1,])
broken_preds<-preds;broken_preds[[2]]<-list(kind='invalid')
assembled<-mppi_fit_simulation_nuisance(cf,sc,dat,broken_preds,list())
x<-fitall(assembled)
check(x$PPI[[1]]$status=='success' && x$PPI[[2]]$status=='failed' && x$PPI[[3]]$status=='success','Nuisance assembly is not isolated')
# Deliberately inaccurate linear solve must be rejected by the score check.
comp<-mppi_parameter_components(dat,nu,'mean'); original_solve<-mppi_spd_solve
mppi_spd_solve<-function(a,b,...) matrix(100,nrow(a),ncol(as.matrix(b)))
rejected<-fails(mppi_solve_parameter(comp,c(1,0,0,0)));mppi_spd_solve<-original_solve
check(rejected,'Inaccurate parameter solution was accepted')
check(normal$MPPI$diagnostics$score_passed && normal$MPPI$diagnostics$score_norm<=1/80,'Score acceptance diagnostics missing')
# Candidate budgets are nested, and the accepted-step residual uses actual eta.
check(identical(mppi_dirichlet_candidates(10,4,9),mppi_dirichlet_candidates(20,4,9)[1:10,]),'Candidate prefixes are not nested')
obj<-list(evaluate=function(w,gradient=TRUE){v<-sum(w^2);if(gradient)list(value=v,gradient=2*w)else v})
ctl<-list(max_iter=1L,eps_w=1e-12,eta0=.1,rho=.5,armijo=.1,max_backtrack=20L)
pg<-mppi_projected_gradient(obj,c(.9,.1),ctl)
check(abs(pg$residual-sqrt(sum((pg$w-c(.9,.1))^2))/pg$eta_last)<1e-12,'Wrong accepted-step residual')
check(pg$status=='budget_stop','Budget stop was mislabeled convergence')
# Manifest, missing-task counts, deterministic resume, version/signature rejection.
cfg<-mppi_simulation_config(B=2L);cfg$code_manifest<-mppi_code_manifest(root)
tmp<-tempfile('mppi-framework-test-');dir.create(tmp)
a<-mppi_run_scenario(cfg,1L,1L,tmp,1L);b<-mppi_run_scenario(cfg,1L,1L,tmp,1L)
check(b[[1]]$reused && b[[1]]$status=='completed','Resume did not reuse completed result')
rows<-mppi_collect_results(tmp)$rows
check(nrow(rows)==12*2*9 && sum(rows$status=='success')==9 && sum(rows$status=='missing')==207,'Expected task roster missing')
changed<-cfg;changed$alpha<-.1
check(fails(mppi_run_scenario(changed,1L,2L,tmp,1L)),'Changed run signature was accepted')
path<-a[[1]]$path;saved<-readRDS(path);tampered<-saved;tampered$run_signature<-'old-signature';saveRDS(tampered,path)
check(fails(mppi_collect_results(tmp)),'Mixed signature was aggregated');saveRDS(saved,path)
# Persisted partial records retain all method rows and expose their denominators.
saved$methods$PPI[[2]]<-list(method='PPI source 2',status='failed',message='injected persisted failure');saved$status<-'partial';saveRDS(saved,path)
summary<-mppi_write_summaries(tmp)$table
z<-summary[summary$scenario_id==1 & summary$method=='PPI source 2',]
check(z$B_expected==2 && z$B_failed==1 && z$B_missing==1 && z$B_valid==0,'Failure denominators are wrong')
# Whole-run failures and completely missing runs remain summarizable.
unlink(path);original_runner<-mppi_run_replication;mppi_run_replication<-function(...)stop('injected whole-rep failure')
invisible(mppi_run_scenario(cfg,1L,1L,tmp,1L));mppi_run_replication<-original_runner
z<-mppi_write_summaries(tmp)$table
check(sum(z$B_failed)==9 && sum(z$B_valid)==0,'All-failed run summary failed')
cat('FRAMEWORK_TESTS_PASSED\n')
