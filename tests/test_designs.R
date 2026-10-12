#!/usr/bin/env Rscript
root <- normalizePath('.', mustWork=TRUE)
source(file.path(root,'R/utils.R'));mppi_source_files(root)
check <- function(ok,msg) if(!isTRUE(ok)) stop(msg,call.=FALSE)
# Five-point Gaussian quadrature is exact through polynomial degree nine.
nodes <- c(-sqrt(5+sqrt(10)),-sqrt(5-sqrt(10)),0,sqrt(5-sqrt(10)),sqrt(5+sqrt(10)))
weights <- c((7-2*sqrt(10))/60,(7+2*sqrt(10))/60,8/15,(7+2*sqrt(10))/60,(7-2*sqrt(10))/60)
index <- expand.grid(i=1:5,j=1:5)
x <- cbind(nodes[index$i],nodes[index$j]); w <- weights[index$i]*weights[index$j]
h <- mppi_h3(nodes)
check(max(abs(c(sum(weights*h),sum(weights*h*nodes),sum(weights*h*(nodes^2-1)),sum(weights*h^2)-1)))<1e-12,'H3 mean/variance/orthogonality failed.')
g <- mppi_g_domain_mean(x,'nonlinear')
check(abs(sum(w*g)-1)<1e-12 && abs(sum(w*(g-1)^2)-2)<1e-12,'Incorrect target mean or signal variance.')
cfg <- mppi_simulation_config(B=2)
check(identical(cfg$scenarios$id,seq_len(12L)),'Wrong scenario roster.')
check(identical(cfg$domain_A,mppi_simulation_config()$domain_A) && identical(cfg$domain_b,mppi_simulation_config()$domain_b),'Affine shift was modified.')

# Verify the sensitivity of the actual affine map, including translation and full covariance.
for(s in 1:2) {
  mu_s <- cfg$domain_b[[s]]; sigma_s <- tcrossprod(cfg$domain_A[[s]])
  eta <- c(mu_s,mppi_vech(sigma_s),rep(0,2),mppi_vech(diag(2)))
  pred <- structure(list(kind='domain',target='mean',dgp='nonlinear',source=s,
      preprocessing=list(center=mu_s,scale=diag(cfg$domain_A[[s]])),
      coefficients=c(1,0,1)),class='mppi_predictor')
  average_prediction <- function(e) sum(w*mppi_predict(pred,mppi_transport_S(x,e,2)))
  deriv <- vapply(seq_along(eta),function(j){ep<-em<-eta;ep[j]<-ep[j]+1e-5;em[j]<-em[j]-1e-5;(average_prediction(ep)-average_prediction(em))/2e-5},numeric(1))
  check(max(abs(deriv))<1e-7,'Oracle cubic predictor has nonzero affine calibration sensitivity.')
}
# Compute the oracle minimum with the implemented simplex solver.
Q <- matrix(c(3,2,2,3,2,2.5,1,2,2,1,2.25,2,3,2,2,3),4,4)
oracle <- mppi_scalar_qp(Q,1e-10)
check(oracle$certified && max(abs(oracle$w-c(0,5/11,6/11,0)))<1e-8 && abs(oracle$value-(1+15/22))<1e-10,'Incorrect oracle simplex benchmark.')


w <- weights[index$i]*weights[index$j]
g <- mppi_g_complementary_regression(x, 'nonlinear')
check(max(abs(colSums(x*(w*g))-c(1.2,-0.8))) < 1e-12, 'New DGP changed the target projection.')
r <- g - as.numeric(x %*% c(1.2,-0.8)); score <- x*r
Q <- crossprod(score, score*w) + diag(2)
check(max(abs(Q-13*diag(2))) < 1e-11, 'Classic population score covariance is not 13 I.')
for (s in 1:3) {
 z <- mppi_feature_complementary(x,s)
 check(ncol(z) == if (s==3) 3L else 4L, 'Wrong source feature dimension.')
 if(s<3) check(max(abs(z[,4]-(x[,s]^2-1)))<1e-12, 'Sources are not using complementary coordinates.')
}

# Full sample sizes and all nine methods in every formal scenario.
for (id in cfg$scenarios$id) {
  sc <- as.list(cfg$scenarios[cfg$scenarios$id==id,])
  pred <- mppi_fit_auxiliary_predictors(cfg,sc$shift,sc$target,sc$dgp,sc$id)
  fit <- mppi_run_replication(cfg,sc,1L,pred)
  fits <- mppi_method_list(fit$methods)
  check(all(vapply(fits,function(z)z$status=="success",logical(1))),paste("Method failure in scenario",id))
  check(all(vapply(fits,function(z)isTRUE(z$diagnostics$score_passed),logical(1))),"Parameter score acceptance failed.")
  cat("PASS formal design",id,"\n")
}
cat("DESIGN_TESTS_PASSED\n")
