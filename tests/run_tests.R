#!/usr/bin/env Rscript
root <- normalizePath(if (file.exists("R/utils.R")) "." else "..", mustWork = TRUE)
source(file.path(root, "R", "utils.R")); mppi_source_files(root)

check <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
near <- function(x, y, tol = 1e-7) isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE))
tests <- list()
test <- function(name, code) {
  force(code); tests[[length(tests) + 1L]] <<- name; cat("PASS", name, "\n")
}

test("empirical covariance uses divisor N", {
  x <- matrix(c(1, 2, 3, 5, 8, 13), 3, 2)
  xc <- sweep(x, 2, colMeans(x), "-")
  check(near(mppi_emp_cov(x), crossprod(xc) / 3), "Wrong covariance divisor.")
})

test("confirmed DGP values", {
  z <- matrix(c(1, 2), 1, 2)
  check(near(mppi_g_domain_mean(z, "linear"), 1 + 1.2 + 2.4), "Domain mean linear DGP changed.")
  check(near(mppi_g_standard(z, "linear"), 1.2 - 1.6), "Standard linear DGP changed.")
  check(near(mppi_g_domain_mean(z, "nonlinear"), 1 + sum(mppi_h3(z))), "Nonlinear domain mean DGP changed.")
})

test("domain mean predictor is fixed after auxiliary training", {
  set.seed(1); x <- matrix(rnorm(200), 100, 2); y <- rnorm(100)
  fit <- mppi_fit_predictor(x, y, "domain", "mean", "nonlinear", 1L)
  point <- matrix(c(0.2, -0.3), 1, 2)
  p1 <- mppi_predict(fit, point)
  p2 <- mppi_predict(fit, rbind(point, matrix(rnorm(40), 20, 2)))[1]
  check(near(p1, p2), "Evaluation-sample restandardization changed a fixed predictor.")
})

test("scalar PSD QP reaches global-gap certificate", {
  B <- matrix(c(1, 0, 1, 0, 1, 1), 2, 3); Q <- crossprod(B)
  ans <- mppi_scalar_qp(Q, 1e-9)
  check(ans$certified && all(ans$w >= -1e-12) && abs(sum(ans$w) - 1) < 1e-12,
        "Scalar QP certificate failed.")
  W <- matrix(rexp(40000 * 3), 40000, 3); W <- W / rowSums(W)
  check(ans$value <= min(rowSums((W %*% Q) * W)) + 1e-8, "Scalar QP missed a sampled feasible value.")
})

test("determinant gradient agrees with finite differences", {
  set.seed(2); p <- 2; K <- 4; Z <- matrix(rnorm(K * p * 7), K * p, 7)
  Q <- tcrossprod(Z); A <- crossprod(matrix(rnorm(12), 6, 2)) + diag(2)
  obj <- mppi_weight_objective(A, Q); w <- c(.15, .25, .35, .25)
  g <- obj$evaluate(w)$gradient; h <- 1e-6
  gn <- numeric(K)
  for (j in seq_len(K)) {e <- numeric(K); e[j] <- h; gn[j] <- (obj$evaluate(w + e, FALSE) - obj$evaluate(w - e, FALSE)) / (2 * h)}
  check(max(abs(g - gn)) < 2e-5, "Determinant gradient mismatch.")
})

test("transport inverse and sensitivity step stability", {
  set.seed(3); xs <- cbind(rnorm(80, 1, 1.2), rnorm(80, -.5, .8)); x0 <- matrix(rnorm(160), 80, 2)
  eta <- mppi_fit_transport_eta(xs, x0)
  recovered <- mppi_transport_S(mppi_transport_T(xs, eta, 2), eta, 2)
  check(sqrt(mean((recovered - xs)^2)) < 1e-9, "Affine maps are not mutual inverses.")
  pred <- mppi_fit_predictor(xs, rnorm(80), "domain", "regression", "linear", 1L)
  args <- list(x_source = xs[1:20,], x_target = x0[1:20,], predictor = pred,
               target = "regression", Ns_total = 80, N0_total = 80, rel_tol = 1e-10)
  d1 <- mppi_transport_derivative(eta, 1, args, list(relative_step=2e-5,min_step=1e-7,spd_rel_tol=1e-10))
  d2 <- mppi_transport_derivative(eta, 1, args, list(relative_step=1e-5,min_step=1e-7,spd_rel_tol=1e-10))
  check(max(abs(d1$dH-d2$dH), abs(d1$dc-d2$dc)) < 2e-5, "Transport derivative is step-unstable.")
})

test("full Q equals independent rowwise construction", {
  set.seed(4); x0 <- matrix(rnorm(200), 100, 2); xs <- list(matrix(rnorm(240),120,2), matrix(rnorm(260),130,2), matrix(rnorm(280),140,2))
  preds <- lapply(1:3, function(s) mppi_fit_predictor(matrix(rnorm(400),200,2), rnorm(200), "standard"))
  nu <- mppi_none_nuisance(x0, xs, preds, "mean"); dat <- list(x0=x0,y0=rnorm(100),xs=xs)
  comp <- mppi_parameter_components(dat,nu,"mean"); q <- mppi_build_Q(0.2,dat,nu,comp)
  G <- mppi_score(0.2,comp$z0,dat$y0)
  D <- lapply(1:3,function(s) G-mppi_score(0.2,comp$z0,nu$pred_target[[s]]))
  manual <- mppi_emp_cov(do.call(cbind,c(list(G),D)))
  for(s in 1:3){j<-s+1; manual[j,j]<-manual[j,j]+100/nrow(xs[[s]])*mppi_emp_cov(mppi_score(0.2,comp$zs[[s]],nu$pred_sources[[s]]))}
  check(near(q$Q,manual,1e-10),"Q block construction mismatch.")
})

test("fixed weights reproduce Classic and source PPI endpoints", {
  set.seed(5); x0<-matrix(rnorm(160),80,2); xs<-lapply(c(90,100,110),function(n)matrix(rnorm(2*n),n,2));
  preds<-lapply(1:3,function(s)mppi_fit_predictor(matrix(rnorm(300),150,2),rnorm(150),"standard")); dat<-list(x0=x0,y0=rnorm(80),xs=xs)
  nu<-mppi_none_nuisance(x0,xs,preds,"mean"); comp<-mppi_parameter_components(dat,nu,"mean")
  a<-mppi_fit_weight(c(1,0,0,0),dat,nu,comp,.05,0,"a"); b<-mppi_fit_weight(c(1,0,0,0),dat,nu,comp,.05,0,"b")
  c1<-mppi_fit_weight(c(0,1,0,0),dat,nu,comp,.05,0,"c"); c2<-mppi_fit_weight(c(0,1,0,0),dat,nu,comp,.05,0,"d")
  check(near(a$theta,b$theta)&&near(a$Sigma,b$Sigma)&&near(c1$theta,c2$theta)&&near(c1$Sigma,c2$Sigma),"Endpoint mismatch.")
})

test("replicate seed is execution-order invariant", {
  cfg<-mppi_simulation_config(B=1); cfg$N0<-60L; cfg$Ns<-c(70L,80L,90L); cfg$K<-2L; cfg$optimizer$M_n<-32L; cfg$optimizer$max_iter<-100L
  sc<-as.list(cfg$scenarios[1,,drop=FALSE][1,]); pred<-mppi_fit_auxiliary_predictors(cfg,sc$shift,sc$target,sc$dgp,sc$id)
  a<-mppi_run_replication(cfg,sc,7L,pred); runif(100); b<-mppi_run_replication(cfg,sc,7L,pred)
  check(near(a$methods$MPPI$theta,b$methods$MPPI$theta)&&near(a$methods$MPPI$w,b$methods$MPPI$w),"Replication changed with execution order.")
})

cat("ALL_TESTS_PASSED", length(tests), "\n")
