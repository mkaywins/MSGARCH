testthat::context("Test Risk")

data("SMI", package = "MSGARCH")
spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                            distribution.spec = list(distribution = c("norm")),
                            switch.spec = list(do.mix = FALSE, K = 2))
par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566, 
         0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894)

testthat::test_that("In Sample Risk", {
  
  est.Risk <- MSGARCH::Risk(object = spec,par = par, data = SMI, do.its = TRUE)
  
  tol <- 0.05
  est.VaR <- est.Risk$VaR[2000,]
  est.ES  <- est.Risk$ES[2000,]
  exp.VaR <- c(-4.9667729046257314, -3.5094969664406950)
  exp.ES  <- c(-5.7302877512616339, -4.4191996116552694)
  
  testthat::expect_true(max(abs(cbind(est.VaR, est.ES) - c(exp.VaR, exp.ES))) < tol)
})

testthat::test_that("Out of sample Risk", {
  set.seed(1234)
  est.Risk <- MSGARCH::Risk(object = spec, par = par, data = SMI, do.its = FALSE, nahead = 2)
  
  tol <- 0.05
  est.VaR <- est.Risk$VaR[, 1]
  est.ES  <- est.Risk$ES[, 2]
  exp.VaR <- c(-2.4300333085258528, -2.4536175976412151)
  exp.ES  <- c(-2.1709270732079351, -2.2054183166214902 )
  
  testthat::expect_true(max(abs(cbind(est.VaR, est.ES) - c(exp.VaR, exp.ES))) < tol)
})

## test TVP functionalities
spec <- MSGARCH::CreateSpec(switch.spec = list(do.tvp=TRUE))
Z    <- as.matrix(data.frame(x1=1, x2=lag(SMI,1)))
fit <- MSGARCH::FitML(spec, data = SMI[1:2499], Z = Z)


testthat::test_that("In Sample Risk TVP", {
  
  est.Risk <- MSGARCH::Risk(fit$spec, par = fit$par, data = SMI[2:2500], Z = Z)
  
  tol <- 0.05
  est.VaR <- est.Risk$VaR
  est.ES  <- est.Risk$ES
  exp.VaR <- c(-1.96224, -1.386536)
  exp.ES  <- c(-2.264561, -1.75527)
  
  testthat::expect_true(max(abs(cbind(est.VaR, est.ES) - c(exp.VaR, exp.ES))) < tol)
})
