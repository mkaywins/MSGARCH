testthat::context("Test Volatility")

data("SMI", package = "MSGARCH")
spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                            distribution.spec = list(distribution = c("norm")),
                            switch.spec = list(do.mix = FALSE, K = 2))
par <- c(0.021631876185, 0.087024443479, 0.881493722371, 0.020659831566, 
         0.005396009353, 0.994040728662, 0.978348086740, 0.998703301894)

testthat::test_that("Forecast", {
  
  tol <- 0.05
  set.seed(1234)
  est.forecast <- predict(object = spec, par = par, newdata = SMI, nahead = 2)$vol
  exp.forecast <- c(1.0304257211510406, 1.0340222685323162)
  
  testthat::expect_true(max(abs(est.forecast - exp.forecast)) < tol)
  
})

testthat::test_that("Conditional Vol", {
  
  tol <- 0.05
  est.Vol <- Volatility(object = spec, par = par, data = SMI)[2000]
  exp.Vol <- c(2.1321725800180471)
  
  testthat::expect_true(max(abs(est.Vol - exp.Vol)) < tol)
  
})

## test TVP functionalities
spec <- MSGARCH::CreateSpec(switch.spec = list(do.tvp=TRUE))
Z    <- as.matrix(data.frame(x1=1, x2=lag(SMI,1)))
fit <- MSGARCH::FitML(spec, data = SMI[1:2499], Z = Z)

testthat::test_that("Forecast TVP", {
  
  tol <- 0.05
  set.seed(1234)
  est.forecast <- predict(object = fit$spec, par = fit$par, newdata = SMI[1:2499], newZ = Z, nahead = 2)$vol
  exp.forecast <- c(0.8604516, 0.8516941)
  
  testthat::expect_true(max(abs(est.forecast - exp.forecast)) < tol)
  
})

testthat::test_that("Conditional Vol MSGARCH_ML_FIT TVP", {
  
  tol <- 0.05
  est.Vol <- Volatility(fit)[2000]
  exp.Vol <- c(1.872766)
  
  testthat::expect_true(max(abs(est.Vol - exp.Vol)) < tol)
  
})

testthat::test_that("Conditional Vol MSGARCH_SPEC TVP", {
  
  tol <- 0.05
  est.Vol <- Volatility(fit$spec, par = fit$par, data = SMI[2:2500], Z = Z)[2000]
  exp.Vol <- c(1.841298)
  
  testthat::expect_true(max(abs(est.Vol - exp.Vol)) < tol)
  
})


