testthat::context("Test Extract State Fit")


## test non-TVP functionalities
data("SMI", package = "MSGARCH")
spec <- MSGARCH::CreateSpec(variance.spec = list(model = c("sGARCH")),
                            distribution.spec = list(distribution = c("norm")),
                            switch.spec = list(do.mix = FALSE, K = 2))
fit <- MSGARCH::FitML(spec, data = SMI)

testthat::test_that("ExtractStateFit", {
 
  est.len = length(MSGARCH::ExtractStateFit(fit))
  exp.len = 2L
  
  testthat::expect_true(est.len == exp.len)
})

## test TVP functionalities
data("SMI", package = "MSGARCH")
spec <- MSGARCH::CreateSpec(switch.spec = list(do.tvp=TRUE))
Z    <- as.matrix(data.frame(x1=1, x2=lag(SMI,1)))
fit <- MSGARCH::FitML(spec, data = SMI[1:2499], Z = Z)

testthat::test_that("ExtractStateFit TVP", {
  
  est.len = length(MSGARCH::ExtractStateFit(fit))
  exp.len = 2L
  
  testthat::expect_true(est.len == exp.len)
})
