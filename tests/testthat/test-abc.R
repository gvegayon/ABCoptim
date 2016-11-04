context("Single function")

test_that("Precision", {
  set.seed(213)
  fw <- function (x)
    10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80
  
  ans0 <- abc_optim(50, fw, lb=-100, ub=100, criter=100)
  ans1 <- abc_cpp(50, fw, lb=-100, ub=100, criter=100)
  expect_equal(ans0$par, -15.81515,  tolerance = .0001, scale = 1)
  expect_equal(ans1$par, -15.81515,  tolerance = .0001, scale = 1)
})