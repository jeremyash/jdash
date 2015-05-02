context("sample size")

test_that("vector_corr returns a warning when n < 64 ", {

  x <- matrix(rnorm(100, 3, 1), ncol=2)
  y <- matrix(rnorm(100, 2, 1), ncol=2)
  
  expect_warning(vector_corr(x, y))

})