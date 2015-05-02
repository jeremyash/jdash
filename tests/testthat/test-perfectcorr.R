context("perfect correlation")

test_that("vector_corr returns a perfect correlation (= 2)", {

  perfect_corr <- 2
  set.seed(1001)
  x <- matrix(rnorm(200, 3, 1), ncol=2)
  x_x <- vector_corr(x, x)

  expect_equal(perfect_corr, unname(x_x[1]))

})