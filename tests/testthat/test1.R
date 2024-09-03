test_that("simple emulation test", {
  cat('simple emulation test')
  f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
  n <- 200
  p <- 3
  X <- matrix(runif(n*p), ncol=p)
  y <- apply(X, 1, f) + rnorm(n, 0, 0.1)
  fit <- bayes_chaos(X, y)
  yhat <- colMeans(predict(fit))
  d1 <- sqrt(mean((y-yhat)^2))/sd(y)
  expect_that(d1, is_less_than(0.2))
})
