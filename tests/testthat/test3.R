test_that("simple ordinal_khaos test", {
  cat("simple ordinal_khaos test\n")

  set.seed(1)

  n <- 200
  p <- 2
  X <- matrix(runif(n * p), ncol = p)

  # latent function
  f <- function(x) 2 * x[1] - x[2]
  eta <- apply(X, 1, f)

  # generate latent z and discretize into K=3 classes
  z <- eta + rnorm(n, 0, 0.3)

  # choose thresholds to ensure all classes appear
  brks <- quantile(z, probs = c(0, 1/3, 2/3, 1))
  y <- cut(z, breaks = brks, include.lowest = TRUE, labels = FALSE)

  # sanity check (should always pass with quantiles)
  expect_true(length(unique(y)) == 3)

  fit <- ordinal_khaos(
    X, y,
    nmcmc = 500,
    nburn = 250,
    thin = 2,
    verbose = FALSE
  )

  # predicted probabilities
  probs <- predict(fit, newdata = X, type = "prob", aggregate = TRUE)

  # predicted class
  yhat <- max.col(probs)

  # classification accuracy
  acc <- mean(yhat == y)

  # not too strict — just ensure model learns something
  expect_true(acc > 0.5)
})
