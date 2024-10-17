library(duqling)
library(lhs)

# Train on piston function
n <- 500
p <- 8
X <- randomLHS(n, p)
y <- apply(X, 1, duqling::piston, scale01=TRUE) + rnorm(n, 0, 0.1)
fit <- khaos::sparse_khaos(X, y)
plot(fit)

# Predict on test set
Xt <- randomLHS(n, p)
yt <- apply(Xt, 1, duqling::piston, scale01=TRUE) + rnorm(n, 0, 0.1)
preds <- predict(fit, Xt)
plot(yt, colMeans(preds))
abline(0,1)
