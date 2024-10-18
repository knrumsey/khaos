library(duqling)
library(lhs)

# Train on piston function
n <- 500
p <- 8
X <- randomLHS(n, p)
y <- apply(X, 1, duqling::piston, scale01=TRUE) + rnorm(n, 0, 0.1)

## SPARSE KHAOS
fit <- khaos::sparse_khaos(X, y, degree=c(2,2,14), order=c(1,1,3), max_basis=3e5)
plot(fit)

# Predict on test set
Xt <- randomLHS(n, p)
yt <- apply(Xt, 1, duqling::piston, scale01=TRUE) + rnorm(n, 0, 0.1)
preds <- predict(fit, Xt)
plot(yt, colMeans(preds))
abline(0,1)

## ADAPTIVE KHAOS

