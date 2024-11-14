library(duqling)
library(lhs)
library(khaos)

# Train on piston function
set.seed(214098)
n <- 500
p <- 5
X <- randomLHS(n, p)
y <- apply(X, 1, duqling::dms_additive, scale01=TRUE)

## SPARSE KHAOS
fit <- khaos::sparse_khaos(X, y, degree=c(2,2,14), order=c(1,1,3), max_basis=3e5)
plot(fit)

## ADAPTIVE KHAOS
fit2 <- khaos::adaptive_khaos(X, y, degree = 10, order=3)
plot(fit2)

# PREDICT ON TEST SET
Xt <- randomLHS(n, p)
yt <- apply(Xt, 1, duqling::dms_additive, scale01=TRUE)

par(mfrow=c(1,2))
preds <- predict(fit, Xt)
plot(yt, colMeans(preds), main="Sparse")

preds2 <- predict(fit2, Xt)
plot(yt, colMeans(preds2), main="Adaptive")



