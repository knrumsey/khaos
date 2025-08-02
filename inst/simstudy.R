library(duqling)
library(BART)
library(laGP)
library(khaos)
library(orthopolynom)
library(lhs)

# BART
fit_BART <- function(X, y){
  fit <- wbart(X, y)
}
pred_BART <- function(fit, XX){
  predict(fit, XX)
}

# laGP
fit_laGP <- function(X, y){
  return(list(X=X, y=y))
}
pred_laGP <- function(fit, XX){
  n <- nrow(XX)
  preds <- matrix(NA, nrow=n, ncol=1000)
  for(i in 1:n){
    tmp <- laGP(XX[i,,drop=FALSE], start=6, end=max(10, ceiling(sqrt(n))),
                X = fit$X, Z = matrix(fit$y, ncol=1),
                g=garg(list(mle=TRUE), fit$y))
    sigma <- sqrt(tmp$s2 + tmp$mle$g)
    preds[i,] <- rnorm(1000, tmp$mean, sigma)
  }
  return(t(preds))
}

# khaos
fit_skhaos <- function(X, y){
  sparse_khaos(X, y, verbose=FALSE)
}

pred_skhaos <- function(fit, XX, verbose=FALSE){
  predict(fit, XX)
}

fit_akhaos <- function(X, y){
  adaptive_khaos(X, y)
}

pred_akhaos <- function(fit, XX){
  predict(fit, XX)
}



# Function to generate Legendre polynomials up to a specified degree
generate_legendre <- function(X, degree) {
  n <- ncol(X)  # Number of variables
  legendre_polys <- list()

  for (j in 1:n) {
    legendre_polys[[j]] <- legendre.polynomials(degree, normalized = TRUE)
  }

  # Generate polynomial features for each column of X
  poly_features <- lapply(1:n, function(j) {
    polys <- sapply(0:degree, function(d) {
      predict(legendre_polys[[j]][[d + 1]], X[, j])
    })
    polys
  })

  # Combine features across all dimensions
  poly_features_combined <- Reduce(function(a, b) cbind(a, b), poly_features)
  colnames(poly_features_combined) <- NULL  # Remove column names for simplicity
  return(poly_features_combined)
}

# Fit a PCE model with forward selection
fit_pce <- function(X, y, max_degree = 10) {
  aic <- Inf
  for(degree in 1:max_degree) {
    # Generate polynomial features
    poly_features <- generate_legendre(X, degree)
    colnames(poly_features) <- paste0("V", 1:ncol(poly_features))

    # Fit linear model using the polynomial features without an intercept
    pce_model <- lm(y ~ . - 1, data = as.data.frame(poly_features))

    # Perform stepwise selection based on AIC
    selected_model <- MASS::stepAIC(pce_model, direction = "both", trace = 0)

    curr <- AIC(selected_model)
    if(curr < aic) {
      aic <- curr
      best_model <- selected_model
      best_degree <- degree

      # Store the selected terms' indices
      selected_terms <- as.character(attr(selected_model$terms, "variables"))[-c(1, 2)]  # Remove response variable
      selected_indices <- which(colnames(poly_features) %in% selected_terms)
    }
  }

  # Return the fitted model, selected terms, and their indices
  return(list(model = best_model, degree = best_degree, sigma = sd(residuals(selected_model)), selected_indices = selected_indices))
}

# Predict using the fitted PCE model
pred_pce <- function(obj, X_new, n_samples = 1000) {
  # Generate polynomial features for the new data
  poly_features_new <- generate_legendre(X_new, obj$degree)
  colnames(poly_features_new) <- paste0("V", 1:ncol(poly_features_new))

  # Use only the selected polynomial features based on the stored indices
  poly_features_new_selected <- poly_features_new[, obj$selected_indices, drop = FALSE]

  # Use the selected model to make predictions
  y_pred <- predict(obj$model, newdata = as.data.frame(poly_features_new_selected))

  # Extract coefficients and residual variance
  coefficients <- coef(obj$model)
  sigma_sq <- obj$sigma^2
  coef_var <- vcov(obj$model)  # Variance-covariance matrix of coefficients

  # Number of new observations
  n_new <- nrow(X_new)

  # Create a matrix to store samples
  samples <- matrix(0, nrow = n_new, ncol = n_samples)

  for (i in 1:n_new) {
    # Predictive variance
    pred_var <- sigma_sq +
      t(poly_features_new_selected[i, ]) %*% coef_var %*%
      poly_features_new_selected[i, ]

    # Generate samples from the predictive distribution
    samples[i, ] <- rnorm(n_samples, mean = y_pred[i], sd = sqrt(pred_var))
  }

  return(t(samples))
}



my_plot <- function(yy, yhat){
  yhat <- t(yhat)
  row_means <- rowMeans(yhat)
  lower_ci <- apply(yhat, 1, function(x) quantile(x, 0.025))
  upper_ci <- apply(yhat, 1, function(x) quantile(x, 0.975))

  # Create the plot
  plot(yy, row_means,
       xlab = "True Values (yy)",
       ylab = "Predicted Means (rowMeans(yhat))",
       pch = 19,
       col = "dodgerblue",
       main = "Predicted Means vs True Values with 95% CI")

  # Add vertical lines for the confidence intervals
  arrows(yy, lower_ci, yy, upper_ci,
         angle = 90, code = 3, length = 0.1, col ="orange")
}


# RUN DUQLING SIM STUDY
fit_list <- pred_list <- list()
fit_list[[1]] <- fit_BART
fit_list[[2]] <- fit_laGP
fit_list[[3]] <- fit_skhaos
fit_list[[4]] <- fit_akhaos
#fit_list[[4]] <- fit_laGP
fit_list[[5]] <- fit_pce

pred_list <- list()
pred_list[[1]] <- pred_BART
pred_list[[2]] <- pred_laGP
pred_list[[3]] <- pred_skhaos
pred_list[[4]] <- pred_akhaos
#pred_list[[4]] <- pred_laGP
pred_list[[5]] <- pred_pce


results <- duqling::run_sim_study(fit_list, pred_list,
                              fnames="ishigami", n_train=c(200, 2000),
                              NSR=c(0, 0.01, 0.1),
                              replications=10,
                              mc_cores=2)


for(i in 1:10){
  X <- lhs::randomLHS(200, 3)
  Xt <- lhs::randomLHS(500, 3)
  y <- apply(X, 1, duqling::ishigami, scale01=TRUE)

  fit <- fit_list[[4]](X, y)
}


  preds <- pred_list[[i]](fit, Xt)
  print(i)
  print(dim(preds))
  print(which(is.na(preds)))
  print("")



# THERE SEEMS TO BE A PROBLEM WHEN order == p
X <- randomLHS(1000, 3)
y <- apply(X, 1, ishigami, scale01=TRUE)
fit <- khaos::adaptive_khaos(X, y, degree=10, order=3)
plot(fit)




library(conforest)
fit <- rfok(X, y)
yhat <- apply(predict(fit, X), 1, mean)

library(gplite)

fit_fitc <- function(X, y){
  n <- length(y)
  ni <- 4*floor(sqrt(n))
  ni <- max(5, ni)
  ni <- min(ni, 100, round(n/2))
  gp <- gp_init(method=method_fitc(num_inducing = ni))
  gp <- gp_optim(gp, X, y, restarts=2)
}

pred_fitc <- function(obj, Xt){
  t_pred <- gp_draw(obj, Xt, draws=1000)
  pred <- t(t_pred)
  return(pred)
}

fit_rffgp <- function(X, y){
  n <- length(y)
  nb <- 8*floor(sqrt(n))
  nb <- max(10, nb)
  nb <- min(400, nb)
  gp <- gp_init(method=method_rf(num_basis=nb))
  gp <- gp_optim(gp, X, y, restarts=2)
}

pred_rffgp <- function(obj, Xt){
  t_pred <- gp_draw(obj, Xt, draws=1000)
  pred <- t(t_pred)
  return(pred)
}

fit1 <- fit_fitc(X3, y3)
pred1 <- pred_fitc(fit1, X3)
plot(y3, colMeans(pred1))
abline(0,1,col='red')

fit2 <- fit_rffgp(X3, y3)
pred2 <- pred_rffgp(fit2, X3)
plot(y3, colMeans(pred2))
abline(0,1,col='red')

