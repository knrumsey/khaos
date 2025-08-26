library(BART)
library(cbass)
library(khaos)
library(ordinal)

qwk <- function(truth, pred, K = NULL) {
  # Ensure inputs are integers
  truth <- as.integer(truth)
  pred <- as.integer(pred)

  if (is.null(K)) {
    K <- max(c(truth, pred))
  }

  # Create observed matrix
  O <- table(factor(truth, levels = 1:K),
             factor(pred, levels = 1:K))
  O <- as.matrix(O)

  # Normalize to frequencies
  O <- O / sum(O)

  # Expected matrix (outer product of marginals)
  row_marg <- rowSums(O)
  col_marg <- colSums(O)
  E <- outer(row_marg, col_marg)

  # Weight matrix
  W <- outer(1:K, 1:K, FUN = function(i, j) (i - j)^2 / (K - 1)^2)

  # Compute QWK
  kappa <- 1 - sum(W * O) / sum(W * E)
  return(kappa)
}

mean_qwk <- function(sample_matrix, truth) {
  n_draws <- nrow(sample_matrix)
  qwk_vals <- numeric(n_draws)
  for (i in 1:n_draws) {
    qwk_vals[i] <- qwk(truth, sample_matrix[i, ])
  }
  mean(qwk_vals)
}

crossval_compare <- function(X, y, k = 5, ndraws = 1000) {
  folds <- caret::createFolds(y, k = k, list = TRUE)
  results <- list()
  eqwk_all <- list()

  for (fold in seq_along(folds)) {
    test_idx <- folds[[fold]]
    train_idx <- setdiff(seq_along(y), test_idx)

    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_test  <- X[test_idx, , drop = FALSE]
    y_test  <- y[test_idx]

    y_fact <- factor(y_train, ordered = TRUE)

    # Fit models
    fit_cbass <- fit.cbass(X_train, y_train)
    fit_khaos <- ordinal_khaos(X_train, y_train, h1=40, h2=0.01/length(y), verbose=FALSE, nmcmc=20000, nburn=18000, thin=2)
    df_train <- data.frame(y = y_fact, X_train)
    fit_clm <- ordinal::clm(y ~ ., data = df_train)

    #fit_bart <- mbart2(X_train, y_train, X_test)
    fit_bart <- mbart(X_train, y_train, X_test)

    # Predict
    # 1. cbass (fix: convert latent z to class probs, then sample)
    predc_chain <- pred.cbass(fit_cbass, X_test)  # [nmcmc × n × latent_z]
    n_mcmc <- dim(predc_chain)[1]
    n_test <- dim(predc_chain)[2]
    n_class <- length(unique(y_train))  # or set explicitly, e.g. 6 or 7

    # Allocate predictive class draw array
    predc_samps <- array(0, dim = c(n_mcmc, n_test))
    for(i in 1:n_test){
      for(j in 1:n_mcmc){
        prob_j <- pmax(0, p.mu(predc_chain[j, i, ]))
        predc_samps[j,i] <- sample(1:n_class, 1, prob=prob_j)
      }
    }

    # 2. khaos
    predp_array <- predict(fit_khaos, X_test, type = "prob")
    predp_mean <- apply(predp_array, c(2, 3), mean)
    predp_class <- apply(predp_mean, 1, which.max)

    n_mcmc <- 1000
    predp_samps <- array(0, dim = c(n_mcmc, n_test))
    for(i in 1:n_test){
      for(j in 1:n_mcmc){
        prob_j <- predp_array[j, i, ]
        predp_samps[j,i] <- sample(1:n_class, 1, prob=prob_j)
      }
    }

    # 3. clm
    df_test <- data.frame(X_test)
    pred_clm <- predict(fit_clm, newdata = df_test, type = "prob")
    pred_clm_class <- apply(pred_clm$fit, 1, which.max)

    # Get the class probabilities from clm (shape: [n_obs, n_classes])
    prob_matrix <- pred_clm$fit
    n_obs <- nrow(prob_matrix)
    n_class <- ncol(prob_matrix)

    n_mcmc <-ndraws
    predclm_samps <- array(0, dim = c(n_mcmc, n_test))
    for(i in 1:n_test){
      for(j in 1:n_mcmc){
        prob_j <- pred_clm$fit[i,]
        predclm_samps[j,i] <- sample(1:n_class, 1, prob=prob_j)
      }
    }

    # BART
    n_mcmc <- nrow(fit_bart$prob.test)
    predb_samps <- matrix(NA, nrow = n_mcmc, ncol = n_test)

    # Loop over test observations
    for (j in 1:n_test) {
      # Extract class probabilities for each MCMC sample for observation j
      idxs <- (1 + (j-1)*n_class):(j*n_class)
      probs_mat <- fit_bart$prob.test[, idxs]  # [ndraws × n_class]

      # Sample from categorical distribution
      for (i in 1:n_mcmc) {
        prob_i <- probs_mat[i, ]
        predi <- sample(1:n_class, 1, prob = prob_i)
        predb_samps[i, j] <- predi
      }
    }

    # Store raw predictive arrays + classes
    results[[fold]] <- list(
      truth = y_test,
      predc_array = predc_array,
      predc_class = predc_class,
      predp_array = predp_array,
      predp_class = predp_class,
      pred_clm_probs = pred_clm_array,
      pred_clm_class = pred_clm_class
    )

    # Compute EQWK
    q_khaos <- mean_qwk(predp_samps, y_test)
    q_cbass <- mean_qwk(predc_samps, y_test)
    q_clm   <- mean_qwk(predclm_samps, y_test)
    q_bart  <- mean_qwk(predb_samps, y_test)
    eq <- c(q_khaos, q_cbass, q_clm, q_bart)
    names(eq) <- c("khaos", "cbass", "clm", "bart")
    print(length(y_test))
    print(eq)


    results[[fold]]$eqwk <- eq
  }

  return(results = results)
}

# === Run analysis ===
data <- read.csv("~/Desktop/wine+quality/winequality-red.csv", sep = ";")
X <- data[ , -12]
X01 <- as.matrix(apply(X, 2, BASS:::scale_range), ncol = ncol(X))
y <- as.numeric(data[, 12]) - 2

cv_out5 <- crossval_compare(X01, y, k = 5, ndraws = 1000)
cv_out10 <- crossval_compare(X01, y, k = 10, ndraws = 1000)

mat5 <- do.call(rbind, lapply(cv_out5, function(l) l$eqwk))
mat10 <- do.call(rbind, lapply(cv_out10, function(l) l$eqwk))

colMeans(mat5)
colMeans(mat10)
save(cv_out5, cv_out10, file="data/cv_wine.rda")


dat <- duqling::get_emulation_data("Z_machine_max_vel1")
XX <- apply(dat$X, 2, BASS:::scale_range)
yy <- unlist(lapply(dat$y, function(yy) which.min(yy > quantile(dat$y, probs=(0:5)/6))))
fit <- ordinal_khaos(XX, yy)
plot(fit)

#cv_out5z <- crossval_compare(XX, yy, k = 5, ndraws = 1000)
cv_out10z <- crossval_compare(XX, yy, k = 10, ndraws = 1000)

#mat5z <- do.call(rbind, lapply(cv_out5z, function(l) l$eqwk))
mat10z <- do.call(rbind, lapply(cv_out10z, function(l) l$eqwk))

colMeans(mat5z)
colMeans(mat10z)
save(cv_out5z, cv_out10z, file="data/cv_zmach.rda")
