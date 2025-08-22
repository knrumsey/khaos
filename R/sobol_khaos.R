#' Compute Sobol Sensitivity Indices for khaos Models
#'
#' Computes first-order, total-effect, and full interaction Sobol indices
#' from an adaptive KHAOS model fit. Returns MCMC samples of sensitivity
#' indices, including all active interaction terms used in the model.
#'
#' @param obj An object of class \code{"adaptive_khaos"} returned by \code{adaptive_khaos()} or of class \code{"sparse_khaos"} returned by \code{sparse_khaos()}.
#' @param plot_it Logical; if `TRUE`, boxplots of the sensitivity indices are shown.
#' @param samples How many samples from the posterior (only used for sparse class).
#' @return A list with the following components:
#' \describe{
#'   \item{S}{A data frame of Sobol indices for all main effects and interactions (MCMC samples).}
#'   \item{T}{A data frame of total-effect Sobol indices (MCMC samples).}
#'   \item{leftover}{A numeric vector of unexplained variance proportions (residuals).}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- adaptive_khaos(X, y)
#' sob <- sobol_khaos(fit)
#' head(sob$S)
#' }
#'
#' @export
sobol_khaos <- function(obj, plot_it = TRUE, samples=1000) {
  if(class(obj) == "sparse_khaos"){
    return(sobol_sparse_khaos(obj, samples, plot_it))
  }else{
    if(!(class(obj) %in% c("ordinal_khaos", "adaptive_khaos"))){
      stop("obj must be an object of class sparse_khaos or adaptive_khaos.")
    }
  }

  beta <- obj$beta
  vars <- obj$vars
  degs <- obj$degs
  nint <- obj$nint
  X    <- obj$X
  if(class(obj) == "ordinal_khaos"){
    s2 <- obj$s2z
  }else{
    s2 <- obj$s2
  }

  n_iter <- dim(beta)[1]
  M <- dim(beta)[2] - 1
  p <- ncol(X)
  var_names <- paste0("x", seq_len(p))

  first_order_mat <- matrix(NA, nrow = n_iter, ncol = p)
  total_effect_mat <- matrix(NA, nrow = n_iter, ncol = p)
  leftover <- rep(NA, n_iter)

  interaction_list <- vector("list", n_iter)
  all_labels_set <- new.env(hash = TRUE)

  for (iter in 1:n_iter) {
    b <- beta[iter, -1]  # exclude intercept
    b2 <- b^2
    explained <- sum(b2, na.rm = TRUE)
    residual <- s2[iter]
    total_var <- explained + residual
    if (total_var < 1e-12 || is.na(total_var)) next

    alpha_mat <- matrix(0, nrow = M, ncol = p)
    interact_map <- list()

    for (m in 1:M) {
      k <- nint[iter, m]
      if (is.na(k)) next
      if (k == 0) next
      idxs <- vars[iter, m, 1:k]
      degrees <- degs[iter, m, 1:k]
      alpha_mat[m, idxs] <- degrees

      idxs_sorted <- sort(idxs)
      label <- paste0("x", idxs_sorted, collapse = ":")

      # Store label globally
      all_labels_set[[label]] <- TRUE

      if (label %in% names(interact_map)) {
        interact_map[[label]] <- interact_map[[label]] + b2[m]
      } else {
        interact_map[[label]] <- b2[m]
      }
    }

    # Normalize interaction values
    interaction_list[[iter]] <- sapply(interact_map, function(val) val / total_var)

    for (j in 1:p) {
      only_j <- which(alpha_mat[, j] > 0 & rowSums(alpha_mat[, -j, drop = FALSE]) == 0)
      any_j  <- which(alpha_mat[, j] > 0)
      first_order_mat[iter, j] <- sum(b2[only_j], na.rm = TRUE) / total_var
      total_effect_mat[iter, j] <- sum(b2[any_j], na.rm = TRUE) / total_var
    }

    leftover[iter] <- residual / total_var
  }

  # Construct full set of labels
  all_labels <- ls(all_labels_set)
  # Sort by interaction order, then alphabetically
  label_order <- order(sapply(strsplit(all_labels, ":"), length), all_labels)
  all_labels <- all_labels[label_order]

  # Construct full matrix for S
  S_mat <- matrix(0, nrow = n_iter, ncol = length(all_labels))
  colnames(S_mat) <- all_labels

  for (i in 1:n_iter) {
    terms <- interaction_list[[i]]
    if (length(terms) > 0) {
      matched <- match(names(terms), all_labels)
      S_mat[i, matched] <- terms
    }
  }

  # Output
  colnames(first_order_mat) <- var_names
  colnames(total_effect_mat) <- var_names

  out <- list(
    S = as.data.frame(S_mat),
    T = as.data.frame(total_effect_mat),
    leftover = leftover
  )

  if (plot_it) {
    op <- par(mfrow=c(1,2))
    # Add plotting later
    tmp <- cbind(out$S, leftover)
    boxplot(tmp, main="Sensitivity")
    tmp <- cbind(out$T, leftover)
    boxplot(tmp, main="Total Sensitivity")
    par(op)
  }

  return(out)
}


sobol_sparse_khaos <- function(obj, samples = 1000, plot_it = TRUE) {
  X <- obj$X
  y <- obj$y
  beta_hat <- obj$beta_hat
  G <- obj$G
  BtBi <- obj$BtBi
  s2 <- obj$s2
  prior <- obj$prior
  nbasis <- obj$nbasis
  vars <- obj$vars
  p <- ncol(X)
  n <- nrow(X)

  # Generate basis matrix phi [n x nbasis]
  phi <- matrix(NA, nrow = n, ncol = nbasis)
  for (i in 1:nbasis) {
    curr <- rep(1, n)
    for (j in 1:p) {
      curr <- curr * khaos:::ss_legendre_poly(X[, j], vars[i, j])
    }
    phi[, i] <- curr
  }

  # Setup for sampling from posterior over coefficients
  v0 <- prior[2]
  ntrain <- length(y)
  shape <- (ntrain + v0 - p) / 2
  coeff_samples <- matrix(NA, nrow = samples, ncol = nbasis)
  sigma2_samples <- numeric(samples)

  for (i in 1:samples) {
    sigma2 <- 1 / stats::rgamma(1, shape, s2 * shape)
    sigma2_samples[i] <- sigma2
    cov_mat <- sigma2 * (G * BtBi)
    coeff_samples[i, ] <- stats::rnorm(nbasis, mean = beta_hat, sd = sqrt(diag(cov_mat)))
  }

  # Prepare for Sobol decomposition
  var_names <- paste0("x", seq_len(p))
  first_order_mat <- matrix(NA, nrow = samples, ncol = p)
  total_effect_mat <- matrix(NA, nrow = samples, ncol = p)
  leftover <- numeric(samples)
  interaction_list <- vector("list", samples)
  all_labels_set <- new.env(hash = TRUE)

  for (iter in 1:samples) {
    b <- coeff_samples[iter, ]
    b2 <- b^2
    explained <- sum(b2)
    residual <- sigma2_samples[iter]
    total_var <- explained + residual
    if (total_var < 1e-12 || is.na(total_var)) next

    alpha_mat <- vars
    alpha_mat[is.na(alpha_mat)] <- 0

    interact_map <- list()

    for (m in 1:nbasis) {
      idxs <- which(alpha_mat[m, ] > 0)
      if (length(idxs) == 0) next

      label <- paste0("x", sort(idxs), collapse = ":")
      all_labels_set[[label]] <- TRUE

      if (label %in% names(interact_map)) {
        interact_map[[label]] <- interact_map[[label]] + b2[m]
      } else {
        interact_map[[label]] <- b2[m]
      }
    }

    interaction_list[[iter]] <- sapply(interact_map, function(val) val / total_var)

    for (j in 1:p) {
      only_j <- which(alpha_mat[, j] > 0 & rowSums(alpha_mat[, -j, drop = FALSE]) == 0)
      any_j  <- which(alpha_mat[, j] > 0)
      first_order_mat[iter, j] <- sum(b2[only_j]) / total_var
      total_effect_mat[iter, j] <- sum(b2[any_j]) / total_var
    }

    leftover[iter] <- residual / total_var
  }

  all_labels <- ls(all_labels_set)
  label_order <- order(sapply(strsplit(all_labels, ":"), length), all_labels)
  all_labels <- all_labels[label_order]

  S_mat <- matrix(0, nrow = samples, ncol = length(all_labels))
  colnames(S_mat) <- all_labels

  for (i in 1:samples) {
    terms <- interaction_list[[i]]
    if (length(terms) > 0) {
      matched <- match(names(terms), all_labels)
      S_mat[i, matched] <- terms
    }
  }

  colnames(first_order_mat) <- var_names
  colnames(total_effect_mat) <- var_names

  out <- list(
    S = as.data.frame(S_mat),
    T = as.data.frame(total_effect_mat),
    leftover = leftover
  )

  if (plot_it) {
    op <- par(mfrow = c(1, 2))
    tmp <- cbind(out$S, leftover)
    boxplot(tmp, main = "Sensitivity")
    tmp <- cbind(out$T, leftover)
    boxplot(tmp, main = "Total Sensitivity")
    par(op)
  }

  return(out)
}
