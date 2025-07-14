#' Bayesian Sparse Polynomial Chaos Expansion (Modified)
#'
#' A Bayesian implementation of sparse Polynomial Chaos Expansion (PCE), based on the framework of Shao et al. (2017), with modifications to improve flexibility, numerical stability, and support for Bayesian model evidence evaluation.
#'
#' @param X A data frame or matrix of input predictors, scaled to be between 0 and 1.
#' @param y A numeric response vector of length \code{nrow(X)}.
#' @param degree A vector of length 3 specifying the degree schedule: \code{c(d_init, d_increment, d_max)}. The degree is increased by \code{d_increment} if the current model includes a term of maximal degree, up to \code{d_max}.
#' @param order A vector of length 3 specifying the interaction order schedule: \code{c(q_init, q_increment, q_max)}.
#' @param sigma_prior A vector \code{c(v0, s0)} specifying the inverse-gamma prior on \eqn{sigma^2}: degrees of freedom (\code{v0}) and scale (\code{s0}). The prior is \eqn{sigma^2 \sim \text{Inv-Gamma}(v_0/2, v_0 s_0^2 / 2)}.
#' @param g_prior A vector \code{c(n0, m0, xi)} where \code{n0} is the prior strength on \eqn{g_0^2}, \code{m0} is the prior mean, and \code{xi} controls the penalty on complex basis functions in the modified \(g\)-prior. The prior is \eqn{g_0^2 \sim \text{Inv-Gamma}(n_0/2, n_0 m_0 / 2)}.
#' @param sM A numeric value controlling the model complexity prior: \eqn{p(\mathcal{M}) \propto |M|^{-s_M}}. Set \code{sM = 0} for a uniform prior over model size (default).
#' @param regularize Logical. If \code{TRUE}, applies LASSO-based pre-screening to reduce the size of the candidate basis set.
#' @param enrichment An integer between 0 and 3 specifying the enrichment strategy to use. See Details.
#' @param evidence Character string specifying how marginal likelihood evaluation should be handled. One of \code{"auto"}, \code{"fast"}, or \code{"full"}. See Details.
#' @param control  named \code{list()} of secondary tuning parameters, see
#'        \emph{Control list} below.
#' @param verbose Logical. If \code{TRUE}, progress updates and messages are printed.
#'
#' @section Control list:
#' \describe{
#'   \item{\code{grid_size}}{number of quantile points for numerically
#'         integrating over \(g_0^2\) (default 25).}
#'   \item{\code{max_basis_enrichment}}{threshold below which the
#_         full active-variable library is rebuilt each round
#'         (default 4\,000).}
#'   \item{\code{max_basis}}{hard upper bound on total candidate columns
#'         (default \(3X10^5\)).}
#'   \item{\code{patience}}{early-stop counter for the forward search
#'         (default \code{Inf}, i.e.\ evaluate all \(k\)).}
#'   \item{\code{orth_test}}{vector \code{c(L, significance_level)} for the automatic
#'         orthogonality diagnostic when \code{evidence = "auto"}
#'         (defaults \code{c(1000, 0.01)}).}
#' }
#' @details
#' The function builds the model incrementally, starting with polynomial degree \code{degree[1]} and interaction order \code{order[1]}. If the current model includes any terms of maximal degree or order, the values are increased (up to \code{degree[3]} and \code{order[3]}, respectively), and the process restarts.
#'
#' The modified \(g\)-prior used here includes an additional shrinkage term controlled by \code{xi}, which penalizes complex basis functions based on their total degree and interaction order. Priors over \(sigma^2\) and \(g_0^2\) are specified using interpretable parameters derived from standard inverse-gamma formulations.
#'
#' The \code{evidence} argument controls how marginal likelihoods are evaluated during model selection. If \code{"auto"}, the algorithm performs a one-time diagnostic of the input distribution using Kolmogorov–Smirnov and simulation-based Frobenius norm tests to determine whether the orthogonality approximation can be safely used. The options \code{"fast"} and \code{"full"} override this and force either the approximation or the full method, respectively.
#'
#' Five enrichment strategies are implemented (from fastest to slowest):
#' \describe{
#'   \item{0}{Local enrichment only, using all four enrichment operators.}
#'   \item{1}{Original Shao method (discard inactives, full rebuild).}
#'   \item{2}{Build all active-variable terms and enrich only around previously accepted basis functions.}
#'   \item{3}{Build all active-variable terms and enrich with inactives.}
#'   \item{4}{Full rebuild}
#' }
#' Note that the sull rebuild is always used unless the full number of candidates exceeds `code{max_basis_enrichment}`.
#'
#' @references
#' Shao, Q., Younes, A., Fahs, M., & Mara, T. A. (2017). Bayesian sparse polynomial chaos expansion for global sensitivity analysis. \emph{Computer Methods in Applied Mechanics and Engineering}, 318, 474–496.
#'
#' Robert, C. P. (2007). \emph{The Bayesian Choice: From Decision-Theoretic Foundations to Computational Implementation}. Springer.
#'
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391 * ((x[1] - 0.4) * (x[2] - 0.6) + 0.36)
#' y <- apply(X, 1, f) + rnorm(100, 0, 0.1)
#' fit <- sparse_khaos2(X, y)
#'
#' @export
sparse_khaos2 <- function(X, y,
                         degree=c(2,2,16), order=c(1,1,4),
                         sigma_prior = c(1,1), g_prior = c(1, 1, 1), sM=0,
                         regularize=TRUE,
                         enrichment=2,
                         evidence="auto",
                         control = list(grid_size = 25,
                                        max_basis_enrichment = 4000,
                                        max_basis = 3e5,
                                        patience = 3,
                                        orth_test = c(1000, 0.01)),
                         verbose=TRUE){

  # Extract control parameters
  grid_size <- control$grid_size
  max_basis <- control$max_basis
  max_basis_enrichment <- control$max_basis_enrichment
  patience <- control$patience
  orth_test <- control$orth_test

  # Check X matrix
  if(max(X) > 1 || min(X) < 0){
    warning("X matrix should have values between 0 and 1.\n")
  }
  if(!is.matrix(X)){
    X <- matrix(X, ncol=1)
  }
  n <- nrow(X)
  p <- ncol(X)

  # Test for inputs independent and U(0, 1) if evidence = "auto"
  if(evidence == "auto"){
    if(verbose){
      cat("Testing if orthogonality assumption is reasonable.\n")
    }
    if(test_iid_unif(X)){
      evidence <- "fast"
    }else{
      evidence <- "full"
    }
  }
  if(!(evidence %in% c("full", "fast"))){
    stop("evidence must be one of auto, full, or fast.")
  }

  d_curr <- degree[1]
  d_inc  <- degree[2]
  d_max  <- degree[3]

  o_curr <- order[1]
  o_inc  <- order[2]
  o_max  <- min(p, order[3])

  mu_y <- mean(y)
  sig_y <- stats::sd(y)
  y <- (y - mu_y)/sig_y

  A_curr <- NULL

  res <-    sparse_khaos_wrapper2(X, y, n, p, d_curr, d_inc, d_max, o_curr, o_inc, o_max, mu_y, sig_y, max_basis, sigma_prior, g_prior, sM, grid_size, patience, enrichment, max_basis_enrichment, evidence, regularize, verbose, A_curr)
  return(res)
}



sparse_khaos_wrapper2 <- function(X, y, n, p, d_curr, d_inc, d_max, o_curr, o_inc, o_max, mu_y, sig_y, max_basis, sigma_prior, g_prior, sM, grid_size, patience, enrichment, max_basis_enrichment, evidence, regularize, verbose, A_curr){
  # Create a list of p sequences from 1 to n
  if(verbose){
    cat("Starting model with max degree = ", d_curr, " and max order = ", o_curr, "\n", sep="")
  }
  A_num <- A_size(p, d_curr, o_curr)
  if(verbose) cat("\tFound ", A_num, " possible basis functions.\n", sep="")
  if(A_num > max_basis){
    stop("Too many basis functions. Increase max_basis or decrease initial degree/order. Consider dimension reduction approaches?")
  }
  if(A_num > max_basis_enrichment){
    if(is.null(A_curr)){
      stop("max_basis_enrichment is smaller than the size of the initial candidate set (size = ", A_num, ").\n",
           "Increase max_basis_enrichment or start with smaller (d, q) settings.")
    }

    # Determine how to do enrichment based on `enrichment` argument
    enrichment_func <- get(paste0("enrichment_", enrichment))
    A_set <- enrichment_func(p, d_curr, o_curr, A_curr)
  }else{
    if(verbose){
      if(A_num > 1000){
        t_est <- 4.258369e-7 * A_num ^ 1.781556
        if(t_est > 60){
          cat("\tComputing initial phi matrix (~", round(t_est/60, 2), " minutes)\n",sep="")
        }else{
          cat("\tComputing initial phi matrix (~", round(t_est, 2), " seconds)\n",sep="")
        }
      }else{
        cat("\tComputing initial phi matrix\n")
      }
    }
    A_set <- generate_A(p, d_curr, o_curr)
  }

  A_deg <- apply(A_set, 1, sum)
  A_ord <- apply(A_set, 1, function(aa) sum(aa > 0))
  N_alpha <- nrow(A_set)
  if(N_alpha != A_num) warning("This warning shouldn't happen, but it's also not that big of a deal. (:\n")
  phi <- matrix(NA, nrow=n, ncol=N_alpha)
  rr <- rep(NA, N_alpha)
  for(i in 1:N_alpha){
    curr <- rep(1, n)
    for(j in 1:p){
      curr <- curr * ss_legendre_poly(X[,j], A_set[i,j])
    }
    phi[,i] <- curr
    rr[i] <- stats::cor(curr, y)
  }
  ord <- rev(order(rr^2))
  A_set <- A_set[ord,]
  A_deg <- A_deg[ord]
  A_ord <- A_ord[ord]
  phi <- phi[,ord]
  rr <- rr[ord]

  ### MODEL SELECTION
  # Add coefficient column in
  phi <- cbind(rep(1, n), phi)
  A_set <- rbind(rep(0, p), A_set)

  # Get K_trunc
  #rho <- 1e-7
  rho <- 0
  K_trunc <- max(which(rr^2 >= rho)) + 1

  # Get parameters
  nu0   <- sigma_prior[1]               # df
  s0_sq <- sigma_prior[2]^2             # scale²
  a_sig <- nu0 / 2
  b_sig <- nu0 * s0_sq / 2

  n0    <- g_prior[1]                   # prior strength
  m0    <- g_prior[2]                   # prior mean
  xi    <- g_prior[3]                   # complexity penalty (ξ)
  a_g   <- n0 / 2
  b_g   <- n0 * m0 / 2

  # Handle the g stuff
  g_diag <- c(1, (1 + A_ord * (A_ord + A_deg - 2))^(-xi/2))
  g0_grid  <- 1 / stats::qgamma((1:grid_size)/(grid_size+1), a_g, b_g)  # Inv-Gamma quantiles
  log_mean_exp <- function(z) {               # stable log-mean-exp
    m <- max(z);  m + log(mean(exp(z - m)))
  }

  # Initialize BtB etc (for intercept-only model)
  n     <- length(y)
  yty   <- sum(y^2)
  BtB   <- matrix(n, 1, 1)
  BtBi  <- 1 / n
  Bty   <- matrix(sum(y), 1, 1)
  v     <- Bty #/ n

  # Bookkeeping
  best <- list(k=0, lpost=-Inf)
  LPOST <- rep(NA, K_trunc)
  momentum <- 0

  # Start loop over nested models
  for(k in 1:K_trunc){
    if (k > 1){
      ## rank-1 update BtB / BtBi / Bty
      bnew  <- phi[, k, drop = FALSE]
      up    <- update_BtB(phi[, 1:k, drop = FALSE], BtB, BtBi, bnew)
      BtB   <- up$BtB
      BtBi  <- up$BtBi
      Bty   <- rbind(Bty, crossprod(bnew, y))
      v     <- c(v, crossprod(bnew, y) ) #/ n)
    }

    # Integrate log evidence over g0^2 prior
    logev_vec <- vapply(g0_grid, function(g0_sq)
      log_evidence(evidence,
                   yty = yty, v = v,
                   BtB = BtB, BtBi = BtBi, Bty = Bty,
                   g_diag = g_diag[1:k], g0_sq = g0_sq,
                   nu0 = nu0, s0_sq = s0_sq, n = n),
      numeric(1))
    logev <- log_mean_exp(logev_vec)


    # Get log posterior
    LPOST[k] <- logev - sM * log(k) # k includes intercept

    # Update best model if improved
    if (LPOST[k] > best$lpost) {
      best$k     <- k
      best$lpost <- LPOST[k]
      best$BtB   <- BtB
      best$BtBi  <- BtBi
      best$G     <- g_diag[1:k]
      momentum <- 0
    }else{
      momentum <- momentum + 1
    }
    if(momentum > patience){
      break
    }
  }

  # Next pass logic
  max_degree_curr    <- max(A_deg[1:(best$k - 1)])
  max_order_curr     <- max(A_ord[1:(best$k - 1)])
  at_capacity_degree <- (max_degree_curr == d_curr)
  at_capacity_order  <- (max_order_curr == o_curr)
  over_max_degree    <- (d_curr >= d_max)
  over_max_order     <- (o_curr >= o_max)

  found_final_model <- (!at_capacity_degree | over_max_degree) &
    (!at_capacity_order  | over_max_order)

  next_basis_too_big <- (A_size(p, d_curr + d_inc,
                                o_curr + o_inc) > max_basis)

  over_max_model <- over_max_degree & over_max_order

  if (found_final_model | next_basis_too_big | over_max_model){
    if (next_basis_too_big)
      cat("Note: max_basis was reached. Consider dimension reduction?\n")
    if (over_max_model)
      cat("Note: Maximum degree and order reached. Consider increasing them.\n")

    obj <- list(B       = phi[, 1:best$k, drop = FALSE],
                nbasis  = best$k,
                vars    = A_set[1:best$k, , drop = FALSE],

                X       = X,
                y       = y * sig_y + mu_y,

                mu_y    = mu_y,
                sigma_y = sig_y,

                BtB     = best$BtB,
                BtBi    = best$BtBi,
                G       = best$G,

                prior   = list(sigma = sigma_prior,
                               g     = g_prior,
                               sM    = sM),

                lpost = LPOST[1:k]) # k = best$k + momentum ?

    class(obj) <- "sparse_khaos2"
    return(obj)

  } else {
    d_next <- min(d_curr + d_inc, d_max)
    o_next <- min(o_curr + o_inc, o_max)
    res <- sparse_khaos_wrapper2(X, y, n, p,
                                 d_next, d_inc, d_max,
                                 o_next, o_inc, o_max,
                                 mu_y, sig_y,
                                 max_basis,
                                 sigma_prior, g_prior, sM,
                                 grid_size, patience, enrichment,
                                 max_basis_enrichment,
                                 evidence, regularize, verbose,
                                 A_set[1:best$k, , drop = FALSE])

    return(res)
  }

}



#' Predict Method for class sparse_khaos2
#'
#' See \code{sparse_khaos2()} for details.
#'
#' @param object An object returned by the \code{sparse_khaos()} function.
#' @param newdata A dataframe of the same dimension as the training data.
#' @param samples How many posterior samples should be taken at each test point? If 0 or FALSE, then the MAP estimate is returned.
#' @param ... Additional arguments to predict
#' @details Predict function for sparse_khaos object.
#' @references Shao, Q., Younes, A., Fahs, M., & Mara, T. A. (2017). Bayesian sparse polynomial chaos expansion for global sensitivity analysis. Computer Methods in Applied Mechanics and Engineering, 318, 474-496.
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- sparse_khaos(X, y)
#' predict(fit)
#'
#' @export
predict.sparse_khaos2 <- function(object, newdata=NULL, samples=1000, ...){
  if(is.null(newdata)){
    newdata <- object$X
  }
  XX <- newdata
  n <- nrow(XX)
  p <- ncol(XX)
  N_alpha <- nrow(object$vars)
  phi <- matrix(NA, nrow=n, ncol=N_alpha)
  for(i in 1:N_alpha){
    curr <- rep(1, n)
    for(j in 1:p){
      curr <- curr * ss_legendre_poly(XX[,j], object$vars[i,j])
    }
    phi[,i] <- curr
  }

  v0_sigma <- object$prior$sigma[1]
  s0_sigma <- object$prior$sigma[2]

  n0_g <- object$prior$g[1]
  m0_g <- object$prior$g[2]
  xi_g <- object$prior$g[3]

  G <- matrix(1, nrow=N_alpha, ncol=N_alpha)

  n0 <- object$prior[1]
  v0 <- object$prior[2]
  s0 <- object$prior[3]
  ntrain <- length(object$y)
  if(samples){
    pred <- matrix(NA, nrow=samples, ncol=n)
    for(i in 1:samples){
      shape <-  (ntrain+v0-p)/2
      sigma2 <- 1/stats::rgamma(1, shape, object$s2*shape)
      a_Sigma <- sigma2 * (object$G * object$BtBi)
      a_hat <- object$coeff
      coeff <- stats::rnorm(a_hat, a_hat, diag(a_Sigma))
      #noise <- sqrt(1/stats::rgamma(1, (n+2)/2, scale=2/(n*object$s2)))
      y_hat <- phi%*%coeff + stats::rnorm(n, 0, sqrt(sigma2))
      pred[i,] <- y_hat
    }
  }else{
    pred <- phi%*%object$coeff
  }
  pred <- object$mu_y + object$sigma_y * pred
  return(pred)
}

#' Plot Method for class sparse_khaos
#'
#' See \code{sparse_khaos()} for details.
#'
#' @param x An object returned by the \code{sparse_khaos()} function.
#' @param ... additional arguments passed to \code{plot}
#' @details Plot function for sparse_khaos object.
#' @references Shao, Q., Younes, A., Fahs, M., & Mara, T. A. (2017). Bayesian sparse polynomial chaos expansion for global sensitivity analysis. Computer Methods in Applied Mechanics and Engineering, 318, 474-496.
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- sparse_khaos(X, y)
#' plot(fit)
#' @export
plot.sparse_khaos2 <- function(x, ...){
  pred <- stats::predict(x, x$X, samples=1000)
  yhat <- colMeans(pred)
  plot(x$y, yhat, ...)
  graphics::abline(0, 1, lwd=2, col='orange')

  ci <- apply(pred, 2, function(yy) stats::quantile(yy, c(0.025, 0.975)))
  for(i in 1:ncol(ci)){
    graphics::segments(x$y[i], ci[1,i], x$y[i], ci[2,i])
  }
}
