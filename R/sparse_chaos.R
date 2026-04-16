#' Bayesian Sparse Polynomial Chaos Expansion
#'
#' A Bayesian sparse polynomial chaos expansion (PCE) method based on the
#' forward-selection framework of Shao et al. (2017), with optional enrichment
#' strategies that adapt the candidate basis set across degree and interaction-order
#' expansions.
#'
#' @param X A data frame or matrix of predictors scaled to lie between 0 and 1.
#' @param y A numeric response vector of length \code{nrow(X)}.
#' @param degree Integer vector of length 3 giving the polynomial-degree schedule:
#'   \code{c(d_init, d_increment, d_max)}.
#' @param order Integer vector of length 3 giving the interaction-order schedule:
#'   \code{c(q_init, q_increment, q_max)}.
#' @param prior Numeric vector of prior hyperparameters
#'   \code{c(n0, v0, s0)}, where \code{n0} controls shrinkage of the regression
#'   coefficients, \code{v0} is the degrees-of-freedom parameter for the inverse-gamma
#'   prior on \eqn{\sigma^2}, and \code{s0} is the scale parameter for that prior.
#' @param lambda Nonnegative penalty parameter used in the KIC calculation. Larger
#'   values favor sparser models.
#' @param rho Basis functions are discarded after partial-correlation screening if
#'   their squared partial correlation is smaller than \code{rho}.
#' @param regularize Logical; if \code{TRUE}, a weighted LASSO pre-screen is used
#'   to reduce the candidate basis set before partial-correlation ranking.
#' @param max_basis Maximum number of candidate basis functions allowed in any
#'   fitting round. This acts as a hard cap on computational cost.
#' @param enrichment Either an integer in \code{0:4} specifying a fixed enrichment
#'   strategy, or the character string \code{"adaptive"}, which allows the strategy
#'   to become more restrictive when the next candidate set would otherwise be too
#'   large. See Details.
#' @param adaptive_ladder Numeric vector of candidate-size thresholds used only when
#'   \code{enrichment = "adaptive"}. The ladder should have length between 1 and 4,
#'   be sorted in increasing order, and satisfy \code{max(adaptive_ladder) <= max_basis}.
#'   Internally, the ladder is used to determine when the algorithm should fall back
#'   to a more restrictive enrichment strategy before attempting the next recursion
#'   step.
#' @param verbose Logical; if \code{TRUE}, progress messages are printed.
#'
#' @details
#' For fixed maximum degree and interaction order, the method constructs a candidate
#' library of polynomial basis functions, ranks them using marginal correlations and
#' sequential partial correlations, and evaluates nested models using the Kashyap
#' information criterion (KIC). The best-ranked model is retained at that stage.
#'
#' If the selected model contains a basis function at the current maximum degree or
#' interaction order, the algorithm attempts another fitting round with increased
#' degree and/or interaction order. In later rounds, the candidate basis set can be
#' rebuilt using one of several enrichment strategies rather than constructing the
#' full basis library from scratch.
#'
#' The \code{enrichment} argument controls how the candidate set is generated after
#' the first fitting round:
#' \describe{
#'   \item{\code{0} (local enrichment)}{
#'     Starts from the previously selected basis functions and applies local
#'     modifications to each term. These include deletion, simplification,
#'     complication, and promotion moves. This is typically the fastest strategy,
#'     but it is also the most local.}
#'   \item{\code{1} (active-variable rebuild)}{
#'     Rebuilds the full candidate library using only variables that were active
#'     in the previously selected model. This is computationally efficient, but
#'     inactive variables cannot re-enter once dropped.}
#'   \item{\code{2} (full active rebuild plus sparse inactive enrichment)}{
#'     Rebuilds the candidate library on the currently active variables, then
#'     enriches around the previously selected basis functions by allowing inactive
#'     variables to re-enter through complication and promotion moves. This is a
#'     compromise between speed and flexibility.}
#'   \item{\code{3} (full active rebuild plus inactive enrichment)}{
#'     Rebuilds the candidate library on the active variables, then enriches that
#'     larger active-variable library by allowing inactive variables to re-enter
#'     through complication and promotion moves. This is more expansive than
#'     \code{enrichment = 2}, but can be substantially more expensive.}
#'   \item{\code{4} (full rebuild)}{
#'     Rebuilds the full candidate library over all variables at each enrichment
#'     step. This is the most exhaustive strategy.}
#'   \item{\code{"adaptive"}}{
#'     Begins with a relatively expansive enrichment strategy and, before the next
#'     recursion step, uses an upper bound on the candidate-set size to decide
#'     whether a more restrictive strategy should be used instead. This recovers
#'     the computational benefit of avoiding construction of candidate sets that
#'     are known in advance to be too large.}
#' }
#'
#' When \code{enrichment = "adaptive"}, the function uses \code{adaptive_ladder}
#' together with an upper bound on the next candidate-set size. If the current
#' enrichment choice appears too expensive, the algorithm steps down to a more
#' restrictive enrichment rule before proceeding. If even the most restrictive
#' admissible next step is predicted to exceed the allowed budget, the current
#' fitted model is returned.
#'
#' Strategies \code{2}, \code{3}, and \code{"adaptive"} are designed to address a
#' limitation of active-variable-only rebuilds: variables that are inactive in one
#' round are allowed to re-enter the model in later rounds.
#'
#' @references
#' Shao, Q., Younes, A., Fahs, M., and Mara, T. A. (2017).
#' Bayesian sparse polynomial chaos expansion for global sensitivity analysis.
#' \emph{Computer Methods in Applied Mechanics and Engineering}, \strong{318},
#' 474--496.
#'
#' Rumsey, K. N., Francom, D., Gibson, G., Tucker, J. D., and Huerta, G. (2026).
#' Bayesian adaptive polynomial chaos expansions.
#' \emph{Stat}, \strong{15}(1), e70151.
#'
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391 * ((x[1] - 0.4) * (x[2] - 0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#'
#' fit <- sparse_khaos(X, y)
#'
#' @export
sparse_khaos <- function(X, y,
                        degree=c(2,2,16), order=c(1,1,5),
                        prior=c(1, 1, 1), lambda=1, rho=0, regularize=TRUE,
                        max_basis=1e6,
                        enrichment="adaptive",
                        adaptive_ladder = c(1e4, 1e5),
                        verbose=TRUE){
  if(max(X) > 1 || min(X) < 0){
    warning("X matrix should have values between 0 and 1.\n")
  }
  if(!is.matrix(X)){
    if(is.data.frame(X)){
      X <- as.matrix(X)
    }else{
      # Assuming it's a 1-d vector, try to coerce
      if(is.numeric(X)){
        X <- matrix(X, ncol=1)
      }else{
        stop("X should be a matrix, data.frame, or vector.")
      }
    }
  }
  n <- nrow(X)
  p <- ncol(X)

  d_curr <- degree[1]
  d_inc  <- degree[2]
  d_max  <- degree[3]

  o_curr <- order[1]
  o_inc  <- order[2]
  o_max  <- min(p, order[3])

  mu_y <- mean(y)
  sig_y <- stats::sd(y)
  y <- (y - mu_y)/sig_y

  if(enrichment == "adaptive"){
    # check ladder for validity
    if(length(adaptive_ladder) > 4 || length(adaptive_ladder) < 1){
      stop("adaptive_ladder has the wrong length, see documentation")
    }
    if(max(adaptive_ladder) > max_basis){
      stop("adaptive_ladder cannot exceed max_basis, see documentation")
    }
    enrichment <- length(adaptive_ladder)
    adaptive_ladder <- c(sort(adaptive_ladder), max_basis)
  }else{
    adaptive_ladder = NA
  }
  if(!(enrichment %in% 0:4)){
    stop("enrichment must be in 0:4")
  }

  A_curr <- NULL
  res <- sparse_khaos_wrapper(X, y, n, p, d_curr, d_inc, d_max, o_curr, o_inc, o_max, mu_y, sig_y, max_basis, prior, lambda, rho, regularize, verbose, enrichment, adaptive_ladder, A_curr)
  return(res)
}

sparse_khaos_wrapper <- function(X, y, n, p, d_curr, d_inc, d_max, o_curr, o_inc, o_max, mu_y, sig_y, max_basis, prior, lambda, rho, regularize, verbose, enrichment, adaptive_ladder, A_curr){
  # Create a list of p sequences from 1 to n
  if(verbose){
    cat("Starting model with max degree = ", d_curr, " and max order = ", o_curr, "\n", sep="")
  }

  if(is.null(A_curr)){
    A_num <- A_size(p, d_curr, o_curr)
    if(verbose) cat("\tFound ", A_num, " possible basis functions.\n", sep="")
    if(A_num > max_basis){
      stop("Too many basis functions. Increase max_basis or decrease initial degree/order. Consider dimension reduction approaches?")
    }
    if(verbose){
      if(A_num > 1000){
        t_est <- 4.258369e-7 * A_num ^ 1.781556
        if(t_est > 60){
          cat("\tComputing initial phi matrix (~", round(t_est/60, 2), " minutes)\n", sep = "")
        }else{
          cat("\tComputing initial phi matrix (~", round(t_est, 2), " seconds)\n", sep = "")
        }
      }else{
        cat("\tComputing initial phi matrix\n")
      }
    }
    A_set <- generate_A(p, d_curr, o_curr)
  }else{
    if(verbose){
      cat("\tBuilding enriched candidate set using strategy ", enrichment, "\n", sep = "")
    }
    A_set <- build_candidate_set(p, d_curr, o_curr, A_curr, enrichment, max_basis)
    A_num <- nrow(A_set)
    if(verbose) cat("\tFound ", A_num, " candidate basis functions after enrichment.\n", sep = "")
    if(A_num > max_basis){
      stop("Too many basis functions after enrichment. Increase max_basis or use a more restrictive enrichment strategy.")
    }
  }

  A_deg <- rowSums(A_set)
  A_ord <- rowSums(A_set > 0)
  N_alpha <- nrow(A_set)

  phi <- matrix(NA, nrow = n, ncol = N_alpha)
  rr <- rep(NA, N_alpha)
  for(i in 1:N_alpha){
    curr <- rep(1, n)
    for(j in 1:p){
      curr <- curr * ss_legendre_poly(X[,j], A_set[i, j])
    }
    phi[, i] <- curr
    rr[i] <- stats::cor(curr, y)
  }

  ord <- rev(order(rr^2))
  A_set <- A_set[ord, , drop = FALSE]
  A_deg <- A_deg[ord]
  A_ord <- A_ord[ord]
  phi <- phi[, ord, drop = FALSE]
  rr <- rr[ord]

  # Do LASSO?
  if(regularize){
    if(verbose) cat("\tRunning weighted LASSO...", sep="")
    lfit <- glmnet::glmnet(phi, y, penalty.factor=1-rr^2, dfmax=n-2)
    count_nonzero <- apply(lfit$beta, 2, function(bb) sum(bb != 0))
    valid_lambda <- which(count_nonzero < n)
    lambda_indx <- if(length(valid_lambda)) max(valid_lambda) else length(count_nonzero)
    alpha_indx <- which(lfit$beta[,lambda_indx] != 0)
    if(length(alpha_indx) == 0) alpha_indx <- 1
    N_alpha <- length(alpha_indx)
    A_set <- A_set[alpha_indx,,drop=FALSE]
    A_deg <- A_deg[alpha_indx]
    A_ord <- A_ord[alpha_indx]
    phi <- phi[,alpha_indx,drop=FALSE]
    rr <- rr[alpha_indx]
    if(verbose) cat(" Keeping ", N_alpha, " basis functions\n", sep="")
  }

  # Get partial correlation coefficients
  if(verbose) cat("\tComputing partial correlation coefficients\n")
  if(verbose) cat("\t\tFitting linear models: 0/", N_alpha, ", ", sep="")
  for(i in 2:N_alpha){
    if(N_alpha > 20 && ((i %% round(N_alpha/5)) == 0)){
      if(verbose) cat(i, "/", N_alpha, ", ",sep="")
    }
    eps_y <- stats::lm(y ~ phi[,1:(i-1)])$residuals
    eps_p <- stats::lm(phi[,i] ~ phi[,1:(i-1)])$residuals
    rr[i] <- stats::cor(eps_y, eps_p)
  }

  if(any(abs(rr) > 1)){
    rr <- rr/max(abs(rr))
  }
  ord <- rev(order(rr^2))
  A_set <- A_set[ord,,drop=FALSE]
  A_deg <- A_deg[ord]
  A_ord <- A_ord[ord]
  phi <- phi[,ord]
  rr <- rr[ord]

  # Get KIC for various models
  if(verbose) cat("\n\tRanking models based on KIC\n")
  # Add coefficient column in
  phi <- cbind(rep(1, n), phi)
  A_set <- rbind(rep(0, p), A_set)
  K_trunc <- 1 + max(c(0, which(rr^2 >= rho)))

  if(verbose){
    if(rho > 0) cat("\t\t Throwing out bases with low partial correlation. ", K_trunc, " remain.\n", sep="")
  }

  # Get g prior info
  Gdiag <- c(1, 1/(1 + A_ord*(A_ord + A_deg - 2)))
  n0 <- prior[1]
  G1 <- n*Gdiag/(n*Gdiag + n0)
  G2 <- n0/(n0 + n*Gdiag)

  # Get current BtB and BtBinv values (B is the phi matrix)
  B    <- phi[,1,drop=FALSE]
  BtB  <- crossprod(phi[,1,drop=FALSE])
  BtBi <- solve(BtB)
  Bty  <- crossprod(phi[,1,drop=FALSE], y)

  # Get MAP estimates
  a_mle  <- BtBi%*%Bty
  y_mle  <- B%*%a_mle
  s2_mle <- crossprod(y - y_mle)
  a_map  <- a_mle*G1[1]
  peff   <- 1
  s2_map <- (s2_mle + prior[3]^2 + G2[1]*crossprod(y_mle))/(n + 2 + prior[2] - 1)

  # Evaluate KIC for intercept-only model
  Sigma_prior <- sqrt(diag(BtBi) * Gdiag[1] * prior[3]^2 / prior[2] * n / n0)
  Sigma_post  <- (G1[1] * BtBi) * as.numeric(s2_mle)

  KIC <- rep(NA, K_trunc)
  KIC[1] <- -2 * sum(stats::dnorm(y, y_mle, sqrt(s2_map), log = TRUE)) -
    2 * sum(stats::dnorm(as.numeric(a_map), 0, Sigma_prior, log = TRUE)) -
    lambda * (peff + 1) * log(2 * pi) + log(det(Sigma_post))

  best <- list(
    k     = 1,
    KIC   = KIC[1],
    coeff = a_map,
    s2    = s2_map,
    BtB   = BtB,
    BtBi  = BtBi,
    G     = Gdiag[1]
  )
  if(K_trunc >= 2){
    for(k in 2:K_trunc){
      # Update BtB quantities
      B    <- phi[,1:k,drop=FALSE]
      bnew <- phi[,k,drop=FALSE]
      BtB_curr <- update_BtB(B, BtB, BtBi, bnew)
      BtB  <- BtB_curr$BtB
      BtBi <- BtB_curr$BtBi
      Bty  <- rbind(Bty, sum(bnew*y))

      # Get MAP estimates
      a_mle  <- BtBi%*%Bty
      y_mle  <- B%*%a_mle
      s2_mle <- crossprod(y - y_mle)
      a_map  <- a_mle*G1[1:k]
      y_reg  <- B%*%(a_mle * sqrt(G2[1:k]))
      peff   <- k
      s2_map <- (s2_mle + prior[3]^2 + crossprod(y_reg))/(n + 2 + prior[2] - peff)

      # Evaluate KIC
      yhat_k <- B%*%a_map
      Sigma_prior <- sqrt(diag(BtBi)*Gdiag[1:k]*prior[3]^2/prior[2]*n/n0) # Take the diagonal for simplicity
      Sigma_post  <- (G1[1:k] * BtBi) * as.numeric(s2_mle) # equivalent to (diag(G1) %*% BtBi)

      KIC[k] <- -2*sum(stats::dnorm(y, yhat_k, sqrt(s2_map), log=TRUE)) -
        2*sum(stats::dnorm(as.numeric(a_map), 0, Sigma_prior, log=TRUE)) -
        lambda*(k+1)*log(2*pi) + log(det(Sigma_post))

      if(KIC[k] <= best$KIC){
        best$k     <- k
        best$KIC   <- KIC[k]
        best$coeff <- a_map
        best$s2    <- s2_map
        best$BtB   <- BtB
        best$BtBi  <- BtBi
        best$G     <- Gdiag[1:k]
      }
    }
  }

  if(best$k <= 1){
    max_degree_curr <- 0
    max_order_curr <- 0
  } else {
    max_degree_curr <- max(A_deg[1:(best$k-1)])
    max_order_curr <- max(A_ord[1:(best$k-1)])
  }
  at_capacity_degree  <- (max_degree_curr == d_curr) #Note, Shao et al. suggest a +1 on the LHS of ineq.
  at_capacity_order   <- (max_order_curr == o_curr)
  over_max_degree     <- (d_curr >= d_max)
  over_max_order      <- (o_curr >= o_max)
  found_final_model   <- (!at_capacity_degree | over_max_degree) & (!at_capacity_order | over_max_order) # Neither degree nor order is at capacity (unless they're at the max)
  order_bigger_than_p <- (o_curr + o_inc > p) & (d_inc == 0)
  over_max_model      <- over_max_degree & over_max_order
  # Check whether the next candidate set would be too large.
  # For adaptive enrichment, try decreasing enrichment before giving up.
  d_next <- min(d_curr + d_inc, d_max)
  o_next <- min(o_curr + o_inc, o_max)

  if(best$k > 1){
    A_next <- A_set[2:best$k, , drop = FALSE]
  } else {
    A_next <- matrix(0, nrow = 0, ncol = p)
  }

  enrichment_next <- enrichment
  next_basis_too_big <- FALSE

  if(!all(is.na(adaptive_ladder))){
    p_active <- sum(colSums(A_next) > 0)
    A_next_num <- nrow(A_next)

    A_num_ub <- A_size_ub(p, d_next, o_next, p_active, A_next_num, enrichment_next)

    # adaptive_ladder has already had max_basis appended to the end
    # Example: c(1e4, 1e5, max_basis)
    # Then enrichment 2 -> 1e4, enrichment 1 -> 1e5, enrichment 0 -> max_basis
    while(enrichment_next > 0 &&
          A_num_ub > adaptive_ladder[length(adaptive_ladder) - enrichment_next]){
      enrichment_next <- enrichment_next - 1
      A_num_ub <- A_size_ub(p, d_next, o_next, p_active, A_next_num, enrichment_next)
    }

    if(A_num_ub > adaptive_ladder[length(adaptive_ladder) - enrichment_next]){
      next_basis_too_big <- TRUE
    }

    if(verbose && enrichment_next < enrichment){
      cat("\tAdaptive enrichment: reducing from ",
          enrichment, " to ", enrichment_next,
          " (new estimated size = ", format(A_num_ub, scientific = TRUE),
          ")\n", sep = "")
    }
  } else {
    p_active <- sum(colSums(A_next) > 0)
    A_next_num <- nrow(A_next)
    A_num_ub <- A_size_ub(p, d_next, o_next, p_active, A_next_num, enrichment_next)
    next_basis_too_big <- (A_num_ub > max_basis)
  }
  if(found_final_model | over_max_model | order_bigger_than_p | next_basis_too_big){
    if(verbose){
      if(next_basis_too_big) cat("Note: max_basis was reached. Consider dimension reduction?\n")
      if(over_max_model) cat("Note: Maximum degree and order were both reached. Consider increasing these values?\n")
    }
    obj <- list(B        = phi[, 1:best$k, drop = FALSE],
                nbasis   = best$k,
                vars     = A_set[1:best$k, , drop = FALSE],
                s2       = best$s2,
                beta_hat = best$coeff,

                X       = X,
                y       = y * sig_y + mu_y,

                mu_y    = mu_y,
                sigma_y = sig_y,

                BtB     = best$BtB,
                BtBi    = best$BtBi,
                G       = best$G,

                prior   = prior,

                KIC = KIC[1:best$k])

    class(obj) <- "sparse_khaos"
  }else{
    res <- sparse_khaos_wrapper(X, y, n, p, d_next, d_inc, d_max, o_next, o_inc, o_max, mu_y, sig_y, max_basis, prior, lambda, rho, regularize, verbose, enrichment_next, adaptive_ladder, A_next)
    return(res)
  }
  return(obj)
}

#' Predict Method for class sparse_khaos
#'
#' See \code{sparse_khaos()} for details.
#'
#' @param object An object returned by the \code{sparse_khaos()} function.
#' @param newdata A dataframe of the same dimension as the training data.
#' @param samples How many posterior samples should be taken at each test point? If 0 or FALSE, then the MAP estimate is returned.
#' @param nugget logical; should predictions include error? If FALSE, predictions will be for mean.
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
predict.sparse_khaos <- function(object, newdata=NULL, samples=1000, nugget=FALSE, ...){
  if(is.null(newdata)){
    newdata <- object$X
  }
  XX <- newdata
  n <- nrow(XX)
  p <- ncol(XX)
  N_alpha <- object$nbasis
  phi <- matrix(NA, nrow=n, ncol=N_alpha)
  for(i in 1:N_alpha){
    curr <- rep(1, n)
    for(j in 1:p){
      curr <- curr * ss_legendre_poly(XX[,j], object$vars[i,j])
    }
    phi[,i] <- curr
  }

  n0 <- object$prior[1]
  v0 <- object$prior[2]
  s0 <- object$prior[3]
  a_hat <- object$beta_hat
  peff <- length(a_hat)
  ntrain <- length(object$y)
  if(samples){
    pred <- matrix(NA, nrow=samples, ncol=n)
    for(i in 1:samples){
      shape <-  (ntrain + v0 - peff)/2
      sigma2 <- 1/stats::rgamma(1, shape=shape, rate=object$s2*shape)
      a_Sigma <- sigma2 * (object$G * object$BtBi)
      coeff <- stats::rnorm(peff, a_hat, sqrt(diag(a_Sigma)))
      y_hat <- phi%*%coeff
      if(nugget){
        y_hat <- y_hat + stats::rnorm(n, 0, sqrt(sigma2))
      }
      pred[i,] <- y_hat
    }
  }else{
    pred <- phi%*%object$beta_hat
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
plot.sparse_khaos <- function(x, ...){
  dots <- list(...)

  nugget <- if("nugget" %in% names(dots)) dots$nugget else FALSE
  samples <- if("samples" %in% names(dots)) dots$samples else 1000

  dots$nugget <- NULL
  dots$samples <- NULL

  pred <- predict(x, x$X, samples = samples, nugget = nugget)
  yhat <- colMeans(pred)

  do.call(graphics::plot, c(list(x = x$y, y = yhat, xlab="observed", ylab="predicted"), dots))
  graphics::abline(0, 1, lwd = 2, col = "orange")

  ci <- apply(pred, 2, function(yy) stats::quantile(yy, c(0.025, 0.975)))
  for(i in 1:ncol(ci)){
    graphics::segments(x$y[i], ci[1, i], x$y[i], ci[2, i])
  }
}

#' Print Method for class sparse_khaos
#'
#' See \code{sparse_khaos()} for details.
#'
#' @param x An object returned by the \code{sparse_khaos()} function.
#' @param ... additional arguments passed to or from other methods.
#' @details Print function for sparse_khaos object.
#' @references Shao, Q., Younes, A., Fahs, M., & Mara, T. A. (2017). Bayesian sparse polynomial chaos expansion for global sensitivity analysis. Computer Methods in Applied Mechanics and Engineering, 318, 474-496.
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- sparse_khaos(X, y)
#' print(fit)
#' @export
print.sparse_khaos <- function(x, ...){
  vars_used <- colSums(x$vars)
  names(vars_used) <- paste("V", seq_along(vars_used), sep="")

  bf_degrees <- table(rowSums(x$vars))
  names(bf_degrees) <- paste("deg", seq_along(bf_degrees), sep="")

  bf_orders <- table(rowSums(x$vars > 0))
  names(bf_orders) <- paste("ord", seq_along(bf_orders), sep="")

  cat("Fitted model contains", x$nbasis, "basis functions\n\n")
  cat("KIC =", round(min(x$KIC), 3), "\n\n")

  cat("Variables used:\n")
  print(vars_used)

  cat("\nBasis function degrees:\n")
  print(bf_degrees)

  cat("\nBasis function orders:\n")
  print(bf_orders)
}
