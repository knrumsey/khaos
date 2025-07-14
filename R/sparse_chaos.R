#' Bayesian Sparse Polynomial Chaos Expansion
#'
#' The emulation approach of Shao et al. (2017), with some minor modifications (see README for details).
#'
#' @param X A dataframe or matrix of predictors scaled to be between 0 and 1
#' @param y a reponse vector
#' @param degree Maximum polynomial degree parameters. c(Initial, Increment, Maximum)
#' @param order Maximum order of interaction parameters. c(Initial, Increment, Maximum)
#' @param prior A vector with prior values: n0 (strength of prior for coefficients), v0 (degrees of freedom for IG prior for sigma2), s0 (scale for IG prior for sigma2).
#' @param lambda A parameter fed to the KIC calculations. Larger values will lead to sparser models.
#' @param rho Basis functions are only kept if their square partial correlation coefficient is bigger than rho.
#' @param regularize Logical. Should LASSO be used to pre-sparsify the library of basis functions?
#' @param max_basis The maximum number of candidate basis functions to consider (bounds the computation time)
#' @param verbose Logical. Should progress information be printed?
#' @details Implements a modified version of the Bayesian sparse PCE method described by Shao et al. (2017).
#' Degree (equivalently for order) is incremented by \code{degree[2]} if the current model contains any terms of maximal degree, and the process is started again from the top. Cannot be incremented above \code{degree[3]}.
#' Fitting process terminates when either (i) max_basis is exceeded or (ii) current model fit does not contain any terms with degree and order not equal to the current maximum value.
#' @references Shao, Q., Younes, A., Fahs, M., & Mara, T. A. (2017). Bayesian sparse polynomial chaos expansion for global sensitivity analysis. Computer Methods in Applied Mechanics and Engineering, 318, 474-496.
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- sparse_khaos(X, y)
#' @export
sparse_khaos <- function(X, y,
                        degree=c(2,2,16), order=c(1,1,4),
                        prior=c(1, 1, 1), lambda=1, rho=0, regularize=TRUE,
                        max_basis=3e5,
                        verbose=TRUE){
  if(max(X) > 1 || min(X) < 0){
    warning("X matrix should have values between 0 and 1.\n")
  }
  if(!is.matrix(X)){
    X <- matrix(X, ncol=1)
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

  res <- sparse_khaos_wrapper(X, y, n, p, d_curr, d_inc, d_max, o_curr, o_inc, o_max, mu_y, sig_y, max_basis, prior, lambda, rho, regularize, verbose)
  return(res)
}

sparse_khaos_wrapper <- function(X, y, n, p, d_curr, d_inc, d_max, o_curr, o_inc, o_max, mu_y, sig_y, max_basis, prior, lambda, rho, regularize, verbose){
  # Create a list of p sequences from 1 to n
  if(verbose){
    cat("Starting model with max degree = ", d_curr, " and max order = ", o_curr, "\n", sep="")
  }
  A_num <- A_size(p, d_curr, o_curr)
  if(verbose) cat("\tFound ", A_num, " possible basis functions.\n", sep="")
  if(A_num > max_basis){
    stop("Too many basis functions. Increase max_basis or decrease initial degree/order. Consider dimension reduction approaches?")
  }
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

  # Do LASSO?
  if(regularize){
    if(verbose) cat("\tRunning weighted LASSO...", sep="")
    lfit <- glmnet::glmnet(phi, y, penalty.factor=1-rr^2, dfmax=n-2)
    count_nonzero <- apply(lfit$beta, 2, function(bb) sum(bb != 0))
    lambda_indx <- max(which(count_nonzero < n))
    lambda_use <- lfit$lambda[lambda_indx]
    alpha_indx <- which(lfit$beta[,lambda_indx] != 0)
    if(length(alpha_indx) == 0) alpha_indx <- 1
    N_alpha <- length(alpha_indx)
    A_set <- A_set[alpha_indx,]
    A_deg <- A_deg[alpha_indx]
    A_ord <- A_ord[alpha_indx]
    phi <- phi[,alpha_indx]
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
  A_set <- A_set[ord,]
  A_deg <- A_deg[ord]
  A_ord <- A_ord[ord]
  phi <- phi[,ord]
  rr <- rr[ord]

  # Get KIC for various models
  if(verbose) cat("\n\tRanking models based on KIC\n")
  # Add coefficient column in
  phi <- cbind(rep(1, n), phi)
  A_set <- rbind(rep(0, p), A_set)
  K_trunc <- max(which(rr^2 >= rho)) + 1
  if(verbose){
    if(rho > 0) cat("\t\t Throwing out bases with low partial correlation. ", K_trunc, " remain.\n", sep="")
  }
  KIC <- rep(NA, K_trunc)
  KIC[1] <- Inf

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
  s2_map <- (s2_mle + prior[3]^2 + G2[1]*crossprod(y_mle))/(n + 2 + prior[2] - p)

  best <- list(k=0, KIC=Inf, map=NULL)
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
    s2_map <- (s2_mle + prior[3]^2 + crossprod(y_reg))/(n + 2 + prior[2] - p)

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

  max_degree_curr     <- max(A_deg[1:(best$k-1)])
  max_order_curr      <- max(A_ord[1:(best$k-1)])
  at_capacity_degree  <- (max_degree_curr == d_curr) #Note, Shao et al. suggest a +1 on the LHS of ineq.
  at_capacity_order   <- (max_order_curr == o_curr)
  over_max_degree     <- (d_curr >= d_max)
  over_max_order      <- (o_curr >= o_max)
  found_final_model   <- (!at_capacity_degree | over_max_degree) & (!at_capacity_order | over_max_order) # Neither degree nor order is at capacity (unless they're at the max)
  next_basis_too_big  <- (A_size(p, d_curr+d_inc, o_curr+o_inc) > max_basis)
  over_max_model      <- over_max_degree & over_max_order
  if(found_final_model | next_basis_too_big | over_max_model){
    #if(next_degree_too_big) cat("Note: max_degree was reached.\n")
    if(next_basis_too_big) cat("Note: max_basis was reached. Consider dimension reduction?\n")
    if(over_max_model) cat("Note: Maximum degree and order were both reached. Consider increasing these values?\n")
    #obj <- list(coeff=best$coeff, s2=best$s2,
    #            phi=phi[,1:best$k,drop=FALSE],
    #            vars=A_set[1:best$k,,drop=FALSE],
    #            mu_y=mu_y, sigma_y=sig_y, KIC=best$KIC, X=X, y=y*sig_y + mu_y,
    #            BtB=best$BtB, BtBi=best$BtBi, prior=prior, G=best$G)


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
    d_next <- min(d_curr + d_inc, d_max)
    o_next <- min(o_curr + o_inc, o_max)
    res <- sparse_khaos_wrapper(X, y, n, p, d_next, d_inc, d_max, o_next, o_inc, o_max, mu_y, sig_y, max_basis, prior, lambda, rho, regularize, verbose)
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
predict.sparse_khaos <- function(object, newdata=NULL, samples=1000, ...){
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
plot.sparse_khaos <- function(x, ...){
  pred <- predict(x, x$X, samples=1000)
  yhat <- colMeans(pred)
  plot(x$y, yhat, ...)
  graphics::abline(0, 1, lwd=2, col='orange')

  ci <- apply(pred, 2, function(yy) stats::quantile(yy, c(0.025, 0.975)))
  for(i in 1:ncol(ci)){
    graphics::segments(x$y[i], ci[1,i], x$y[i], ci[2,i])
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
  vars_used <- colSums(fit$vars)
  names(vars_used) <- paste("V", seq_along(vars_used), sep="")

  bf_degrees <- table(rowSums(fit$vars))
  names(bf_degrees) <- paste("deg", seq_along(bf_degrees), sep="")

  bf_orders <- table(rowSums(fit$vars))
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
