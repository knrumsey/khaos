#' Ordinal Probit Regression with BA-PCE
#'
#' The emulation approach of Francom et al. (2020) for BMARS, modified for polynomial chaos.
#'
#' @param X A data frame or matrix of predictors scaled to be between 0 and 1
#' @param y a response vector
#' @param degree Maximum polynomial degree for each basis function.
#' @param order Maximum order of interaction for each basis function.
#' @param nmcmc Number of MCMC iterations.
#' @param nburn Number of initial MCMC iterations to discard.
#' @param thin Keep every \code{thin} samples.
#' @param max_basis Maximum number of basis functions.
#' @param tau2 Prior variance for coefficients
#' @param sigma_thresh Prior variance for the threshold parameters.
#' @param h1,h2 Shape/scale parameters for the Gamma prior on expected number of basis functions.
#' @param move_probs A 3-vector with probabilities for (i) birth, (ii) death, and (iii) mutation.
#' @param coin_pars A list of control parameters for coinflip proposal
#' @param degree_penalty Increasing this value encourages lower order polynomial terms (0 is no penalization).
#' @param verbose Logical. Should progress information be printed?
#' @details Implements the RJMCMC algorithm described by Francom & Sanso (2020) for BMARS, modifying it for polynomial chaos basis functions.
#' As an alternative the NKD procedure of Nott et al. (2005), we use a coinflipping procedure to identify useful variables. See writeup (coming soon) for details.
#' The coin_pars argument is an ordered list with elements:
#' \enumerate{
#' \item{a function giving proposal probability for q_0, the expected order of interaction,}
#' \item{epsilon, the base weight for each variable (epsilon=Inf corresponds to uniform sampling),}
#' \item{alpha, the exponent-learning rate,}
#' \item{num_passes a numerical parameter only.}
#' }
#'
#' Other proposal notes:
#'    Expected order q0 is chosen from 1:order with weights `coin_pars[[1]](1:order)`
#'    Degree is chosen from q0:degree with weights `(1:(degree-q0+1))^(-degree_penalty)` and a random partition is created (see helper function).
#'
#' Two types of change steps are possible (and equally likely)
#'    A re-ordering of the degrees
#'    A single variable is swapped out with another one
#'
#' @references Francom, Devin, and Bruno Sansó. "BASS: An R package for fitting and performing sensitivity analysis of Bayesian adaptive spline surfaces." Journal of Statistical Software 94.LA-UR-20-23587 (2020).
#'
#' Nott, David J., Anthony YC Kuk, and Hiep Duc. "Efficient sampling schemes for Bayesian MARS models with many predictors." Statistics and Computing 15 (2005): 93-101.
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- adaptive_khaos(X, y)
#' @export
ordinal_khaos <-function(X, y,
                          degree=15, order=5,
                          nmcmc=10000,
                          nburn=9000,
                          thin=1,
                          max_basis=1000,
                          tau2=10^5,
                          sigma_thresh=10,
                          h1=4,h2=40/length(y),
                          move_probs=rep(1/3, 3),
                          coin_pars=list(function(j) 1/j, 1, 2, 3),
                          degree_penalty=0,
                          verbose=TRUE
){
  # Define a safer solve
  # Note: Under iid uniform X, columns of B are uncorrelated.
  #       So the tolerance should be pretty safe unless there's an
  #       extreme amount of correlation in X.
  safe_solve <- function(A, tol = 1e-9) {
    rc <- rcond(A)
    if (is.na(rc) || rc < tol) return(FALSE)
    solve(A)
  }
  is_odd <- function(n) (n %% 2) == 1

  n<-length(y)
  p<-ncol(X)
  if(p < order){
    order <- p
  }
  if(degree < order){
    stop("degree must be at least as large as order")
  }

  if(!all(unlist(lapply(y, function(yy) (yy %% 1) == 0)))){
    stop("y must be integer valued.")
  }
  K <- max(y)
  if(K < 2) stop("Need at least 2 categories")
  if(!all(y %in% 1:K)){
    stop("y must have values 1, 2, ... K")
  }
  thresh_curr <- mu_thresh <- 1:K - (K + is_odd(K)) / 2
  # Initialize z
  z <- thresh_curr[y] + rbeta(n, 5, 5) - 1
  ssz<-sum(z^2)
  thresh_curr[K] <- Inf

  s2z <- rep(NA, nmcmc)
  s2z[1] <- var(z)

  #KR: generalized harmonic number
  J_probs <- coin_pars[[1]](1:order)
  Hjs <- sum(J_probs)
  A_total <- A_size(p, degree, order)

  # storage
  degs<-vars<-array(dim=c(nmcmc,max_basis,order)) # this is the largest possible, our filling will be ragged

  nint<-dtot<-matrix(nrow=nmcmc,ncol=max_basis) # order of interaction J, and total degree, again filling will be ragged
  thresh_mat <- matrix(nrow=nmcmc, ncol=K-1)
  beta<-matrix(nrow=nmcmc,ncol=max_basis+1) # +1 for intercept, again filling will be ragged
  lam<-nbasis<- rep(NA,nmcmc) # error variance, poisson hyperprior, and number of basis functions
  #sum_sq <- rep(NA, nmcmc)
  eta_vec <- rep(0, p) # KR: count how many times each variable is added to the model

  # initialize
  nbasis[1]<-0
  s2 <- 1
  lam[1]<-1
  B.curr<-matrix(rep(1,n)) # matrix of current basis functions, so that yhat = B.curr %*% beta
  BtB.curr <- crossprod(B.curr)
  BtBi.curr <- solve(BtB.curr)

  Vinv.curr<-crossprod(B.curr)+1/tau2 # V^(-1) from DMS
  bhat<-1

  count_accept<-c(0,0,0,0)  # count how many times we accept birth, death, change degree, change vars
  count_propose<-c(0,0,0,0) # count how many times we propose each step
  names(count_accept) <- c("Birth", "Death", "Mutate (degree)", "Mutate (vars)")
  names(count_propose) <- c("Birth", "Death", "Mutate (degree)", "Mutate (vars)")
  beta[1,1]<-bhat

  if(verbose){
    pr<-c('MCMC iteration',0,myTimestamp(),'nbasis:', 0)
    cat(pr,'\n')
  }

  # start MCMC
  for(i in 2:nmcmc){
    ## Reversible jump step

    move.type<-sample(c('birth','death','change'),1,prob=move_probs)
    if(nbasis[i-1]<=1)
      move.type<-'birth'
    if(nbasis[i-1]==max_basis)
      move.type<-sample(c('death','change'),1,prob=move_probs[2:3])

    # set all of this iterations values to last iteration values...we'll change them if we accept a move below
    nbasis[i]<-nbasis[i-1]
    nint[i,]<-nint[i-1,]
    dtot[i,]<-dtot[i-1,]
    #knots[i,,]<-knots[i-1,,]
    #signs[i,,]<-signs[i-1,,]
    degs[i,,]<-degs[i-1,,]
    vars[i,,]<-vars[i-1,,]

    if(move.type=='birth'){
      # KR:
      q0  <- sample(1:order, 1, prob=J_probs)
      wts <- make_weights(eta_vec, q0, coin_pars[[3]], coin_pars[[2]], coin_pars[[4]])
      # Delayed rejection component
      chi.cand <- 0
      delayed_reject_term <- 0
      while(sum(chi.cand) == 0 | sum(chi.cand) > order){
        chi.cand  <- stats::rbinom(p, 1, wts)
        vars.cand <- which(chi.cand == 1)
        res <- 0
        for(jj in 1:order){
          wts.j <- make_weights(eta_vec, jj, coin_pars[[3]], coin_pars[[2]], coin_pars[[4]])
          term_j <- log(J_probs[jj]) - log(Hjs) + sum(log(wts.j[vars.cand])) + sum(log(1 - wts.j[-vars.cand]))
          res <- res + exp(term_j)
        }
        delayed_reject_term <- delayed_reject_term + log(res)
      }
      #vars.cand <- which(chi.cand == 1)
      nint.cand <- length(vars.cand)
      d_probs <- seq_along(nint.cand:degree)^(-degree_penalty)
      d_probs <- d_probs / sum(d_probs)
      dtot.cand <- sample(nint.cand:degree, 1, prob=d_probs)
      degs.cand <- random_partition(dtot.cand, nint.cand)

      if(nint.cand == 0) stop("Why is nint.cand zero? This shouldn't be possible.")
      basis.cand <- make_basis(vars.cand, degs.cand, X)
      B.cand <- cbind(B.curr,basis.cand) # add the new basis function to the basis functions we already have
      BtB.cand <- crossprod(B.cand)
      BtBi.cand <- tryCatch({
        safe_solve(BtB.cand)
      }, error = function(e) {
        FALSE
      })

      if(!isFALSE(BtBi.cand)){ # Otherwise reject immediately
        Vinv.cand<-crossprod(B.cand)+diag(nbasis[i-1]+2)/tau2 # +2: one for intercept and one for birth

        b.curr <- crossprod(B.curr, z)
        q.curr <- as.numeric(crossprod(b.curr, solve(Vinv.curr, b.curr)))
        logdet.curr <- as.numeric(determinant(Vinv.curr)$mod)

        b.cand <- crossprod(B.cand, z)
        q.cand <- as.numeric(crossprod(b.cand, solve(Vinv.cand, b.cand)))
        logdet.cand <- as.numeric(determinant(Vinv.cand)$mod)

        # calculate the log likelihood ratio (sort of) after integrating out beta and s2
        llik.alpha <- 0.5 * (q.cand - q.curr) - 0.5 * (logdet.cand - logdet.curr) + 0.5 * log(1/tau2)

        # log prior ratio
        lprior.alpha <- (
          log(lam[i-1]) - log(nbasis[i-1] + 1) # number of basis functions
          - log(A_total)                         # Prior over vars/degs
          + log(nbasis[i-1] + 1)                 # Account for ordering (a historical artifact?)
        )

        # log proposal ratio (proposed to current via death) - (current to proposed)
        lprop.alpha <- (
          (log(move_probs[2])                    # probability of selection a death step
           + log(1/(nbasis[i-1]+1)))              # probability that this basis function is selected to kill
          - (log(move_probs[1])                  # probability of selecting a birth step
             + log(d_probs[1+dtot.cand-nint.cand])  # KR: Probability of total degree
             - lchoose(dtot.cand, nint.cand)        # KR: Probability of degree partition
             + delayed_reject_term)                 # KR: probability of nint and vars based on coin flipping (delayed reject, see above)
        )

        alpha <- llik.alpha + lprior.alpha + lprop.alpha


        if(log(stats::runif(1))<alpha){ # accept, update current values
          B.curr     <-B.cand
          BtB.curr   <- BtB.cand
          BtBi.curr  <- BtBi.cand
          Vinv.curr<-Vinv.cand

          nbasis[i]<-nbasis[i-1]+1
          nint[i,nbasis[i]]<-nint.cand
          dtot[i,nbasis[i]]<-dtot.cand
          vars[i,nbasis[i],1:nint.cand]<-vars.cand
          degs[i,nbasis[i],1:nint.cand]<-degs.cand

          eta_vec[vars.cand] <- eta_vec[vars.cand] + 1
          count_accept[1] <- count_accept[1] + 1
        }
      }
      # Auto reject birth
      count_propose[1] <- count_propose[1] + 1

    } else if(move.type=='death'){
      tokill<-sample(nbasis[i-1],1) # which basis function we will delete
      B.cand<-B.curr[,-(tokill+1)]  # + 1 to skip the intercept
      BtB.cand <- crossprod(B.cand)
      # There really shouldn't be any issues here,
      # but for some reason we still need tryCatch
      BtBi.cand <- tryCatch({
        safe_solve(BtB.cand)
      }, error = function(e) {
        FALSE
      })

      if(!isFALSE(BtBi.cand)){ #Otherwise reject immediately
        Vinv.cand<-crossprod(B.cand)+diag(nbasis[i-1])/tau2

        b.curr <- crossprod(B.curr, z)
        q.curr <- as.numeric(crossprod(b.curr, solve(Vinv.curr, b.curr)))
        logdet.curr <- as.numeric(determinant(Vinv.curr)$mod)

        b.cand <- crossprod(B.cand, z)
        q.cand <- as.numeric(crossprod(b.cand, solve(Vinv.cand, b.cand)))
        logdet.cand <- as.numeric(determinant(Vinv.cand)$mod)

        # KR:
        vars.cand <- vars[i-1, tokill, 1:nint[i-1,tokill]]
        degs.cand <- degs[i-1, tokill, 1:nint[i-1,tokill]]
        dtot.cand <- sum(degs.cand)
        nint.cand <- length(vars.cand)
        d_probs <- seq_along(nint.cand:degree)^(-degree_penalty)
        d_probs <- d_probs / sum(d_probs)

        eta.cand <- eta_vec
        eta.cand[vars.cand] <- eta.cand[vars.cand] - 1
        proposal_nint   <- 0 # Probability of proposing the sum that we did
        proposal_nint_0 <- 0 # Probability of proposing a sum that equals 0
        proposal_nint_q <- 0 # Probability of proposing a sum that exceeds q=order
        for(jj in 1:order){
          wts.j <- make_weights(eta.cand, jj, coin_pars[[3]], coin_pars[[2]], coin_pars[[4]])
          term_j <- log(J_probs[jj]) - log(Hjs) + sum(log(wts.j[vars.cand])) + sum(log(1 - wts.j[-vars.cand]))
          proposal_nint <- proposal_nint + exp(term_j)

          # Include the delayed rejection stuff
          # P(sum(chi) = 0) = P(sum(chi) = 0 | q0=1)P(q0=1) + ... + P(sum(chi) = 0 | q0=q)P(q0=q)
          term_0_j <- log(J_probs[jj]) - log(Hjs) + sum(log(1 - wts.j))
          proposal_nint_0 <- proposal_nint_0 + exp(term_0_j)

          # if order == p, then it is not possible to overshoot
          term_q_j <- -Inf
          if(order < p){
            # Use a normal approx for this part
            mu_chi   <- sum(wts.j)
            sig_chi  <- sqrt(sum(wts.j*(1-wts.j)))
            term_q_j <- log(J_probs[jj]) - log(Hjs) +
              stats::pnorm(order+1/2,mu_chi,sig_chi,lower.tail=FALSE,log.p=TRUE)
          }
          proposal_nint_q <- proposal_nint_q + exp(term_q_j)
        }
        proposal_reject <- min(1, proposal_nint_0 + proposal_nint_q) # The 1 should never trigger, but I added it to be safe since we're using a normal approximation


        # calculate the log likelihood ratio (sort of) after integrating out beta and s2
        llik.alpha <- 0.5 * (q.cand - q.curr) - 0.5 * (logdet.cand - logdet.curr) - 0.5 * log(1/tau2)

        # log prior ratio
        lprior.alpha <- (
          - log(lam[i-1]) + log(nbasis[i-1])     # number of basis functions
          + log(A_total)                         # Prior over vars/degs
          - log(nbasis[i-1] + 1)                 # Account for ordering (a historical artifact?)
        )

        # log proposal ratio (proposed to current via birth) - (current to proposed)
        lprop.alpha <- (
          (log(move_probs[1])                               # probability of selection a birth step
           + log(proposal_nint) - log(1 - proposal_reject)   # probability that we proposed the variables we did (accounting for rejection shenanigans)
           + log(d_probs[1+dtot.cand-nint.cand])             # probability that we proposed the degree that we did
           - lchoose(dtot.cand, nint.cand))                  # probability that we proposed the degree partition that we did
          - (log(move_probs[2])                             # probability of selecting a death step
             + log(1/nbasis[i-1]))                             # probability of selecting this particular basis function to kill
        )

        alpha <- llik.alpha + lprior.alpha + lprop.alpha

        if(log(stats::runif(1))<alpha){ # accept, update
          B.curr<-B.cand
          BtB.curr <- BtB.cand
          BtBi.curr <- BtBi.cand
          Vinv.curr<-Vinv.cand

          nbasis[i]<-nbasis[i-1]-1
          nint[i,]<-dtot[i,]<-NA
          vars[i,,]<-degs[i,,]<-NA
          if(nbasis[i] > 0){
            nint[i,1:nbasis[i]]<-nint[i-1,(1:nbasis[i-1])[-tokill]]
            dtot[i,1:nbasis[i]]<-dtot[i-1,(1:nbasis[i-1])[-tokill]]

            vars[i,1:nbasis[i],]<-vars[i-1,(1:nbasis[i-1])[-tokill],]
            degs[i,1:nbasis[i],]<-degs[i-1,(1:nbasis[i-1])[-tokill],]
          }
          eta_vec <- eta.cand
          count_accept[2] <- count_accept[2] + 1
        }
      }
      count_propose[2] <- count_propose[2] + 1

    } else{ # KR: mutate

      # Adapt which mutation should get used.
      if(p <= 3){
        # Second mutation type isn't really needed for small p
        mutate_type <- 1
      }else{
        mutate_eps <- 0.1 # hard coded (worst case probability for either type)
        success_rates <- count_accept[3:4]/(0.01 + count_propose[3:4])
        mutate_prob <- mutate_eps + (1 - 2*mutate_eps) * stats::pnorm(-diff(pmin(5,pmax(-5,stats::qnorm(success_rates)))))
        mutate_type <- stats::rbinom(1, 1, mutate_prob)
      }

      # MUTATION TYPE 1: Re-sample the degrees
      if(mutate_type == 1){
        tochange<-sample(nbasis[i-1],1)                        # which basis function we will change

        # Get current info
        vars.curr <- vars[i-1,tochange,1:nint[i-1,tochange]]
        degs.curr <- degs[i-1,tochange,1:nint[i-1,tochange]]
        dtot.curr <- sum(degs.curr)    # equivalently dtot[i-1,tochange]
        nint.curr <- length(vars.curr) # equivalently nint[i-1,tochange]
        d_probs <- seq_along(nint.curr:degree)^(-degree_penalty)
        d_probs <- d_probs / sum(d_probs)

        # Generate new degree partition
        dtot.cand <- sample(nint.curr:degree, 1, prob=d_probs)
        degs.cand <- random_partition(dtot.cand, nint.curr)

        # Generate candidate basis function
        basis <- make_basis(vars.curr, degs.cand, X)
        B.cand <- B.curr
        B.cand[,tochange+1] <- basis # replace with our new basis function (+1 for intercept)
        BtB.cand <- crossprod(B.cand)
        BtBi.cand <- tryCatch({
          safe_solve(BtB.cand)
        }, error = function(e) {
          FALSE
        })

        if(!isFALSE(BtBi.cand)){ #Otherwise reject immediately
          Vinv.cand<-crossprod(B.cand)+diag(nbasis[i-1]+1)/tau2

          b.curr <- crossprod(B.curr, z)
          q.curr <- as.numeric(crossprod(b.curr, solve(Vinv.curr, b.curr)))
          logdet.curr <- as.numeric(determinant(Vinv.curr)$mod)

          b.cand <- crossprod(B.cand, z)
          q.cand <- as.numeric(crossprod(b.cand, solve(Vinv.cand, b.cand)))
          logdet.cand <- as.numeric(determinant(Vinv.cand)$mod)

          # Compute acceptance ratio
          # Assume that "interaction" between the two mutation types is negligible
          llik.alpha <- 0.5 * (q.cand - q.curr) - 0.5 * (logdet.cand - logdet.curr)

          # Variation in the degree proposal
          lprop.alpha <- log(d_probs[1+dtot.curr-nint.curr]) - log(d_probs[1+dtot.cand-nint.curr]) +
            (-lchoose(dtot.curr, nint.curr) + lchoose(dtot.cand, nint.curr))

          # Uniform priors + no dimension change <=> no term for the prior
          alpha <- llik.alpha + lprop.alpha

          if(log(stats::runif(1))<alpha){ # accept, update
            B.curr<-B.cand
            BtB.curr <- BtB.cand
            BtBi.curr <- BtBi.cand
            Vinv.curr<-Vinv.cand

            vars[i,tochange,1:nint[i-1,tochange]] <- vars.curr # no change
            degs[i,tochange,1:nint[i-1,tochange]] <- degs.cand
            nint[i,tochange] <- nint.curr                      # no change
            dtot[i,tochange] <- dtot.cand

            count_accept[3] <- count_accept[3] + 1
          }
        }
        count_propose[3] <- count_propose[3] + 1

      }else{
        # MUTATION TYPE 2: Swap out a variable
        tochange<-sample(nbasis[i-1],1)           # which basis function we will change

        # Get current info
        vars.curr <- vars[i-1,tochange,1:nint[i-1,tochange]]
        degs.curr <- degs[i-1,tochange,1:nint[i-1,tochange]]
        dtot.curr <- sum(degs.curr)    # equivalently dtot[i-1,tochange]
        nint.curr <- length(vars.curr) # equivalently nint[i-1,tochange]

        # Choose a variable to delete
        ind_tochange <- sample(nint.curr, 1)
        var_tochange <- vars.curr[ind_tochange]
        newvar_probs.curr <- 1e-9 + coin_pars[[2]] + eta_vec
        newvar_probs.curr[var_tochange] <- 0
        newvar_probs.curr <- newvar_probs.curr / sum(newvar_probs.curr)
        newvar.curr <- sample(p, 1, prob=newvar_probs.curr)

        vars.cand <- vars.curr
        vars.cand[ind_tochange] <- newvar.curr

        # Construct probabilities for the reverse case
        eta.cand <- eta_vec
        eta.cand[var_tochange] <- eta.cand[var_tochange] - 1
        eta.cand[newvar.curr] <- eta.cand[newvar.curr] + 1
        newvar_probs.cand <- 1e-9 + coin_pars[[2]] + eta.cand
        newvar_probs.cand[newvar.curr] <- 0
        newvar_probs.cand <- newvar_probs.cand / sum(newvar_probs.cand)
        newvar.cand <- var_tochange

        # Generate candidate basis function
        basis <- make_basis(vars.cand, degs.curr, X)
        B.cand <- B.curr
        B.cand[,tochange+1] <- basis # replace with our new basis function (+1 for intercept)
        BtB.cand <- crossprod(B.cand)
        BtBi.cand <- tryCatch({
          safe_solve(BtB.cand)
        }, error = function(e) {
          FALSE
        })

        if(!isFALSE(BtBi.cand)){ # Otherwise reject immediately
          Vinv.cand<-crossprod(B.cand)+diag(nbasis[i-1]+1)/tau2

          b.curr      <- crossprod(B.curr, z)
          q.curr      <- as.numeric(crossprod(b.curr, solve(Vinv.curr, b.curr)))
          logdet.curr <- as.numeric(determinant(Vinv.curr)$mod)

          b.cand      <- crossprod(B.cand, z)
          q.cand      <- as.numeric(crossprod(b.cand, solve(Vinv.cand, b.cand)))
          logdet.cand <- as.numeric(determinant(Vinv.cand)$mod)

          # Compute acceptance ratio
          # Assume that "interaction" between the two mutation types is negligible
          llik.alpha <- 0.5 * (q.cand - q.curr) - 0.5 * (logdet.cand - logdet.curr)

          # Difference in probability of the variables.
          lprop.alpha <- log(newvar_probs.cand[newvar.cand]) - log(newvar_probs.curr[newvar.curr])

          # Uniform priors + no dimension change <=> no term for the prior
          alpha <- llik.alpha + lprop.alpha

          if(log(stats::runif(1))<alpha){ # accept, update
            B.curr<-B.cand
            BtB.curr <- BtB.cand
            BtBi.curr <- BtBi.cand
            Vinv.curr<-Vinv.cand

            vars[i,tochange,1:nint[i,tochange]] <- vars.cand
            degs[i,tochange,1:nint[i,tochange]] <- degs.curr # no change
            nint[i,tochange] <- nint.curr                    # no change
            dtot[i,tochange] <- dtot.cand                    # no change

            eta_vec <- eta.cand

            count_accept[4] <- count_accept[4] + 1
          }
        }
        count_propose[4] <- count_propose[4] + 1
      }
    }

    ## Gibbs steps

    # SAMPLE COEFFICIENTS
    Lambda_n <- BtB.curr + diag(nbasis[i] + 1)/tau2
    Lambda_i_n <- solve(Lambda_n)
    beta_ls <- BtBi.curr %*% crossprod(B.curr, z)
    mu_n <- Lambda_i_n %*% (BtB.curr %*% beta_ls)
    Sigma_n <- s2 * Lambda_i_n
    beta_curr <- t(mvtnorm::rmvnorm(1, mu_n, Sigma_n))
    beta[i,1:(nbasis[i]+1)] <- beta_curr
    zhat_curr <- B.curr %*% beta_curr
    resid <- z-zhat_curr

    # SAMPLE LAMBDA
    lam[i] <- stats::rgamma(1,h1+nbasis[i],h2+1)                                               # update lambda

    # SAMPLE LATENT Zs
    for(j in 1:n){
      ez <- zhat_curr[j]
      a <- max(-Inf, thresh_curr[y[j]-1], na.rm=TRUE)
      b <- thresh_curr[y[j]]  # min(thresh_curr[y[i]], Inf, na.rm=TRUE)
      u <- runif(1, pnorm(a-ez), pnorm(b-ez))
      z[j] <- ez + qnorm(u)
    }
    ssz <- sum(z^2)
    s2z[i] <- (ssz - sum(z)^2/n) / (n-1)

    # SAMPLE THRESHOLDS
    for(k in 2:(K-1)){
      a <- max(z[which(y == k)])
      b <- min(z[which(y == k+1)])
      u <- runif(1,
                 pnorm(a, mu_thresh[k], sigma_thresh),
                 pnorm(b, mu_thresh[k], sigma_thresh))
      thresh_curr[k] <- mu_thresh[k] + sigma_thresh * qnorm(u)
    }
    thresh_mat[i,] <- thresh_curr[1:(K-1)]

    if(verbose & i%%1000 == 0){
      pr<-c('MCMC iteration',i,myTimestamp(),'nbasis:', nbasis[i])
      cat(pr,'\n')
    }
  }

  # Trim down data structures
  mcmc_iter <- seq(nburn+1, nmcmc, by=thin)
  basis_high <- 1:max(nbasis)
  inter_high <- 1:max(nint, na.rm=TRUE)

  vars <- vars[mcmc_iter, basis_high, inter_high, drop=FALSE]
  degs <- degs[mcmc_iter, basis_high, inter_high, drop=FALSE]
  nint <- nint[mcmc_iter, basis_high, drop=FALSE]
  dtot <- dtot[mcmc_iter, basis_high, drop=FALSE]

  beta <- beta[mcmc_iter, 1:(1+max(nbasis))]
  nbasis <- nbasis[mcmc_iter]
  s2z <- s2z[mcmc_iter]
  #s2 <- s2[mcmc_iter]
  lam <- lam[mcmc_iter]
  thresh_mat <- thresh_mat[mcmc_iter,,drop=FALSE]

  out <- list(B=B.curr,
              vars=vars,degs=degs,
              nint=nint,dtot=dtot,
              thresh=thresh_mat,
              nbasis=nbasis,beta=beta,
              lam=lam,
              eta=eta_vec,
              s2z=s2z,
              count_accept=count_accept,
              count_propose=count_propose,
              X=X, y=y)
  class(out) <- "ordinal_khaos"
  return(out)
}

#' Predict Method for class \code{ordinal_khaos}
#'
#' Compute posterior predictive quantities for an ordinal probit BA-PCE model.
#' For each selected MCMC draw, the method builds the polynomial–chaos basis,
#' forms the linear predictor \eqn{f(x) = B(x)\beta}, reconstructs the ordered
#' cutpoints, and returns either category probabilities, MAP classes, or latent
#' means. Results can be averaged across draws or returned per-draw.
#'
#' @param object An object of class \code{ordinal_khaos} returned by \code{ordinal_khaos()}.
#' @param newdata A data frame or matrix of predictors on the same scale and with
#'   the same columns as used for fitting. If \code{NULL}, uses \code{object$X}.
#' @param mcmc.use Integer vector of posterior draw indices to use. If \code{NULL},
#'   uses all stored draws (i.e., \code{seq_along(object$nbasis)}).
#' @param type Character; one of \code{"prob"}, \code{"class"}, or \code{"latent"}.
#'   \code{"prob"} returns category probabilities \eqn{P(Y=k \mid x)} for
#'   \eqn{k=1,\dots,K}. \code{"class"} returns MAP classes. \code{"latent"} returns
#'   the latent mean \eqn{f(x)} (i.e., \eqn{B(x)\beta}) for each draw.
#' @param aggregate Logical; if \code{TRUE} (default), averages results across the
#'   selected MCMC draws (probabilities are averaged; classes are aggregated by
#'   simple voting; latent means are averaged). If \code{FALSE}, returns per-draw
#'   arrays.
#' @param ... Unused; included for method compatibility.
#'
#' @details
#' Let \eqn{Z = f_\beta(X) + \varepsilon} with \eqn{\varepsilon \sim \mathcal{N}(0,1)} and
#' ordered cutpoints \eqn{\gamma_0=-\infty < \gamma_1 < \cdots < \gamma_{K-1} < \gamma_K=\infty}.
#' For a new input \eqn{x}, the category probabilities are
#' \deqn{P(Y=k\mid x) = \Phi(\gamma_k - f(x)) - \Phi(\gamma_{k-1} - f(x)).}
#' The function reconstructs the full cutpoint vector internally. If the fitted
#' object stored \eqn{\gamma_1,\dots,\gamma_{K-1}} then these are used directly;
#' if it stored only \eqn{\gamma_2,\dots,\gamma_{K-1}} (with \eqn{\gamma_1 \equiv 0})
#' the missing fixed cutpoint is inserted automatically.
#'
#' @return
#' If \code{type = "prob"} and \code{aggregate = TRUE}, an \code{nrow(newdata) × K}
#' numeric matrix of posterior mean category probabilities. If \code{aggregate = FALSE},
#' a 3-D array of dimension \code{(#draws) × nrow(newdata) × K}.
#'
#' If \code{type = "class"} and \code{aggregate = TRUE}, an integer vector of length
#' \code{nrow(newdata)} with MAP classes aggregated by voting across draws. If
#' \code{aggregate = FALSE}, a \code{(#draws) × nrow(newdata)} integer matrix of
#' per-draw MAP classes.
#'
#' If \code{type = "latent"}, returns the latent mean(s) \eqn{f(x)}: a numeric vector
#' of length \code{nrow(newdata)} if \code{aggregate = TRUE}, or a
#' \code{(#draws) × nrow(newdata)} numeric matrix otherwise.
#'
#' @references
#' Hoff, Peter D. (2009). *A First Course in Bayesian Statistical Methods*.
#' Springer. (See ordinal probit data augmentation.)
#'
#' Francom, Devin, and Bruno Sansó (2020). "BASS: An R package for fitting and
#' performing sensitivity analysis of Bayesian adaptive spline surfaces."
#' \emph{Journal of Statistical Software} 94. LA-UR-20-23587.
#'
#' @seealso \code{\link{ordinal_khaos}} for model fitting.
#'
#' @examples
#' \dontrun{
#' # Toy example (X scaled to [0,1])
#' set.seed(1)
#' X <- lhs::maximinLHS(200, 3)
#' # Suppose y has K=4 ordered categories:
#' fit <- ordinal_khaos(X, sample(1:4, 200, replace=TRUE),
#'                      nmcmc = 2000, nburn = 1000, thin = 2, verbose = FALSE)
#'
#' # Posterior mean probabilities for each class
#' p_hat <- predict(fit, newdata = X, type = "prob", aggregate = TRUE)
#'
#' # MAP classes (aggregated across draws)
#' y_hat <- predict(fit, newdata = X, type = "class")
#'
#' # Latent mean f(x)
#' f_hat <- predict(fit, newdata = X, type = "latent")
#' }
#'
#' @export
predict.ordinal_khaos <- function(object, newdata = NULL,
                                  mcmc.use = NULL,
                                  type = c("class","prob","latent"),
                                  aggregate = FALSE, ...) {
  type <- match.arg(type)
  print(type)
  if (is.null(newdata)) newdata <- object$X
  K <- max(object$y)
  n <- nrow(newdata)

  # which posterior draws to use
  if (is.null(mcmc.use)) mcmc.use <- seq_along(object$nbasis)
  I <- length(mcmc.use)

  # helper: build B for a given draw i
  build_B <- function(i) {
    nb <- object$nbasis[i]
    B <- matrix(1, nrow = n, ncol = nb + 1)  # intercept
    if (nb > 0) {
      for (j in seq_len(nb)) {
        J <- object$nint[i, j]
        if (is.na(J) || J == 0) next
        vars_ij <- object$vars[i, j, 1:J]
        degs_ij <- object$degs[i, j, 1:J]
        B[, j + 1] <- make_basis(vars_ij, degs_ij, newdata)
      }
    }
    B
  }

  # helper: extend thresholds to (-Inf, ..., +Inf)
  make_gamma_ext <- function(i) {
    thr <- object$thresh[i, ]
    thr <- as.numeric(thr)
    # stored length could be K-1 (gamma1..gamma_{K-1}) or K-2 (gamma2..gamma_{K-1} with gamma1 fixed at 0)
    if (length(thr) == (K - 1)) {
      # assume thr = (gamma1, ..., gamma_{K-1})
      c(-Inf, thr, Inf)
    } else if (length(thr) == (K - 2)) {
      # assume gamma1 == 0 is fixed
      c(-Inf, 0, thr, Inf)
    } else {
      stop("object$thresh has incompatible number of columns.")
    }
  }

  # allocate outputs
  if (type == "prob") {
    if (aggregate) {
      out <- matrix(0, nrow = n, ncol = K)
    } else {
      out <- array(NA_real_, dim = c(I, n, K))
    }
  } else if (type == "class") {
    out <- integer(n)  # aggregated MAP class
    if (!aggregate) {
      out <- matrix(NA_integer_, nrow = I, ncol = n)
    }
  } else { # latent
    if (aggregate) {
      out <- numeric(n)
    } else {
      out <- matrix(NA_real_, nrow = I, ncol = n)
    }
  }

  # loop over draws
  for (idx in seq_along(mcmc.use)) {
    i <- mcmc.use[idx]
    B <- build_B(i)
    beta_i <- object$beta[i, 1:ncol(B)]
    f <- as.numeric(B %*% beta_i)

    if (type == "latent") {
      if (aggregate) {
        out <- out + f / I
      } else {
        out[idx, ] <- f
      }
      next
    }

    # probs from cutpoints
    gamma_ext <- make_gamma_ext(i)  # length K+1
    # P(y=k|x) = Phi(gamma_k - f) - Phi(gamma_{k-1} - f)
    # vectorize over k by recycling f
    pk <- matrix(NA_real_, nrow = n, ncol = K)
    for (k in 1:K) {
      pk[, k] <- stats::pnorm(gamma_ext[k + 1] - f) - stats::pnorm(gamma_ext[k] - f)
    }

    if (type == "prob") {
      if (aggregate) {
        out <- out + pk / I
      } else {
        out[idx, , ] <- pk
      }
    } else { # class
      cls <- max.col(pk, ties.method = "first")
      if (aggregate) {
        # running vote: accumulate probabilities then decide at end is better,
        # but for simplicity we vote by MAP per draw
        if (idx == 1) vote <- matrix(0, nrow = n, ncol = K)
        vote[cbind(seq_len(n), cls)] <- vote[cbind(seq_len(n), cls)] + 1
        if (idx == I) out <- max.col(vote, ties.method = "first")
      } else {
        out[idx, ] <- cls
      }
    }
  }

  return(out)
}


#' Plot Method for class \code{ordinal_khaos}
#'
#' Produces a 2×2 diagnostic plot: (1) trace of the number of basis functions;
#' (2) jittered predicted vs. actual classes with point size indicating posterior
#' support; (3) posterior mean class probabilities on \code{x$X} (lines if 1D,
#' otherwise barplot of averages); (4) histogram of latent means \eqn{f(x)} with
#' posterior-average thresholds overlaid.
#'
#' @param x An object of class \code{ordinal_khaos}.
#'
#' @details
#' Panel (2) uses \code{predict(x, type="class", aggregate=FALSE)} to compute
#' per-draw MAP classes, plotting the modal class vs. observed class; point size
#' reflects the fraction of draws supporting the modal class. Panel (3) uses
#' \code{predict(x, type="prob")} to plot posterior mean class probabilities.
#' Panel (4) uses \code{predict(x, type="latent", aggregate=FALSE)} to pool
#' latent means across draws and overlays posterior-average cutpoints (handles
#' either stored \eqn{\gamma_1,\dots,\gamma_{K-1}} or \eqn{\gamma_2,\dots,\gamma_{K-1}}
#' with \eqn{\gamma_1\equiv 0}).
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(runif(200), ncol = 1)
#' y <- sample(1:4, 200, replace = TRUE)
#' fit <- ordinal_khaos(X, y, nmcmc = 2000, nburn = 1000, thin = 2, verbose = FALSE)
#' plot(fit)
#' }
#' @export
plot.ordinal_khaos <- function(x, ...){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(2, 2))

  K <- max(x$y)

  ## --- (1) Trace of nbasis ---
  plot(seq_along(x$nbasis), x$nbasis, type = "l", lwd = 1,
       xlab = "MCMC iteration", ylab = "Number of basis functions",
       main = "Trace: nbasis")

  ## --- (2) Predicted class (uncertainty) vs actual class ---
  cls_draws <- predict(x, newdata = x$X, type = "class", aggregate = FALSE)
  if (is.vector(cls_draws)) cls_draws <- matrix(cls_draws, nrow = 1)
  n <- ncol(cls_draws)

  modal_class <- integer(n)
  support <- numeric(n)
  for (j in seq_len(n)) {
    tabj <- tabulate(cls_draws[, j], nbins = K)
    modal_class[j] <- which.max(tabj)
    support[j] <- max(tabj) / nrow(cls_draws)
  }
  # map support to point size
  cex_j <- 0.2 + 1.6 * (support - min(support)) /
    max(1e-12, (max(support) - min(support)))

  plot(jitter(x$y, 0.15), jitter(modal_class, 0.15),
       xlab = "Actual class", ylab = "Predicted class (modal)",
       xlim = c(0.5, K + 0.5), ylim = c(0.5, K + 0.5),
       axes = FALSE, pch = 19, col = "black", cex = cex_j,
       main = "Predicted vs actual (size = support)")
  abline(0,1,col='red')
  axis(1, at = 1:K, labels = 1:K); axis(2, at = 1:K, labels = 1:K, las = 1); box()

  ## --- (3) Posterior mean class probabilities on training X ---
  probs <- predict(x, newdata = x$X, type = "prob", aggregate = TRUE)
  if (ncol(x$X) == 1) {
    ord <- order(x$X[, 1])
    matplot(x$X[ord, 1], probs[ord, ], type = "l", lty = 1, lwd = 2,
            xlab = "X", ylab = "Predicted probability",
            main = "Posterior mean class probabilities", col = seq_len(K))
    legend("topright", legend = paste("Class", 1:K),
           col = seq_len(K), lty = 1, lwd = 2, bty = "n")
  } else {
    barplot(colMeans(probs), names.arg = paste("Class", 1:K),
            ylab = "Mean predicted probability",
            main = "Avg class probabilities (train)", col = seq_len(K))
  }

  ## --- (4) Histogram of latent means f(x) with avg thresholds ---
  f_draws <- predict(x, newdata = x$X, type = "latent", aggregate = FALSE)  # I x n
  f_vec <- as.numeric(f_draws)  # pool across draws & cases
  hist(f_vec, breaks = "FD", freq = TRUE, col = "grey85", border = "black",
       xlab = "Latent mean f(x)", main = "Latent means with avg thresholds")

  # compute posterior-average thresholds
  thresh <- x$thresh
  thr_len <- ncol(thresh)
  if (thr_len == (K - 1)) {
    # assume thresh = (gamma1..gamma_{K-1})
    gamma_mean <- c(-Inf, colMeans(thresh), Inf)
  } else if (thr_len == (K - 2)) {
    # assume gamma1 == 0 is fixed
    gamma_mean <- c(-Inf, 0, colMeans(thresh), Inf)
  } else {
    gamma_mean <- c(-Inf, 0, Inf)  # fallback to avoid errors
  }

  # draw interior thresholds only (exclude +/-Inf)
  if (length(gamma_mean) >= 3) {
    g_int <- gamma_mean[is.finite(gamma_mean)]
    abline(v = g_int, col = "red", lwd = 2, lty = 2)
    legend("topright", legend = "Avg thresholds", lty = 2, lwd = 2, col = "red", bty = "n")
  }

  invisible(NULL)
}
