#' Bayesian Adaptive Polynomial Chaos Expansion
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
#' @param s2_lower Lower bound on process variance (numerically useful for deterministic functions).
#' @param g1,g2 Shape/scale parameters for the IG prior on process variance (default is Jefffrey's prior)
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
adaptive_khaos <-function(X, y,
                          degree=15, order=5,
                          nmcmc=10000,
                          nburn=9000,
                          thin=1,
                          max_basis=1000,
                          tau2=10^5,
                          s2_lower=0,
                          g1=0,g2=0,
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

  n<-length(y)
  p<-ncol(X)
  ssy<-sum(y^2)
  if(p < order){
    order <- p
  }
  if(degree < order){
    stop("degree must be at least as large as order")
  }
  #KR: generalized harmonic number
  J_probs <- coin_pars[[1]](1:order)
  Hjs <- sum(J_probs)
  A_total <- A_size(p, degree, order)

  # storage
  degs<-vars<-array(dim=c(nmcmc,max_basis,order)) # this is the largest possible, our filling will be ragged

  nint<-dtot<-matrix(nrow=nmcmc,ncol=max_basis) # order of interaction J, and total degree, again filling will be ragged
  beta<-matrix(nrow=nmcmc,ncol=max_basis+1) # +1 for intercept, again filling will be ragged
  s2<-lam<-nbasis<-rep(NA,nmcmc) # error variance, poisson hyperprior, and number of basis functions
  sum_sq <- rep(NA, nmcmc)
  eta_vec <- rep(0, p) # KR: count how many times each variable is added to the model

  # initialize
  nbasis[1]<-0
  s2[1]<-stats::var(y)
  lam[1]<-1
  B.curr<-matrix(rep(1,n)) # matrix of current basis functions, so that yhat = B.curr %*% beta
  BtB.curr <- crossprod(B.curr)
  BtBi.curr <- solve(BtB.curr)

  Vinv.curr<-crossprod(B.curr)+1/tau2 # V^(-1) from DMS
  bhat<-solve(Vinv.curr)%*%t(B.curr)%*%y
  d.curr <- g2 + ssy - t(bhat)%*%Vinv.curr%*%bhat # d from DMS

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
      #nint.cand<-sample(order,1) # sample degree of interaction for new basis function (different from DMS)
      #vars.cand<-sample(p,nint.cand,replace = F) # variables to use in new basis function (different from DMS)
      #knots.cand<-stats::runif(nint.cand) # sample knots for new basis function
      #signs.cand<-sample(c(-1,1),nint.cand,replace = T) # signs for new basis function

      # KR:


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
        bhat.cand<-solve(Vinv.cand)%*%t(B.cand)%*%y
        d.cand <- g2 + ssy - t(bhat.cand)%*%Vinv.cand%*%bhat.cand

        # calculate the log likelihood ratio (sort of) after integrating out beta and s2
        llik.alpha <- (
                      0.5*log(1/tau2) +
                      (determinant(Vinv.curr)$mod/2 - determinant(Vinv.cand)$mod/2) +
                      (g1+n/2)*(log(d.curr) - log(d.cand))
        )

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
          B.curr<-B.cand
          BtB.curr <- BtB.cand
          BtBi.curr <- BtBi.cand
          bhat<-bhat.cand
          Vinv.curr<-Vinv.cand
          d.curr<-d.cand
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
        bhat.cand<-solve(Vinv.cand)%*%crossprod(B.cand,y)
        d.cand <- g2 + ssy - crossprod(bhat.cand,Vinv.cand%*%bhat.cand)

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
        llik.alpha <- (
          -0.5*log(1/tau2) +
            (determinant(Vinv.curr)$mod/2 - determinant(Vinv.cand)$mod/2) +
            (g1+n/2)*(log(d.curr) - log(d.cand))
        )

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
          bhat<-bhat.cand
          Vinv.curr<-Vinv.cand
          d.curr<-d.cand
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
          bhat.cand<-solve(Vinv.cand)%*%crossprod(B.cand,y)
          d.cand <- g2 + ssy - crossprod(bhat.cand,Vinv.cand%*%bhat.cand)

          # Compute acceptance ratio
          # Assume that "interaction" between the two mutation types is negligible
          llik.alpha <- determinant(Vinv.curr)$mod/2 - determinant(Vinv.cand)$mod/2 +
                        (g1+n/2)*(log(d.curr) - log(d.cand))

          # Variation in the degree proposal
          lprop.alpha <- log(d_probs[1+dtot.curr-nint.curr]) - log(d_probs[1+dtot.cand-nint.curr]) +
                         (-lchoose(dtot.curr, nint.curr) + lchoose(dtot.cand, nint.curr))

          # Uniform priors + no dimension change <=> no term for the prior
          alpha <- llik.alpha + lprop.alpha

          if(log(stats::runif(1))<alpha){ # accept, update
            B.curr<-B.cand
            BtB.curr <- BtB.cand
            BtBi.curr <- BtBi.cand
            bhat<-bhat.cand
            Vinv.curr<-Vinv.cand
            d.curr<-d.cand

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
          bhat.cand<-solve(Vinv.cand)%*%crossprod(B.cand,y)
          d.cand <- g2 + ssy - crossprod(bhat.cand,Vinv.cand%*%bhat.cand)

          # Compute acceptance ratio
          # Assume that "interaction" between the two mutation types is negligible
          llik.alpha <- determinant(Vinv.curr)$mod/2 - determinant(Vinv.cand)$mod/2 +
            (g1+n/2)*(log(d.curr) - log(d.cand))

          # Difference in probability of the variables.
          lprop.alpha <- log(newvar_probs.cand[newvar.cand]) - log(newvar_probs.curr[newvar.curr])

          # Uniform priors + no dimension change <=> no term for the prior
          alpha <- llik.alpha + lprop.alpha

          if(log(stats::runif(1))<alpha){ # accept, update
            B.curr<-B.cand
            BtB.curr <- BtB.cand
            BtBi.curr <- BtBi.cand
            bhat<-bhat.cand
            Vinv.curr<-Vinv.cand
            d.curr<-d.cand

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
    Lambda_n <- BtB.curr + diag(nbasis[i] + 1)/tau2
    Lambda_i_n <- solve(Lambda_n)
    beta_ls <- BtBi.curr %*% crossprod(B.curr, y)
    mu_n <- Lambda_i_n %*% (BtB.curr %*% beta_ls)
    Sigma_n <- s2[i-1] * Lambda_i_n
    beta_curr <- t(mvtnorm::rmvnorm(1, mu_n, Sigma_n))
    beta[i,1:(nbasis[i]+1)] <- beta_curr
    yhat_curr <- B.curr %*% beta_curr
    resid <- y-yhat_curr


    s2[i] <- max(s2_lower, 1/stats::rgamma(1, n/2+g1, rate = g2+.5*sum(resid^2)))
    lam[i] <- stats::rgamma(1,h1+nbasis[i],h2+1)                                               # update lambda
    sum_sq[i] <- sum(resid^2)

    # S_cov <- solve(crossprod(B.curr)/s2[i-1]+diag(nbasis[i]+1)/tau2)                    # covariance matrix for beta update
    # beta[i,1:(nbasis[i]+1)] <- sample_mvn(1, S_cov%*%t(B.curr)%*%y/s2[i-1], S_cov)      # update beta
    # s2[i] <- 1/stats::rgamma(1,n/2+g1,rate=g2+.5*sum((y-B.curr%*%beta[i,1:(nbasis[i]+1)])^2))  # update s2
    # #if(s2[i] > 1000) browser()

    if(verbose & i%%1000 == 0){
      pr<-c('MCMC iteration',i,myTimestamp(),'nbasis:', nbasis[i])
      cat(pr,'\n')
    }
  }

  # Trim down data structures
  mcmc_iter <- seq(nburn+1, nmcmc, by=thin)
  basis_high <- 1:max(nbasis)
  inter_high <- 1:max(nint, na.rm=TRUE)

  vars <- vars[mcmc_iter, basis_high, inter_high]
  degs <- degs[mcmc_iter, basis_high, inter_high]
  nint <- nint[mcmc_iter, basis_high]
  dtot <- dtot[mcmc_iter, basis_high]

  beta <- beta[mcmc_iter, 1:(1+max(nbasis))]
  nbasis <- nbasis[mcmc_iter]
  s2 <- s2[mcmc_iter]
  lam <- lam[mcmc_iter]


  out <- list(B=B.curr,
              vars=vars,degs=degs,
              nint=nint,dtot=dtot,
              nbasis=nbasis,beta=beta,
              s2=s2,lam=lam,sum_sq=sum_sq,
              eta=eta_vec,
              count_accept=count_accept,
              count_propose=count_propose,
              X=X, y=y)
  class(out) <- "adaptive_khaos"
  return(out)
}

#' Predict Method for class adaptive_khaos
#'
#' See \code{adaptive_khaos()} for details.
#'
#' @param object An object returned by the \code{adaptive_khaos()} function.
#' @param newdata A dataframe of the same dimension as the training data.
#' @param mcmc.use Which posterior samples should be used for prediction?
#' @param nugget Logical. Should predictions include error?
#' @param nreps How many predictions should be taken for each mcmc sample (ignored when \code{nugget = FALSE}).
#' @param ... Additional arguments to predict
#' @details Predict function for adaptive_khaos object.
#' @references Francom, Devin, and Bruno Sansó. "BASS: An R package for fitting and performing sensitivity analysis of Bayesian adaptive spline surfaces." Journal of Statistical Software 94.LA-UR-20-23587 (2020).
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- adaptive_khaos(X, y)
#' predict(fit)
#'
#' @export
predict.adaptive_khaos<-function(object, newdata=NULL, mcmc.use=NULL, nugget=FALSE, nreps=1, ...){ # prediction function, gets a prediction for each MCMC iteration
  if(is.null(newdata)){
    newdata <- object$X
  }
  if(is.null(mcmc.use)){
    mcmc.use <- seq_along(object$nbasis)
  }
  if(!nugget){
    nreps <- 1
  }
  pred <- matrix(NA, nrow=length(mcmc.use)*nreps, ncol=nrow(newdata))
  cnt <- 1
  for(i in mcmc.use){
    B <- matrix(1, nrow=nrow(newdata),ncol=object$nbasis[i]+1)
    for(j in 1:object$nbasis[i]){
      B[,j+1] <- make_basis(object$vars[i,j,1:object$nint[i,j]], object$degs[i,j,1:object$nint[i,j]], newdata)
    }
    beta_curr <- object$beta[i,1:(object$nbasis[i]+1)]
    if(nugget){
      mu <- matrix(rep(B %*% beta_curr, each=nreps),
                   nrow=nreps, ncol=nrow(newdata))
      sigma <- sqrt(object$s2[i])
      pred[(1 + (cnt-1)*nreps):(cnt*nreps),] <- mu + stats::rnorm(nrow(newdata)*nreps, 0, sigma)
    }else{
      pred[cnt,] <- B %*% beta_curr
    }
    cnt <- cnt + 1
  }
  return(pred)
}


#' Plot Method for class adaptive_khaos
#'
#' See \code{adaptive_khaos()} for details.
#'
#' @param x An object returned by the \code{adaptive_khaos()} function.
#' @param ... Additional arguments for plotting
#' @details Plot function for adaptive_khaos object.
#' @references Francom, Devin, and Bruno Sansó. "BASS: An R package for fitting and performing sensitivity analysis of Bayesian adaptive spline surfaces." Journal of Statistical Software 94.LA-UR-20-23587 (2020).
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- adaptive_khaos(X, y)
#' plot(fit)
#'
#' @export
plot.adaptive_khaos <- function(x, ...){
  opar <- graphics::par(no.readonly = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0))

  preds <- stats::predict(x, nugget=TRUE, nreps=10)
  yhat <- apply(preds, 2, mean)
  ci <- 2*apply(preds, 2, stats::sd)

  stats::ts.plot(x$nbasis, ylab="nbasis")
  stats::ts.plot(x$s2, ylab="s2")
  plot(x$y, yhat, pch=16, xlab="y")
  graphics::segments(x0=x$y, y0=yhat-ci, y1=yhat+ci, col='orange')

  graphics::points(x$y, yhat, pch=16)
  graphics::abline(0,1,col='dodgerblue')
  rr <- x$y-yhat
  graphics::hist(rr, breaks=ceiling(length(rr)^0.33*diff(range(rr))/(3.5*stats::sd(rr))), freq=F)
  x = seq(range(rr)[1], range(rr)[2], length.out=100)
  graphics::curve(stats::dnorm(x, mean(rr), stats::sd(rr)), add=TRUE, col='orange')

  graphics::par(opar)
  return(NULL)
}

