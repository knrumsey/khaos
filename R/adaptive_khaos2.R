#' Experimental g-prior: Bayesian Adaptive Polynomial Chaos Expansion
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
#' @param a_g,b_g Prior parameters for global g-prior term (on precision scale)
#' @param zeta Hyperparameter for modified g-prior. Setting zeta = 0 reduces to the usual g-prior.
#' @param g2_sample Character string specifying the method used to sample or update the prior scaling parameter \code{g0^2}.
#' @param g2_init Initial calue of g2 (global precision for g-prior). Becomes the only value when \code{g2_sample = "fixed"}.
#' @param s2_lower Lower bound on process variance (numerically useful for deterministic functions).
#' @param a_sigma,b_sigma Shape/scale parameters for the IG prior on process variance (default is Jefffrey's prior)
#' @param a_M,b_M Shape/scale parameters for the Gamma prior on expected number of basis functions.
#' @param move_probs A 3-vector with probabilities for (i) birth, (ii) death, and (iii) mutation.
#' @param coin_pars A list of control parameters for coinflip proposal
#' @param degree_penalty Increasing this value encourages lower order polynomial terms (0 is no penalization).
#' @param legacy Logical. If TRUE, mimics original implementation behavior (from the emulator comparison paper),
#'        which may retain invalid basis function structures across iterations. Defaults to FALSE.
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
#' Sampling method for g0^2 depends on \code{g2_sample}
#' \describe{
#'   \item{\code{"f"}}{Fixed: Do not update \code{g0^2}; carry forward the previous value (\code{g2[i] <- g2[i-1]}).}
#'   \item{\code{"lf"}}{Laplace-Full: Use a Laplace approximation to sample \code{g0^2} under the full marginal likelihood (no orthogonality assumption).}
#'   \item{\code{"lo"}}{Laplace-Orthogonal: Use a Laplace approximation to sample \code{g0^2} assuming the orthogonality approximation applies.}
#'   \item{\code{"mh"}}{Metropolis-Hastings with full Laplace approximation as proposal. Target is the full posterior density of \code{g0^2}.}
#'   \item{\code{"mho"}}{Metropolis-Hastings with orthogonal Laplace approximation as proposal. Target is the full posterior density of \code{g0^2}. under the orthogonality assumption.}
#'   \item{\code{"mhoo"}}{Metropolis-Hastings with orthogonal Laplace approximation as proposal. Target is the posterior density of \code{g0^2} under the orthogonality assumption.}
#' }
#'
#' Other proposal notes:
#'    Expected order q0 is chosen from 1:order with weights `coin_pars[[1]](1:order)`
#'    Degree is chosen from q0:degree with weights `(1:(degree-q0+1))^(-degree_penalty)` and a random partition is created (see helper function).
#'
#' Two types of change steps are possible (and equally likely)
#'    A re-ordering of the degrees
#'    A single variable is swapped out with another one (only used when \code{ncol(X) > 3})
#'
#' @references
#' Francom, Devin, and Bruno Sansó. "BASS: An R package for fitting and performing sensitivity analysis of Bayesian adaptive spline surfaces." Journal of Statistical Software 94.LA-UR-20-23587 (2020).
#'
#' Nott, David J., Anthony YC Kuk, and Hiep Duc. "Efficient sampling schemes for Bayesian MARS models with many predictors." Statistics and Computing 15 (2005): 93-101.
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- adaptive_khaos(X, y)
#' @export
adaptive_khaos2 <-function(X, y,
                          degree=15, order=5,
                          nmcmc=10000,
                          nburn=9000,
                          thin=1,
                          max_basis=1000,
                          a_g = 1e-3, b_g=1e3, zeta=1,
                          g2_sample="mh",
                          g2_init=NULL,
                          s2_lower=0,
                          a_sigma=0,b_sigma=0,
                          a_M=4,b_M=4/length(y),
                          move_probs=rep(1/3, 3),
                          coin_pars=list(function(j) 1/j, 1, 2, 3),
                          degree_penalty=0,
                          legacy=TRUE,
                          verbose=TRUE
){
  # Define a safer solve
  # Note: Under iid uniform X, columns of B are uncorrelated.
  #       So the tolerance should be pretty safe unless there's an
  #       extreme amount of correlation in X.
  # safe_solve <- function(A, tol = 1e-9) {
  #   rc <- rcond(A)
  #   if (is.na(rc) || rc < tol) return(FALSE)
  #   solve(A)
  # }
  #
  safe_solve <- function(A, tol = 1e-9) {
    rc <- rcond(A)
    if (is.na(rc) || rc < tol) return(FALSE)
    chol_A <- tryCatch(chol(A), error = function(e) return(FALSE))
    if (isFALSE(chol_A)) return(FALSE)
    return(chol2inv(chol_A))
  }

  if(max(X) > 1 | min(X) < 0) warning("Inputs are expected to be scaled on (0, 1). Is this intentional?")

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
  s2[1]<-stats::var(y) + 1e-9
  lam[1]<-1
  B.curr<-matrix(rep(1,n)) # matrix of current basis functions, so that yhat = B.curr %*% beta
  BtB.curr <- crossprod(B.curr)
  #BtBi.curr <- solve(BtB.curr)
  #==================================================
  # Modified G-Prior Stuff
  g2 <- rep(NA, nmcmc)
  if(is.null(g2_init)){
    g2[1] <- b_g / (a_g + 1) # Prior mode
  }else{
    g2[1] <- g2_init
  }
  g_vec.curr <- 1          # Intercept only right now
  G.curr <- 1 + tcrossprod(1/g_vec.curr) / g2[1]
  #==================================================
  Sigma.curr <- solve(BtB.curr * G.curr)
  v.curr <- crossprod(B.curr, y)
  quad.curr <- ssy - t(v.curr) %*% Sigma.curr %*% v.curr
  Q.curr <- b_sigma + 0.5 * quad.curr
  ldet.curr <- determinant(Sigma.curr, logarithm = TRUE)$modulus

  #Vinv.curr<-crossprod(B.curr)+1/tau2 # V^(-1) from DMS
  #bhat<-solve(Vinv.curr)%*%t(B.curr)%*%y
  #d.curr <- b_sigma + ssy - t(bhat)%*%Vinv.curr%*%bhat # d from DMS

  count_accept<-c(0,0,0,0)  # count how many times we accept birth, death, change degree, change vars
  count_propose<-c(0,0,0,0) # count how many times we propose each step
  names(count_accept) <- c("Birth", "Death", "Mutate (degree)", "Mutate (vars)")
  names(count_propose) <- c("Birth", "Death", "Mutate (degree)", "Mutate (vars)")
  beta[1,1] <- mean(y)

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

      #==================================================
      # Modified G-Prior Stuff
      qq <- nint.cand
      dd <- dtot.cand
      g_vec.cand <- c(g_vec.curr,
                      (1 + qq*(qq + dd - 2))^(-zeta/2))
      G.cand <- 1 + tcrossprod(1/g_vec.cand) / g2[i-1]

      #==================================================

      if(nint.cand == 0) stop("Why is nint.cand zero? This shouldn't be possible.")
      basis.cand <- make_basis(vars.cand, degs.cand, X)
      B.cand <- cbind(B.curr,basis.cand) # add the new basis function to the basis functions we already have
      BtB.cand <- crossprod(B.cand)
      v.cand <- crossprod(B.cand, y)
      Sigma.cand <- tryCatch({
        safe_solve(G.cand * BtB.cand)
      }, error = function(e) {
        FALSE
      })

      # BtBi.cand <- tryCatch({
      #   safe_solve(BtB.cand)
      # }, error = function(e) {
      #   FALSE
      # })


      if(!isFALSE(Sigma.cand)){ # Otherwise reject immediately
        #Vinv.cand<-crossprod(B.cand)+diag(nbasis[i-1]+2)/tau2 # +2: one for intercept and one for birth
        #bhat.cand<-solve(Vinv.cand)%*%t(B.cand)%*%y
        #d.cand <- b_sigma + ssy - t(bhat.cand)%*%Vinv.cand%*%bhat.cand

        quad.cand <- ssy - t(v.cand) %*% Sigma.cand %*% v.cand
        Q.cand <- b_sigma + 0.5 * quad.cand
        ldet.cand <- determinant(Sigma.cand, logarithm = TRUE)$modulus

        llik.alpha <- as.numeric(
          0.5 * (ldet.cand - ldet.curr) +
            (a_sigma + n/2) * (log(Q.curr) - log(Q.cand))
        )

        # calculate the log likelihood ratio (sort of) after integrating out beta and s2
        #llik.alpha <- (
        #  0.5*log(1/tau2) +
        #    (determinant(Vinv.curr)$mod/2 - determinant(Vinv.cand)$mod/2) +
        #    (a_sigma+n/2)*(log(d.curr) - log(d.cand))
        #)

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
          #BtBi.curr <- BtBi.cand
          #bhat<-bhat.cand
          #Vinv.curr<-Vinv.cand
          #d.curr<-d.cand
          v.curr <- v.cand
          G.curr <- G.cand
          g_vec.curr <- g_vec.cand
          Sigma.curr <- Sigma.cand
          Q.curr <- Q.cand
          ldet.curr <- ldet.cand


          nbasis[i]<-nbasis[i-1]+1
          nint[i,nbasis[i]]<-nint.cand
          dtot[i,nbasis[i]]<-dtot.cand
          if(!legacy){
            vars[i,nbasis[i],] <- NA
            degs[i,nbasis[i],] <- NA
          }
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
      g_vec.cand <- g_vec.curr[-(tokill+1)]
      G.cand <- 1 + tcrossprod(1/g_vec.cand) / g2[i-1]
      v.cand <- crossprod(B.cand, y)
      Sigma.cand <- tryCatch({
        safe_solve(G.cand * BtB.cand)
      }, error = function(e) {
        FALSE
      })

      # BtBi.cand <- tryCatch({
      #   safe_solve(BtB.cand)
      # }, error = function(e) {
      #   FALSE
      # })

      if(!isFALSE(Sigma.cand)){ #Otherwise reject immediately
        #Vinv.cand<-crossprod(B.cand)+diag(nbasis[i-1])/tau2
        #bhat.cand<-solve(Vinv.cand)%*%crossprod(B.cand,y)
        #d.cand <- b_sigma + ssy - crossprod(bhat.cand,Vinv.cand%*%bhat.cand)

        quad.cand <- ssy - t(v.cand) %*% Sigma.cand %*% v.cand
        Q.cand <- b_sigma + 0.5 * quad.cand
        ldet.cand <- determinant(Sigma.cand, logarithm = TRUE)$modulus

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
        #llik.alpha <- (
        #  -0.5*log(1/tau2) +
        #    (determinant(Vinv.curr)$mod/2 - determinant(Vinv.cand)$mod/2) +
        #    (a_sigma+n/2)*(log(d.curr) - log(d.cand))
        #)

        llik.alpha <- as.numeric(
          0.5 * (ldet.cand - ldet.curr) +
            (a_sigma + n/2) * (log(Q.curr) - log(Q.cand))
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
          #BtBi.curr <- BtBi.cand
          #bhat<-bhat.cand
          #Vinv.curr<-Vinv.cand
          #d.curr<-d.cand
          v.curr <- v.cand
          G.curr <- G.cand
          g_vec.curr <- g_vec.cand
          Sigma.curr <- Sigma.cand
          Q.curr <- Q.cand
          ldet.curr <- ldet.cand

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
        tochange<-sample(nbasis[i-1],1) # which basis function we will change

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

        qq <- nint.curr # This is how its defined above
        dd <- dtot.cand
        g_vec.cand <- g_vec.curr
        g_vec.cand[tochange + 1] <- (1 + qq*(qq + dd - 2))^(-zeta/2)
        G.cand <- 1 + tcrossprod(1/g_vec.cand) / g2[i-1]
        v.cand <- crossprod(B.cand, y)
        Sigma.cand <- tryCatch({
          safe_solve(G.cand * BtB.cand)
        }, error = function(e) {
          FALSE
        })

        #BtBi.cand <- tryCatch({
        #  safe_solve(BtB.cand)
        #}, error = function(e) {
        #  FALSE
        #})

        if(!isFALSE(Sigma.cand)){ #Otherwise reject immediately

          quad.cand <- ssy - t(v.cand) %*% Sigma.cand %*% v.cand
          Q.cand <- b_sigma + 0.5 * quad.cand
          ldet.cand <- determinant(Sigma.cand, logarithm = TRUE)$modulus

          # Compute acceptance ratio
          # Assume that "interaction" between the two mutation types is negligible
          llik.alpha <- as.numeric(
            0.5 * (ldet.cand - ldet.curr) +
              (a_sigma + n/2) * (log(Q.curr) - log(Q.cand))
          )

          # Variation in the degree proposal
          lprop.alpha <- log(d_probs[1+dtot.curr-nint.curr]) - log(d_probs[1+dtot.cand-nint.curr]) +
            (-lchoose(dtot.curr, nint.curr) + lchoose(dtot.cand, nint.curr))

          # Uniform priors + no dimension change <=> no term for the prior
          alpha <- llik.alpha + lprop.alpha

          if(log(stats::runif(1))<alpha){ # accept, update
            B.curr<-B.cand
            BtB.curr <- BtB.cand
            #BtBi.curr <- BtBi.cand
            #bhat<-bhat.cand
            #Vinv.curr<-Vinv.cand
            #d.curr<-d.cand
            v.curr <- v.cand
            G.curr <- G.cand
            g_vec.curr <- g_vec.cand
            Sigma.curr <- Sigma.cand
            Q.curr <- Q.cand
            ldet.curr <- ldet.cand

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
        g_vec.cand <- g_vec.curr # No dimension change
        G.cand <- G.curr         # No dimension change
        v.cand <- crossprod(B.cand, y)
        Sigma.cand <- tryCatch({
          safe_solve(G.cand * BtB.cand)
        }, error = function(e) {
          FALSE
        })

        if(!isFALSE(Sigma.cand)){ # Otherwise reject immediately
          quad.cand <- ssy - t(v.cand) %*% Sigma.cand %*% v.cand
          Q.cand <- b_sigma + 0.5 * quad.cand
          ldet.cand <- determinant(Sigma.cand, logarithm = TRUE)$modulus

          # Compute acceptance ratio
          # Assume that "interaction" between the two mutation types is negligible
          llik.alpha <- as.numeric(
            0.5 * (ldet.cand - ldet.curr) +
              (a_sigma + n/2) * (log(Q.curr) - log(Q.cand))
          )

          # Difference in probability of the variables.
          lprop.alpha <- log(newvar_probs.cand[newvar.cand]) - log(newvar_probs.curr[newvar.curr])

          # Uniform priors + no dimension change <=> no term for the prior
          alpha <- llik.alpha + lprop.alpha

          if(log(stats::runif(1))<alpha){ # accept, update
            B.curr<-B.cand
            BtB.curr <- BtB.cand
            #BtBi.curr <- BtBi.cand
            #bhat<-bhat.cand
            #Vinv.curr<-Vinv.cand
            #d.curr<-d.cand
            v.curr <- v.cand
            G.curr <- G.cand
            g_vec.curr <- g_vec.cand
            Sigma.curr <- Sigma.cand
            Q.curr <- Q.cand
            ldet.curr <- ldet.cand

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

    #================
    #  GIBB STEPS
    #================
    # Sample beta
    mu_n <- Sigma.curr %*% v.curr
    beta_curr <- as.vector(mvtnorm::rmvnorm(1,
                                            mu_n,
                                            s2[i-1] * Sigma.curr))
    beta[i,1:(nbasis[i]+1)] <- beta_curr

    # Residuals / variance
    yhat_curr <- B.curr %*% beta_curr
    resid <- y - yhat_curr
    sum_sq[i] <- sum(resid^2)

    # Lambda and s2
    s2[i] <- max(s2_lower,
                 1/stats::rgamma(1, n/2+a_sigma, rate = b_sigma+.5*sum_sq[i]))
    lam[i] <- stats::rgamma(1,a_M+nbasis[i],b_M+1)                                               # update lambda

    # Sample g2
    if(g2_sample == "f"){
      g2[i] <- g2[i-1]
    }else if(g2_sample == "lf"){
      g2[i] <- rg0sq_laplace_full(1, a_g, b_g, g_vec.curr, BtB.curr)
    }else if(g2_sample == "lo"){
      g2[i] <- rg0sq_laplace_orth(1, a_g, b_g, g_vec.curr)
    }else if(g2_sample == "mh"){
      # Metropolis hastings using full Laplace as proposal
      lap_pars <- rg0sq_laplace_full(-1, a_g, b_g, g_vec.curr, BtB.curr)
      g2_cand <- 1/stats::rgamma(1, lap_pars$a, lap_pars$b)

      lp_cand <- log(dgsq_full(g2_cand, a_g, b_g, g_vec.curr, BtB.curr))
      lp_curr <- log(dgsq_full(g2[i-1], a_g, b_g, g_vec.curr, BtB.curr))

      lprop_cand <- dgamma(1 / g2_cand, lap_pars$a, lap_pars$b, log=TRUE) - 2 * log(g2_cand)
      lprop_curr <- dgamma(1 / g2[i-1], lap_pars$a, lap_pars$b, log=TRUE) - 2 * log(g2[i-1])

      log_alpha <- (lp_cand - lp_curr) + (lprop_curr - lprop_cand)
      if(log(runif(1)) < log_alpha){
        g2[i] <- g2_cand
      }else{
        g2[i] <- g2[i-1]
      }
    }else if(g2_sample == "mho"){
      # Metropolis hastings using full Laplace as proposal
      lap_pars <- rg0sq_laplace_orth(-1, a_g, b_g, g_vec.curr)
      g2_cand <- 1/stats::rgamma(1, lap_pars$a, lap_pars$b)

      lp_cand <- log(dgsq_full(g2_cand, a_g, b_g, g_vec.curr, BtB.curr))
      lp_curr <- log(dgsq_full(g2[i-1], a_g, b_g, g_vec.curr, BtB.curr))

      lprop_cand <- dgamma(1 / g2_cand, lap_pars$a, lap_pars$b, log=TRUE) - 2 * log(g2_cand)
      lprop_curr <- dgamma(1 / g2[i-1], lap_pars$a, lap_pars$b, log=TRUE) - 2 * log(g2[i-1])

      log_alpha <- (lp_cand - lp_curr) + (lprop_curr - lprop_cand)
      if(log(runif(1)) < log_alpha){
        g2[i] <- g2_cand
      }else{
        g2[i] <- g2[i-1]
      }
    }else if(g2_sample == "mhoo"){
      # Metropolis hastings using full Laplace as proposal
      lap_pars <- rg0sq_laplace_orth(-1, a_g, b_g, g_vec.curr)
      g2_cand <- 1/stats::rgamma(1, lap_pars$a, lap_pars$b)

      lp_cand <- log(dgsq_orth(g2_cand, a_g, b_g, g_vec.curr))
      lp_curr <- log(dgsq_orth(g2[i-1], a_g, b_g, g_vec.curr))

      lprop_cand <- dgamma(1 / g2_cand, lap_pars$a, lap_pars$b, log=TRUE) - 2 * log(g2_cand)
      lprop_curr <- dgamma(1 / g2[i-1], lap_pars$a, lap_pars$b, log=TRUE) - 2 * log(g2[i-1])

      log_alpha <- (lp_cand - lp_curr) + (lprop_curr - lprop_cand)
      if(log(runif(1)) < log_alpha){
        g2[i] <- g2_cand
      }else{
        g2[i] <- g2[i-1]
      }
    }else{
      warning("g2_sample not recognized. Fixing g2_sample.")
      g2[i] <- g2[i-1]
    }

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
  s2  <- s2[mcmc_iter]
  lam <- lam[mcmc_iter]
  g2  <- g2[mcmc_iter]


  out <- list(B=B.curr,
              vars=vars,degs=degs,
              nint=nint,dtot=dtot,
              nbasis=nbasis,beta=beta,
              s2=s2,lam=lam,g2=g2,sum_sq=sum_sq,
              eta=eta_vec,
              count_accept=count_accept,
              count_propose=count_propose,
              X=X, y=y)
  class(out) <- "adaptive_khaos"
  return(out)
}
