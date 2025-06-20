legendre_poly <- function(x, j){
  if(j == 0){
    return(rep(1, length(x)))
  }
  if(j == 1){
    return(x)
  }
  n <- j - 1
  res <- ((2*n+1)*x*legendre_poly(x, n) - n*legendre_poly(x, n-1))/(n+1)
  return(res)
}

ss_legendre_poly <- function(x, j){
  sqrt(2*j+1)*legendre_poly(2*x-1, j)
}

generate_A <- function(p, d, q) {
  # Initialize the result matrix
  res <- NULL
  #res <- list()

  # Recursive function
  generate <- function(current_set, current_sum, current_nonzero, index){
    if (current_sum > d || current_nonzero > q) return()
    if (index > p){
      if(current_nonzero > 0){
        res <<- rbind(res, current_set)
        #res <- c(res, list(current_set))
      }
      return()
    }

    for (value in 0:(d - current_sum)) {
      new_set <- current_set
      new_set[index] <- value
      generate(new_set,
               current_sum + value,
               current_nonzero + as.numeric(value > 0),
               index + 1)
    }
  }

  # Start the generation process
  generate(integer(p), 0, 0, 1)
  return(res)
}

update_BtB <- function(B, BtB, BtBi, bnew){
  k <- nrow(BtB)
  BtB_new <- BtBi_new <- matrix(NA, nrow=k+1, ncol=k+1)

  # Update BtB
  Btb <- t(B[,1:k])%*%as.matrix(bnew, ncol=1)
  BtB_new[1:k, 1:k] <- BtB
  BtB_new[k+1,1:k]  <- Btb
  BtB_new[1:k,k+1]  <- Btb
  BtB_new[k+1,k+1]  <- sum(bnew^2)

  # Update BtB inverse
  A <- BtBi%*%Btb
  S <- as.numeric(sum(bnew^2) - t(Btb)%*%A)
  BtBi_new[1:k,1:k] <- BtBi + tcrossprod(A)/S
  BtBi_new[1:k,k+1] <- -A/S
  BtBi_new[k+1,1:k] <- -A/S
  BtBi_new[k+1,k+1] <- 1/S

  out <- list(BtB=BtB_new, BtBi=BtBi_new)
  return(out)
}



estimate_map <- function(y, phi, Cinv, sig0=1, tol=1e-3, max_iter=1000){
  n <- length(y)
  p <- ncol(phi)
  sig_curr <- sig0
  Cinv <- Cinv
  flag <- TRUE
  cnt <- 1
  while(flag){
    Cinv_curr <-  crossprod(phi)/sig_curr^2 + Cinv
    tmp <- tryCatch({
      solve(Cinv_curr)
    }, error=function(e){ NULL })
    if(is.null(tmp)){
      #Cinv_curr <- as.matrix(Matrix::nearPD(Cinv_curr)$mat)
      stop("Matrix can't be inverted. Is regularize = TRUE set?")
    }
    C_curr <- solve(Cinv_curr)
    a_curr <- (tcrossprod(C_curr, phi)%*%y)/sig_curr^2
    yhat   <- phi%*%a_curr
    sig_new <- sqrt(mean((y-yhat)^2))
    delta <- abs(sig_new-sig_curr)
    sig_curr <- sig_new

    if(delta < tol){
      flag <- FALSE
    }
    if(cnt >= max_iter){
      flag <- FALSE
    }
    cnt <- cnt + 1
  }
  out <- list(a=a_curr, sig=sig_curr, Caa_inv=Cinv_curr)
  return(out)
}


A_size <- function(p, d, q){
  res <- 0
  for(qq in 1:q){
    for(dd in 1:d){
      res <- res + choose(p, qq)*choose(dd-1, qq-1)
    }
  }
  return(res)
}

# HELPERS FOR ADAPTIVE VERSION

make_weights <- function(eta, p0, epsilon, alpha, num_passes){
  p <- length(eta)
  v <- (eta^alpha + epsilon)/sum((eta^alpha + epsilon))*p0/p
  delta <- 0
  for(i in 1:num_passes){
    beta <- mean((log(p0 + delta) - log(p))/log(v))
    delta <- delta + p0 - sum(v^beta)
  }
  return(v^beta)
}

# Breaks d into q elements d_1, ... d_q >= 1 s.t d1+...+dq = d.
# Probability of getting any particular sequence is 1/choose(d-1, q-1)
random_partition <- function(d, q) {
  if (d < q) stop("d must be greater than or equal to q")
  partition <- rep(1, q)
  remaining <- d - q

  # Intuitive but slower
  #extra     <- sample(factor(1:q), remaining, replace=TRUE)
  #partition <- partition + as.numeric(table(extra))

  # Less intuitive but faster
  cuts <- sort(sample(0:remaining, q - 1, replace = TRUE))
  additions <- diff(c(0, cuts, remaining))
  partition <- partition + additions
  return(partition)
}


make_basis<-function(vv,dd,XX){ # function to make a basis function given vars and degs
  n <- nrow(XX)
  curr <- rep(1, n)
  for(j in seq_along(vv)){
    curr <- curr * ss_legendre_poly(XX[,vv[j]], dd[j])
  }
  return(curr)
}

# Function to sample from multivariate normal distribution
sample_mvn <- function(n, m, S) {
  p <- length(m)
  L <- chol(S)
  Z <- matrix(stats::rnorm(n * p), nrow = n)
  samples <- tcrossprod(Z, L) + matrix(m, nrow = n, ncol = length(m), byrow = TRUE)
  return(samples)
}

myTimestamp <-function(){
  x<-Sys.time()
  paste('#--',format(x,"%b %d %X"),'--#')
}

test_iid_unif <- function(X, Nsims=1000, alpha=0.01){
  # Test uniformity
  p_unif <- min(apply(X, 2, function(xx) stats::ks.test(xx, "punif")$p.value))

  # Test pairwise independence
  nn <- nrow(X)
  pp <- ncol(X)
  NN <- Nsims

  S0 <- diag(rep(1/12), pp)
  Sx <- stats::cov(X)
  Tx <- sum((S0-Sx)^2)
  T0 <- rep(NA, NN)
  for(ii in 1:NN){
    Xi <- matrix(stats::runif(nn*pp), nrow=nn, ncol=pp)
    Si <- stats::cov(Xi)
    T0[ii] <- sum((S0-Si)^2)
  }
  p_ind <- mean(T0 > Tx)


  T_obs <- norm(stats::cov(X) - (1/12) * diag(pp), type = "F")^2
  T_null <- replicate(NN, {
    norm(stats::cov(matrix(stats::runif(nn * pp), nn, pp)) - (1/12) * diag(pp), type = "F")^2
  })
  p_val <- mean(T_null > T_obs)  # one-sided

  if(min(p_ind, p_unif) < alpha){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

