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
  generate <- function(current_set, current_sum, current_nonzero, index) {
    if (current_sum > d || current_nonzero > q) return()
    if (index > p) {
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





