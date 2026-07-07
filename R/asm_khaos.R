#' Compute Active Subspace C Matrices for adaptive_khaos Objects
#'
#' @param mod A fitted adaptive_khaos object.
#' @param prior Optional prior list. Currently only used for native-scale transforms.
#' @param mcmc.use Vector of MCMC sample indices to use.
#' @param use_native_scale Logical. If TRUE, transform C back to native input scale.
#' @param func.use Included for compatibility with C_bass. Currently ignored.
#' @param func.true Included for compatibility with C_bass. Currently ignored.
#' @param pca_method Included for compatibility with C_bass. Currently ignored.
#' @param verbose Logical. Should progress be displayed?
#'
#' @return A list of C matrices, one per selected MCMC sample.
#' @export
C_khaos <- function(
    mod,
    prior = NULL,
    mcmc.use = NULL,
    use_native_scale = FALSE,
    func.use = NULL,
    func.true = NULL,
    pca_method = 2,
    verbose = FALSE
) {

  if (!inherits(mod, "adaptive_khaos")) {
    stop("mod must be an object of class 'adaptive_khaos'.")
  }

  if (is.null(mcmc.use)) {
    mcmc.use <- seq_along(mod$nbasis)
  }

  p <- ncol(mod$X)

  max_degree <- max(mod$degs, na.rm = TRUE)
  if (!is.finite(max_degree)) max_degree <- 0

  legendre_tabs <- make_khaos_integral_tables(max_degree)
  A <- legendre_tabs$A
  B <- legendre_tabs$B

  out <- vector("list", length(mcmc.use))

  for (ii in seq_along(mcmc.use)) {
    s <- mcmc.use[ii]

    if (verbose) {
      message("Computing C for MCMC sample ", ii, " of ", length(mcmc.use))
    }

    C <- C_khaos_one_sample(
      mod = mod,
      s = s,
      p = p,
      A = A,
      B = B
    )

    if (use_native_scale) {
      scale_vec <- get_khaos_native_scale(mod, prior, p)
      D <- diag(scale_vec, nrow = p)
      C <- D %*% C %*% D
    }

    out[[ii]] <- C
  }

  out
}


C_khaos_one_sample <- function(mod, s, p, A, B) {

  M <- mod$nbasis[s]

  C <- matrix(0, nrow = p, ncol = p)

  if (M == 0) {
    return(C)
  }

  alpha <- matrix(0L, nrow = M, ncol = p)

  for (m in seq_len(M)) {
    q <- mod$nint[s, m]

    if (is.na(q) || q == 0) next

    vars_m <- mod$vars[s, m, seq_len(q)]
    degs_m <- mod$degs[s, m, seq_len(q)]

    keep <- !is.na(vars_m) & !is.na(degs_m)

    vars_m <- vars_m[keep]
    degs_m <- degs_m[keep]

    alpha[m, vars_m] <- degs_m
  }

  beta <- mod$beta[s, seq_len(M + 1)]
  beta <- beta[-1]  # remove intercept

  for (i in seq_len(p)) {
    for (m1 in seq_len(M)) {
      for (m2 in seq_len(M)) {

        if (p > 1 && any(alpha[m1, -i] != alpha[m2, -i])) next

        a <- alpha[m1, i]
        b <- alpha[m2, i]

        C[i, i] <- C[i, i] +
          beta[m1] * beta[m2] * B[a + 1, b + 1]
      }
    }
  }

  if (p >= 2) {
    for (i in seq_len(p - 1)) {
      for (j in (i + 1):p) {

        keep <- setdiff(seq_len(p), c(i, j))
        val <- 0

        for (m1 in seq_len(M)) {
          for (m2 in seq_len(M)) {

            if (length(keep) > 0 &&
                any(alpha[m1, keep] != alpha[m2, keep])) next

            ai1 <- alpha[m1, i]
            ai2 <- alpha[m2, i]
            aj1 <- alpha[m1, j]
            aj2 <- alpha[m2, j]

            val <- val +
              beta[m1] * beta[m2] *
              A[ai1 + 1, ai2 + 1] *
              A[aj2 + 1, aj1 + 1]
          }
        }

        C[i, j] <- val
        C[j, i] <- val
      }
    }
  }

  C
}


make_khaos_integral_tables <- function(max_degree) {

  A <- matrix(0, nrow = max_degree + 1, ncol = max_degree + 1)
  B <- matrix(0, nrow = max_degree + 1, ncol = max_degree + 1)

  for (a in 0:max_degree) {
    for (b in 0:max_degree) {

      ia <- a + 1
      ib <- b + 1

      # A_ab = int_0^1 psi_a'(x) psi_b(x) dx
      if (a > b && ((a - b) %% 2 == 1)) {
        A[ia, ib] <- 2 * sqrt((2 * a + 1) * (2 * b + 1))
      }

      # B_ab = int_0^1 psi_a'(x) psi_b'(x) dx
      if (a >= 1 && b >= 1 && ((a - b) %% 2 == 0)) {
        r <- min(a, b)
        B[ia, ib] <- 2 * sqrt((2 * a + 1) * (2 * b + 1)) * r * (r + 1)
      }
    }
  }

  list(A = A, B = B)
}


get_khaos_native_scale <- function(mod, prior, p) {

  if (!is.null(prior)) {
    if (!is.null(prior$lower) && !is.null(prior$upper)) {
      return(1 / (prior$upper - prior$lower))
    }

    if (!is.null(prior$min) && !is.null(prior$max)) {
      return(1 / (prior$max - prior$min))
    }

    if (!is.null(prior$bounds)) {
      return(1 / (prior$bounds[, 2] - prior$bounds[, 1]))
    }
  }

  stop(
    "use_native_scale = TRUE requires a prior list with lower/upper, min/max, or bounds."
  )
}
