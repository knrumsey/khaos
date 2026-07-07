# Joint Bayesian Adaptive Polynomial Chaos Expansion with ridge prior
#
# This is a first/reference implementation of the joint-output KHAOS idea:
#   Y_{ik} = alpha_k + sum_m a_m * 1(k in r_m) * h_m(x_i) + e_{ik},
#   e_{ik} ~ N(0, sigma_k^2),
# with one shared coefficient a_m per non-intercept basis and response-specific
# intercepts alpha_k. Coefficients are collapsed for RJ moves and sampled in a
# Gibbs step. Basis subsets r_m are proposed with an NKD-style sparse hash-map
# proposal.

# -----------------------------------------------------------------------------
# Small helpers kept local / prefixed to avoid clashing with existing code.
# -----------------------------------------------------------------------------

.jk_timestamp <- function() {
  x <- Sys.time()
  paste("#--", format(x, "%b %d %X"), "--#")
}

.jk_logspace_sum <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

.jk_log1mexp <- function(a) {
  # log(1 - exp(a)) for a <= 0
  if (a > 0) stop("a must be <= 0")
  if (a < log(0.5)) log1p(-exp(a)) else log(-expm1(a))
}

.jk_legendre_poly <- function(x, j) {
  if (j == 0L) return(rep(1, length(x)))
  if (j == 1L) return(x)
  n <- j - 1L
  ((2 * n + 1) * x * .jk_legendre_poly(x, n) - n * .jk_legendre_poly(x, n - 1L)) / (n + 1L)
}

.jk_ss_legendre_poly <- function(x, j) {
  sqrt(2 * j + 1) * .jk_legendre_poly(2 * x - 1, j)
}

.jk_make_basis <- function(vv, dd, XX) {
  n <- nrow(XX)
  curr <- rep(1, n)
  if (length(vv) == 0L) return(curr)
  for (j in seq_along(vv)) {
    curr <- curr * .jk_ss_legendre_poly(XX[, vv[j]], dd[j])
  }
  curr
}

.jk_A_size <- function(p, d, q) {
  res <- 0
  for (qq in seq_len(q)) {
    for (dd in seq_len(d)) {
      if (dd >= qq) res <- res + choose(p, qq) * choose(dd - 1, qq - 1)
    }
  }
  res
}

.jk_make_weights <- function(eta, p0, epsilon, alpha, num_passes) {
  p <- length(eta)
  if (is.infinite(epsilon)) {
    return(rep(min(max(p0 / p, 1e-12), 1 - 1e-12), p))
  }
  raw <- eta^alpha + epsilon
  if (!all(is.finite(raw)) || sum(raw) <= 0) raw <- rep(1, p)
  v <- raw / sum(raw) * p0 / p
  v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
  delta <- 0
  beta <- 1
  for (ii in seq_len(num_passes)) {
    denom <- log(v)
    beta <- mean((log(p0 + delta) - log(p)) / denom)
    if (!is.finite(beta)) beta <- 1
    delta <- delta + p0 - sum(v^beta)
  }
  out <- v^beta
  pmin(pmax(out, 1e-12), 1 - 1e-12)
}

.jk_random_partition <- function(d, q) {
  # Uniform positive integer composition of d into q pieces.
  if (d < q) stop("d must be >= q")
  if (q == 1L) return(d)
  cuts <- sort(sample(seq_len(d - 1L), q - 1L, replace = FALSE))
  diff(c(0L, cuts, d))
}

.jk_poibin_valid_prob <- function(w, max_count) {
  # Probability 1 <= sum Bernoulli(w_j) <= max_count.
  max_count <- min(max_count, length(w))
  pmf <- numeric(max_count + 1L)
  pmf[1L] <- 1
  for (prob in w) {
    new <- numeric(max_count + 1L)
    new[1L] <- pmf[1L] * (1 - prob)
    if (max_count >= 1L) {
      for (cc in seq_len(max_count)) {
        new[cc + 1L] <- pmf[cc + 1L] * (1 - prob) + pmf[cc] * prob
      }
    }
    pmf <- new
  }
  sum(pmf[-1L])
}

.jk_log_prob_vars_coin <- function(vars, eta, order, J_probs, coin_pars) {
  p <- length(eta)
  vars <- sort(unique(as.integer(vars)))
  if (length(vars) < 1L || length(vars) > order) return(-Inf)
  if (any(vars < 1L | vars > p)) return(-Inf)

  J_probs <- J_probs / sum(J_probs)
  log_terms <- rep(-Inf, order)
  log_valid_terms <- rep(-Inf, order)

  for (jj in seq_len(order)) {
    wts <- .jk_make_weights(eta, jj, coin_pars[[3]], coin_pars[[2]], coin_pars[[4]])
    wts <- pmin(pmax(wts, 1e-12), 1 - 1e-12)

    not_vars <- setdiff(seq_len(p), vars)
    log_p_subset <- sum(log(wts[vars]))
    if (length(not_vars) > 0L) log_p_subset <- log_p_subset + sum(log1p(-wts[not_vars]))
    log_terms[jj] <- log(J_probs[jj]) + log_p_subset

    valid_prob <- .jk_poibin_valid_prob(wts, order)
    if (valid_prob > 0) log_valid_terms[jj] <- log(J_probs[jj]) + log(valid_prob)
  }

  .jk_logspace_sum(log_terms) - .jk_logspace_sum(log_valid_terms)
}

.jk_basis_log_q <- function(vars, degs, eta, degree, order, J_probs, coin_pars, degree_penalty) {
  vars <- as.integer(vars)
  degs <- as.integer(degs)
  nint <- length(vars)
  if (nint < 1L || nint > order || length(degs) != nint) return(-Inf)
  if (any(degs < 1L)) return(-Inf)
  dtot <- sum(degs)
  if (dtot < nint || dtot > degree) return(-Inf)

  d_probs <- seq_along(nint:degree)^(-degree_penalty)
  d_probs <- d_probs / sum(d_probs)
  idx <- dtot - nint + 1L

  .jk_log_prob_vars_coin(vars, eta, order, J_probs, coin_pars) +
    log(d_probs[idx]) -
    lchoose(dtot - 1L, nint - 1L)
}

.jk_propose_basis <- function(X, eta, degree, order, J_probs, coin_pars, degree_penalty) {
  p <- ncol(X)
  repeat {
    q0 <- sample(seq_len(order), 1L, prob = J_probs)
    wts <- .jk_make_weights(eta, q0, coin_pars[[3]], coin_pars[[2]], coin_pars[[4]])
    chi <- stats::rbinom(p, 1L, wts)
    nint <- sum(chi)
    if (nint >= 1L && nint <= order) break
  }
  vars <- which(chi == 1L)
  nint <- length(vars)
  d_probs <- seq_along(nint:degree)^(-degree_penalty)
  d_probs <- d_probs / sum(d_probs)
  dtot <- sample(nint:degree, 1L, prob = d_probs)
  degs <- .jk_random_partition(dtot, nint)
  basis <- .jk_make_basis(vars, degs, X)
  log_q <- .jk_basis_log_q(vars, degs, eta, degree, order, J_probs, coin_pars, degree_penalty)
  list(vars = vars, degs = degs, nint = nint, dtot = dtot, basis = basis, log_q = log_q)
}

.jk_newvar_probs <- function(eta, forbidden, epsilon) {
  p <- length(eta)
  if (is.infinite(epsilon)) {
    probs <- rep(1, p)
  } else {
    probs <- 1e-9 + epsilon + pmax(eta, 0)
  }
  if (length(forbidden) > 0L) probs[forbidden] <- 0
  if (sum(probs) <= 0) return(NULL)
  probs / sum(probs)
}

.jk_key_from_vec <- function(r) {
  r <- sort(unique(as.integer(r)))
  paste(r, collapse = ",")
}

.jk_vec_from_key <- function(key) {
  if (length(key) == 0L || is.na(key) || !nzchar(key)) return(integer(0))
  as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
}

.jk_subset_counts <- function(rkeys) {
  rkeys <- rkeys[!is.na(rkeys) & nzchar(rkeys)]
  if (length(rkeys) == 0L) return(setNames(integer(0), character(0)))
  tab <- table(rkeys)
  out <- as.integer(tab)
  names(out) <- names(tab)
  out
}

.jk_subset_size_log_probs <- function(K, subset_power) {
  d <- seq_len(K)
  logs <- subset_power * log(d)
  logs - .jk_logspace_sum(logs)
}

.jk_subset_log_q0 <- function(key, K, subset_power_prop) {
  r <- .jk_vec_from_key(key)
  d <- length(r)
  if (d < 1L || d > K) return(-Inf)

  log_size_probs <- .jk_subset_size_log_probs(K, subset_power_prop)

  ## q0(r) = q0(|r|=d) * uniform among subsets of size d
  log_size_probs[d] - lchoose(K, d)
}

.jk_subset_log_prior <- function(key, K, subset_power_prior) {
  r <- .jk_vec_from_key(key)
  d <- length(r)
  if (d < 1L || d > K) return(-Inf)

  log_size_probs <- .jk_subset_size_log_probs(K, subset_power_prior)

  ## pi(r) = pi(|r|=d) * uniform among subsets of size d
  log_size_probs[d] - lchoose(K, d)
}

.jk_subset_log_q_mix <- function(key, K, counts, M, rho, subset_power_prop) {
  log_q0 <- .jk_subset_log_q0(key, K, subset_power_prop)
  if (M <= 0L) return(log_q0)

  cc <- 0L
  if (length(counts) > 0L && key %in% names(counts)) cc <- counts[[key]]

  terms <- c(-Inf, -Inf)
  if (rho > 0) terms[1L] <- log(rho) + log_q0
  if ((1 - rho) > 0 && cc > 0) terms[2L] <- log1p(-rho) + log(cc) - log(M)

  .jk_logspace_sum(terms)
}

.jk_propose_subset <- function(K, counts, M, rho, subset_power_prop) {
  use_base <- (M <= 0L) || (stats::runif(1) < rho) || length(counts) == 0L

  if (use_base) {
    log_size_probs <- .jk_subset_size_log_probs(K, subset_power_prop)
    d <- sample(seq_len(K), 1L, prob = exp(log_size_probs))
    r <- sort(sample(seq_len(K), d, replace = FALSE))
    key <- .jk_key_from_vec(r)
  } else {
    key <- sample(names(counts), 1L, prob = counts)
    r <- .jk_vec_from_key(key)
  }

  log_q <- .jk_subset_log_q_mix(key, K, counts, M, rho, subset_power_prop)

  list(r = r, key = key, log_q = log_q)
}

.jk_build_design <- function(H, rkeys, n, K) {
  M <- if (is.null(H)) 0L else ncol(H)
  Xstk <- matrix(0, nrow = n * K, ncol = K + M)
  for (kk in seq_len(K)) {
    rows <- ((kk - 1L) * n + 1L):(kk * n)
    Xstk[rows, kk] <- 1
  }
  if (M > 0L) {
    for (mm in seq_len(M)) {
      rr <- .jk_vec_from_key(rkeys[mm])
      if (length(rr) > 0L) {
        for (kk in rr) {
          rows <- ((kk - 1L) * n + 1L):(kk * n)
          Xstk[rows, K + mm] <- H[, mm]
        }
      }
    }
  }
  Xstk
}

.jk_log_marginal_ridge <- function(Xstk, yvec, s2, n, K, prior_var, need_posterior = FALSE) {
  N <- length(yvec)
  P <- ncol(Xstk)
  if (length(prior_var) != P) stop("prior_var length does not match number of columns")
  if (any(s2 <= 0) || any(!is.finite(s2))) return(list(loglik = -Inf))

  dinv <- rep(1 / s2, each = n)
  logdetD <- sum(n * log(s2))
  logdetC <- sum(log(prior_var))
  Cinv_diag <- 1 / prior_var

  WX <- sweep(Xstk, 1L, dinv, "*")
  XtDX <- crossprod(Xstk, WX)
  m <- crossprod(Xstk, dinv * yvec)
  Lambda <- XtDX + diag(Cinv_diag, P)

  R <- tryCatch(chol(Lambda), error = function(e) NULL)
  if (is.null(R)) return(list(loglik = -Inf))

  logdetLambda <- 2 * sum(log(diag(R)))
  u <- forwardsolve(t(R), m)
  z <- backsolve(R, u)
  quad <- sum(dinv * yvec^2) - sum(m * z)
  loglik <- -0.5 * (N * log(2 * pi) + logdetD + logdetC + logdetLambda + quad)

  if (!need_posterior) return(list(loglik = as.numeric(loglik)))
  list(loglik = as.numeric(loglik), mu = as.vector(z), R = R)
}

.jk_rmvnorm_precision_chol <- function(mu, R) {
  # R is upper triangular chol(precision), so precision = t(R) %*% R.
  as.vector(mu + backsolve(R, stats::rnorm(length(mu))))
}

.jk_move_prob <- function(type, M, max_basis, move_probs) {
  names(move_probs) <- c("birth", "death", "change")
  if (M <= 0L) return(if (type == "birth") 1 else 0)
  if (M >= max_basis) {
    denom <- move_probs["death"] + move_probs["change"]
    if (type == "birth") return(0)
    return(move_probs[type] / denom)
  }
  move_probs[type] / sum(move_probs)
}

# -----------------------------------------------------------------------------
# Main fitting routine.
# -----------------------------------------------------------------------------

#' Joint adaptive KHAOS with ridge prior and subset-shared bases
#'
#' @param X n by p matrix/data frame of predictors, expected to be scaled to [0,1].
#' @param Y n by K response matrix. A vector is treated as K=1.
#' @param degree Maximum total polynomial degree for a basis.
#' @param order Maximum interaction order for a basis.
#' @param nmcmc Number of MCMC iterations.
#' @param nburn Number of burn-in iterations.
#' @param thin Keep every thin-th post-burn-in sample.
#' @param max_basis Maximum number of non-intercept bases.
#' @param tau2 Prior variance for shared non-intercept coefficients.
#' @param tau2_intercept Prior variance for response-specific intercepts.
#' @param s2_lower Lower bound for sigma_k^2 draws.
#' @param g1,g2 Shape/rate hyperparameters for sigma_k^2 inverse-gamma priors.
#' @param h1,h2 Shape/rate hyperparameters for lambda, the Poisson mean of M.
#' @param move_probs Birth/death/change probabilities.
#' @param subset_rho Exploration probability in q(r)=rho*q0(r)+(1-rho)*c_r/M.
#' @param subset_power_prop Base proposal power for subset size: q0(|r|=d) proportional to d^subset_power_prop.
#' @param subset_power_prior Prior power for subset size: pi(|r|=d) proportional to d^subset_power_prior.
#' @param coin_pars Input-variable coin-flip proposal parameters, as in adaptive_khaos_ridge.
#' @param degree_penalty Increasing this favors lower total degree proposals.
#' @param legacy If FALSE, clears unused storage after accepted moves.
#' @param verbose Print progress messages.
#'
#' @details
#' This function fits a joint-output adaptive polynomial chaos model using a
#' shared basis expansion. For observations \eqn{i = 1,\ldots,n} and outputs
#' \eqn{k = 1,\ldots,K}, the model has the form
#' \deqn{
#'   Y_{ik} = \alpha_k + \sum_{m=1}^M a_m 1(k \in r_m) h_m(x_i)
#'            + \epsilon_{ik},
#'   \qquad \epsilon_{ik} \sim N(0,\sigma_k^2).
#' }
#' Each basis function \eqn{h_m} is assigned to a nonempty subset
#' \eqn{r_m \subseteq \{1,\ldots,K\}}. If \eqn{k \in r_m}, then basis
#' \eqn{m} contributes to output \eqn{k}; otherwise it is inactive for that
#' output. Thus, basis functions may be output-specific, shared by a small
#' group of outputs, or shared globally across all outputs.
#'
#' The coefficients for non-intercept basis functions are shared across all
#' outputs in their assigned subset. This means that a basis assigned to
#' \eqn{r_m = \{1,3\}} contributes the same coefficient-weighted function
#' to outputs 1 and 3. Response-specific intercepts and response-specific
#' variances are still included.
#'
#' Model exploration is performed with reversible-jump MCMC birth, death, and
#' change moves. In a birth move, a new polynomial chaos basis function is
#' proposed together with a subset \eqn{r_m}. Subsets are proposed using a
#' mixture of a base distribution and an adaptive count-based distribution:
#' \deqn{
#'   q(r) = \rho q_0(r) + (1-\rho)\frac{c_r}{M},
#' }
#' where \eqn{c_r} is the number of current basis functions assigned exactly
#' to subset \eqn{r}, \eqn{M} is the current number of non-intercept basis
#' functions, and \eqn{\rho} is controlled by \code{subset_rho}. The count
#' term encourages reuse of subset structures that are already present in
#' the model, while the base proposal ensures continued exploration.
#'
#' The base proposal first samples a subset size \eqn{d = |r|} according to
#' \deqn{
#'   q_0(|r|=d) \propto d^{s_q},
#' }
#' where \eqn{s_q} is \code{subset_power_prop}, and then samples uniformly
#' among all subsets of size \eqn{d}. Thus, the probability of proposing a
#' particular subset \eqn{r} under the base proposal is
#' \deqn{
#'   q_0(r) =
#'   \frac{|r|^{s_q}}{\sum_{d=1}^K d^{s_q}}
#'   {K \choose |r|}^{-1}.
#' }
#'
#' The prior on subset assignments has the same size-based form, but with its
#' own power parameter:
#' \deqn{
#'   \pi(|r|=d) \propto d^{s_\pi},
#' }
#' where \eqn{s_\pi} is \code{subset_power_prior}. Conditional on the size,
#' subsets are again assumed uniform, giving
#' \deqn{
#'   \pi(r) =
#'   \frac{|r|^{s_\pi}}{\sum_{d=1}^K d^{s_\pi}}
#'   {K \choose |r|}^{-1}.
#' }
#'
#' The power parameters may be positive, zero, or negative. Setting a power
#' to zero gives a distribution that is uniform over subset sizes. Positive
#' values favor larger subsets, encouraging more shared structure across
#' outputs. Negative values favor smaller subsets, encouraging more
#' output-specific structure and less sharing. For example,
#' \code{subset_power_prior = 1} mildly favors broadly shared basis functions,
#' while \code{subset_power_prior = -1} favors models with less shared
#' structure.
#'
#' Conditional on the current basis functions, subset assignments, and
#' response variances, the shared coefficients have a Gaussian ridge prior
#' and can be integrated out analytically for birth, death, and change
#' acceptance probabilities. They are then sampled in a Gibbs step for
#' posterior prediction. The residual variances \eqn{\sigma_k^2} are updated
#' independently for each output using conjugate inverse-gamma updates.
#'
#' @return An object of class "adaptive_jkhaos".
#' @export
adaptive_jkhaos_ridge <- function(X, Y,
                                  degree = 15, order = 5,
                                  nmcmc = 10000,
                                  nburn = 9000,
                                  thin = 1,
                                  max_basis = 1000,
                                  tau2 = 1e5,
                                  tau2_intercept = tau2,
                                  s2_lower = 0,
                                  g1 = 0, g2 = 0,
                                  h1 = 4, h2 = 40 / nrow(as.matrix(Y)),
                                  move_probs = rep(1/3, 3),
                                  subset_rho = 0.5,
                                  subset_power_prop = 1,
                                  subset_power_prior = 1,
                                  coin_pars = list(function(j) 1 / j, 1, 2, 3),
                                  degree_penalty = 0,
                                  legacy = TRUE,
                                  verbose = TRUE) {
  subset_prior <- TRUE # deprecated

  if (!is.matrix(X)) {
    if (is.data.frame(X)) {
      X <- as.matrix(X)
    } else if (is.numeric(X)) {
      X <- matrix(X, ncol = 1)
    } else {
      stop("X must be a matrix, data frame, or numeric vector.")
    }
  }
  if (!is.matrix(Y)) {
    if (is.data.frame(Y)) {
      Y <- as.matrix(Y)
    } else if (is.numeric(Y)) {
      Y <- matrix(Y, ncol = 1)
    } else {
      stop("Y must be a matrix, data frame, or numeric vector.")
    }
  }
  storage.mode(X) <- "double"
  storage.mode(Y) <- "double"

  if (max(X) > 1 || min(X) < 0) warning("Inputs are expected to be scaled on [0,1]. Is this intentional?")
  if (nrow(Y) != nrow(X)) stop("nrow(Y) must equal nrow(X).")
  if (nmcmc < 2) stop("nmcmc must be at least 2.")
  if (nburn >= nmcmc) stop("nburn must be smaller than nmcmc.")
  if (thin < 1) stop("thin must be >= 1.")
  if (subset_rho < 0 || subset_rho > 1) stop("subset_rho must be in [0,1].")
  if (!is.finite(subset_power_prop)) {
    stop("subset_power_prop must be finite.")
  }
  if (!is.finite(subset_power_prior)) {
    stop("subset_power_prior must be finite.")
  }
  if (length(move_probs) != 3L || any(move_probs < 0) || sum(move_probs) <= 0) {
    stop("move_probs must be a nonnegative length-3 vector with positive sum.")
  }

  n <- nrow(X)
  p <- ncol(X)
  K <- ncol(Y)
  yvec <- as.vector(Y) # stacks by response: Y[,1], Y[,2], ...

  if (p < order) order <- p
  if (degree < order) stop("degree must be at least as large as order.")
  if (order < 1L) stop("order must be at least 1.")

  J_probs <- coin_pars[[1]](seq_len(order))
  if (length(J_probs) != order || any(J_probs < 0) || sum(J_probs) <= 0) {
    stop("coin_pars[[1]] must return nonnegative weights of length order with positive sum.")
  }
  J_probs <- J_probs / sum(J_probs)
  A_total <- .jk_A_size(p, degree, order)
  log_prior_basis <- -log(A_total)

  # Storage.
  degs <- vars <- array(NA_integer_, dim = c(nmcmc, max_basis, order))
  nint <- dtot <- matrix(NA_integer_, nrow = nmcmc, ncol = max_basis)
  rkeys <- matrix(NA_character_, nrow = nmcmc, ncol = max_basis)

  beta <- matrix(NA_real_, nrow = nmcmc, ncol = K + max_basis)
  s2 <- matrix(NA_real_, nrow = nmcmc, ncol = K)
  colnames(s2) <- colnames(Y)
  lam <- nbasis <- rep(NA_real_, nmcmc)
  sum_sq <- matrix(NA_real_, nrow = nmcmc, ncol = K)
  colnames(sum_sq) <- colnames(Y)
  eta_vec <- rep(0, p)

  count_accept <- c(0, 0, 0, 0)
  count_propose <- c(0, 0, 0, 0)
  names(count_accept) <- c("Birth", "Death", "Mutate (degree)", "Mutate (vars)")
  names(count_propose) <- names(count_accept)

  # Current state.
  nbasis[1] <- 0L
  lam[1] <- 1
  s2_init <- apply(Y, 2L, stats::var)
  s2_init[!is.finite(s2_init) | s2_init <= 0] <- 1e-9
  s2[1, ] <- s2_init + 1e-9
  H.curr <- matrix(nrow = n, ncol = 0L)
  rkeys.curr <- character(0)
  Xstk.curr <- .jk_build_design(H.curr, rkeys.curr, n, K)
  prior_var.curr <- rep(tau2_intercept, K)
  ml.curr <- .jk_log_marginal_ridge(Xstk.curr, yvec, s2[1, ], n, K, prior_var.curr, need_posterior = TRUE)
  if (!is.finite(ml.curr$loglik)) stop("Initial marginal likelihood is not finite.")
  beta[1, 1:K] <- ml.curr$mu

  if (verbose) {
    cat(c("MCMC iteration", 0, .jk_timestamp(), "nbasis:", 0), "\n")
  }

  for (ii in 2:nmcmc) {
    M.prev <- as.integer(nbasis[ii - 1L])
    move.type <- sample(c("birth", "death", "change"), 1L,
                        prob = c(.jk_move_prob("birth", M.prev, max_basis, move_probs),
                                 .jk_move_prob("death", M.prev, max_basis, move_probs),
                                 .jk_move_prob("change", M.prev, max_basis, move_probs)))

    # Carry state forward in storage.
    nbasis[ii] <- M.prev
    lam[ii] <- lam[ii - 1L]
    s2[ii, ] <- s2[ii - 1L, ]
    nint[ii, ] <- nint[ii - 1L, ]
    dtot[ii, ] <- dtot[ii - 1L, ]
    vars[ii, , ] <- vars[ii - 1L, , ]
    degs[ii, , ] <- degs[ii - 1L, , ]
    rkeys[ii, ] <- rkeys[ii - 1L, ]

    if (move.type == "birth") {
      counts.curr <- .jk_subset_counts(rkeys.curr)
      subset.cand <- .jk_propose_subset(K, counts.curr, M.prev, subset_rho, subset_power_prop)
      basis.cand <- .jk_propose_basis(X, eta_vec, degree, order, J_probs, coin_pars, degree_penalty)

      H.cand <- cbind(H.curr, basis.cand$basis)
      rkeys.cand <- c(rkeys.curr, subset.cand$key)
      Xstk.cand <- .jk_build_design(H.cand, rkeys.cand, n, K)
      prior_var.cand <- c(rep(tau2_intercept, K), rep(tau2, M.prev + 1L))
      ml.cand <- .jk_log_marginal_ridge(Xstk.cand, yvec, s2[ii - 1L, ], n, K, prior_var.cand)

      if (is.finite(ml.cand$loglik)) {
        lp_subset <- .jk_subset_log_prior(subset.cand$key, K, subset_power_prior)
        llik.alpha <- ml.cand$loglik - ml.curr$loglik
        lprior.alpha <- log(lam[ii - 1L]) + log_prior_basis + lp_subset
        lprop.alpha <- (log(.jk_move_prob("death", M.prev + 1L, max_basis, move_probs)) - log(M.prev + 1L)) -
          (log(.jk_move_prob("birth", M.prev, max_basis, move_probs)) + basis.cand$log_q + subset.cand$log_q)
        alpha <- llik.alpha + lprior.alpha + lprop.alpha

        if (log(stats::runif(1)) < alpha) {
          H.curr <- H.cand
          rkeys.curr <- rkeys.cand
          Xstk.curr <- Xstk.cand
          ml.curr <- .jk_log_marginal_ridge(Xstk.curr, yvec, s2[ii - 1L, ], n, K, prior_var.cand, need_posterior = TRUE)
          nbasis[ii] <- M.prev + 1L

          nint[ii, nbasis[ii]] <- basis.cand$nint
          dtot[ii, nbasis[ii]] <- basis.cand$dtot
          if (!legacy) {
            vars[ii, nbasis[ii], ] <- NA_integer_
            degs[ii, nbasis[ii], ] <- NA_integer_
          }
          vars[ii, nbasis[ii], seq_len(basis.cand$nint)] <- basis.cand$vars
          degs[ii, nbasis[ii], seq_len(basis.cand$nint)] <- basis.cand$degs
          rkeys[ii, nbasis[ii]] <- subset.cand$key

          eta_vec[basis.cand$vars] <- eta_vec[basis.cand$vars] + 1
          count_accept["Birth"] <- count_accept["Birth"] + 1
        }
      }
      count_propose["Birth"] <- count_propose["Birth"] + 1

    } else if (move.type == "death") {
      tokill <- sample(seq_len(M.prev), 1L)
      keep <- setdiff(seq_len(M.prev), tokill)

      H.cand <- H.curr[, keep, drop = FALSE]
      rkeys.cand <- rkeys.curr[keep]
      Xstk.cand <- .jk_build_design(H.cand, rkeys.cand, n, K)
      prior_var.cand <- c(rep(tau2_intercept, K), rep(tau2, M.prev - 1L))
      ml.cand <- .jk_log_marginal_ridge(Xstk.cand, yvec, s2[ii - 1L, ], n, K, prior_var.cand)

      if (is.finite(ml.cand$loglik)) {
        vars.dead <- vars[ii - 1L, tokill, seq_len(nint[ii - 1L, tokill])]
        degs.dead <- degs[ii - 1L, tokill, seq_len(nint[ii - 1L, tokill])]
        key.dead <- rkeys.curr[tokill]

        eta.cand <- eta_vec
        eta.cand[vars.dead] <- eta.cand[vars.dead] - 1
        log_q_basis_rev <- .jk_basis_log_q(vars.dead, degs.dead, eta.cand,
                                           degree, order, J_probs, coin_pars, degree_penalty)
        counts.cand <- .jk_subset_counts(rkeys.cand)
        log_q_subset_rev <- .jk_subset_log_q_mix(key.dead, K, counts.cand, M.prev - 1L,
                                                 subset_rho, subset_power_prop)
        lp_subset <- .jk_subset_log_prior(key.dead, K, subset_power_prior)

        llik.alpha <- ml.cand$loglik - ml.curr$loglik
        lprior.alpha <- -log(lam[ii - 1L]) - log_prior_basis - lp_subset
        lprop.alpha <- (log(.jk_move_prob("birth", M.prev - 1L, max_basis, move_probs)) +
                          log_q_basis_rev + log_q_subset_rev) -
          (log(.jk_move_prob("death", M.prev, max_basis, move_probs)) - log(M.prev))
        alpha <- llik.alpha + lprior.alpha + lprop.alpha

        if (log(stats::runif(1)) < alpha) {
          H.curr <- H.cand
          rkeys.curr <- rkeys.cand
          Xstk.curr <- Xstk.cand
          ml.curr <- .jk_log_marginal_ridge(Xstk.curr, yvec, s2[ii - 1L, ], n, K, prior_var.cand, need_posterior = TRUE)
          nbasis[ii] <- M.prev - 1L

          nint[ii, ] <- dtot[ii, ] <- NA_integer_
          vars[ii, , ] <- degs[ii, , ] <- NA_integer_
          rkeys[ii, ] <- NA_character_
          if (nbasis[ii] > 0L) {
            nint[ii, seq_len(nbasis[ii])] <- nint[ii - 1L, keep]
            dtot[ii, seq_len(nbasis[ii])] <- dtot[ii - 1L, keep]
            vars[ii, seq_len(nbasis[ii]), ] <- vars[ii - 1L, keep, , drop = FALSE]
            degs[ii, seq_len(nbasis[ii]), ] <- degs[ii - 1L, keep, , drop = FALSE]
            rkeys[ii, seq_len(nbasis[ii])] <- rkeys[ii - 1L, keep]
          }

          eta_vec <- eta.cand
          count_accept["Death"] <- count_accept["Death"] + 1
        }
      }
      count_propose["Death"] <- count_propose["Death"] + 1

    } else if (move.type == "change") {
      if (M.prev <= 0L) next
      tochange <- sample(seq_len(M.prev), 1L)
      nint.curr <- nint[ii - 1L, tochange]
      vars.curr <- vars[ii - 1L, tochange, seq_len(nint.curr)]
      degs.curr <- degs[ii - 1L, tochange, seq_len(nint.curr)]
      dtot.curr <- sum(degs.curr)

      can_swap_var <- p > nint.curr
      mutate_type <- if (!can_swap_var) "degree" else sample(c("degree", "vars"), 1L)

      if (mutate_type == "degree") {
        d_probs <- seq_along(nint.curr:degree)^(-degree_penalty)
        d_probs <- d_probs / sum(d_probs)
        dtot.cand <- sample(nint.curr:degree, 1L, prob = d_probs)
        degs.cand <- .jk_random_partition(dtot.cand, nint.curr)
        basis.new <- .jk_make_basis(vars.curr, degs.cand, X)

        H.cand <- H.curr
        H.cand[, tochange] <- basis.new
        Xstk.cand <- .jk_build_design(H.cand, rkeys.curr, n, K)
        prior_var.curr <- c(rep(tau2_intercept, K), rep(tau2, M.prev))
        ml.cand <- .jk_log_marginal_ridge(Xstk.cand, yvec, s2[ii - 1L, ], n, K, prior_var.curr)

        if (is.finite(ml.cand$loglik)) {
          log_q_forward <- log(d_probs[dtot.cand - nint.curr + 1L]) - lchoose(dtot.cand - 1L, nint.curr - 1L)
          log_q_reverse <- log(d_probs[dtot.curr - nint.curr + 1L]) - lchoose(dtot.curr - 1L, nint.curr - 1L)
          alpha <- (ml.cand$loglik - ml.curr$loglik) + (log_q_reverse - log_q_forward)

          if (log(stats::runif(1)) < alpha) {
            H.curr <- H.cand
            Xstk.curr <- Xstk.cand
            ml.curr <- .jk_log_marginal_ridge(Xstk.curr, yvec, s2[ii - 1L, ], n, K,
                                              prior_var.curr, need_posterior = TRUE)
            degs[ii, tochange, seq_len(nint.curr)] <- degs.cand
            dtot[ii, tochange] <- dtot.cand
            count_accept["Mutate (degree)"] <- count_accept["Mutate (degree)"] + 1
          }
        }
        count_propose["Mutate (degree)"] <- count_propose["Mutate (degree)"] + 1

      } else {
        ind_tochange <- sample(seq_len(nint.curr), 1L)
        oldvar <- vars.curr[ind_tochange]
        forbidden.forward <- setdiff(vars.curr, oldvar)
        probs.forward <- .jk_newvar_probs(eta_vec, forbidden.forward, coin_pars[[3]])
        if (!is.null(probs.forward)) {
          probs.forward[oldvar] <- 0
          probs.forward <- probs.forward / sum(probs.forward)
          newvar <- sample(seq_len(p), 1L, prob = probs.forward)

          vars.cand <- vars.curr
          vars.cand[ind_tochange] <- newvar
          basis.new <- .jk_make_basis(vars.cand, degs.curr, X)

          eta.cand <- eta_vec
          eta.cand[oldvar] <- eta.cand[oldvar] - 1
          eta.cand[newvar] <- eta.cand[newvar] + 1

          forbidden.reverse <- setdiff(vars.cand, newvar)
          probs.reverse <- .jk_newvar_probs(eta.cand, forbidden.reverse, coin_pars[[3]])
          if (!is.null(probs.reverse)) {
            probs.reverse[newvar] <- 0
            probs.reverse <- probs.reverse / sum(probs.reverse)

            H.cand <- H.curr
            H.cand[, tochange] <- basis.new
            Xstk.cand <- .jk_build_design(H.cand, rkeys.curr, n, K)
            prior_var.curr <- c(rep(tau2_intercept, K), rep(tau2, M.prev))
            ml.cand <- .jk_log_marginal_ridge(Xstk.cand, yvec, s2[ii - 1L, ], n, K, prior_var.curr)

            if (is.finite(ml.cand$loglik)) {
              alpha <- (ml.cand$loglik - ml.curr$loglik) +
                (log(probs.reverse[oldvar]) - log(probs.forward[newvar]))

              if (log(stats::runif(1)) < alpha) {
                H.curr <- H.cand
                Xstk.curr <- Xstk.cand
                ml.curr <- .jk_log_marginal_ridge(Xstk.curr, yvec, s2[ii - 1L, ], n, K,
                                                  prior_var.curr, need_posterior = TRUE)
                vars[ii, tochange, seq_len(nint.curr)] <- vars.cand
                eta_vec <- eta.cand
                count_accept["Mutate (vars)"] <- count_accept["Mutate (vars)"] + 1
              }
            }
          }
        }
        count_propose["Mutate (vars)"] <- count_propose["Mutate (vars)"] + 1
      }
    }

    # Gibbs updates: coefficients, sigma_k^2, lambda.
    M.now <- as.integer(nbasis[ii])
    prior_var.now <- c(rep(tau2_intercept, K), rep(tau2, M.now))
    ml.curr <- .jk_log_marginal_ridge(Xstk.curr, yvec, s2[ii - 1L, ], n, K,
                                      prior_var.now, need_posterior = TRUE)
    beta.curr <- .jk_rmvnorm_precision_chol(ml.curr$mu, ml.curr$R)
    beta[ii, seq_len(K + M.now)] <- beta.curr

    fitted.vec <- as.vector(Xstk.curr %*% beta.curr)
    fitted.mat <- matrix(fitted.vec, nrow = n, ncol = K)
    resid <- Y - fitted.mat
    ss <- colSums(resid^2)
    sum_sq[ii, ] <- ss

    for (kk in seq_len(K)) {
      s2[ii, kk] <- max(s2_lower,
                        1 / stats::rgamma(1, shape = n / 2 + g1, rate = g2 + 0.5 * ss[kk]))
    }
    lam[ii] <- stats::rgamma(1, shape = h1 + M.now, rate = h2 + 1)

    # Refresh collapsed state at the new sigma values for the next RJ proposal.
    ml.curr <- .jk_log_marginal_ridge(Xstk.curr, yvec, s2[ii, ], n, K,
                                      prior_var.now, need_posterior = TRUE)

    if (verbose && ii %% 1000 == 0L) {
      cat(c("MCMC iteration", ii, .jk_timestamp(), "nbasis:", M.now), "\n")
    }
  }

  # Trim storage.
  mcmc_iter <- seq(nburn + 1L, nmcmc, by = thin)
  max_basis_used <- max(nbasis, na.rm = TRUE)
  basis_high <- if (is.finite(max_basis_used) && max_basis_used > 0L) seq_len(max_basis_used) else integer(0)
  max_inter_used <- suppressWarnings(max(nint, na.rm = TRUE))
  inter_high <- if (is.finite(max_inter_used) && max_inter_used > 0L) seq_len(max_inter_used) else integer(0)

  vars.out <- vars[mcmc_iter, basis_high, inter_high, drop = FALSE]
  degs.out <- degs[mcmc_iter, basis_high, inter_high, drop = FALSE]
  nint.out <- nint[mcmc_iter, basis_high, drop = FALSE]
  dtot.out <- dtot[mcmc_iter, basis_high, drop = FALSE]
  rkeys.out <- rkeys[mcmc_iter, basis_high, drop = FALSE]
  beta.out <- beta[mcmc_iter, seq_len(K + max_basis_used), drop = FALSE]

  out <- list(
    H = H.curr,
    Xstack = Xstk.curr,
    vars = vars.out,
    degs = degs.out,
    nint = nint.out,
    dtot = dtot.out,
    rkeys = rkeys.out,
    nbasis = as.integer(nbasis[mcmc_iter]),
    beta = beta.out,
    s2 = s2[mcmc_iter, , drop = FALSE],
    lam = lam[mcmc_iter],
    sum_sq = sum_sq[mcmc_iter, , drop = FALSE],
    eta = eta_vec,
    count_accept = count_accept,
    count_propose = count_propose,
    prior_type = "ridge",
    model_type = "joint_khaos_shared_coefficients",
    X = X,
    Y = Y,
    K = K,
    n = n,
    p = p,
    degree = degree,
    order = order,
    tau2 = tau2,
    tau2_intercept = tau2_intercept,
    subset_rho = subset_rho,
    subset_power_prop = subset_power_prop,
    subset_power_prior = subset_power_prior,
    call = match.call()
  )
  class(out) <- "adaptive_jkhaos"
  out
}

# -----------------------------------------------------------------------------
# Prediction method.
# -----------------------------------------------------------------------------

#' Predict method for adaptive_jkhaos
#'
#' @param object Object returned by adaptive_jkhaos_ridge.
#' @param newdata Optional new n by p design matrix.
#' @param mcmc.use Posterior sample indices to use.
#' @param nugget If TRUE, add independent N(0, sigma_k^2) noise.
#' @param nreps Number of nugget replicates per MCMC draw.
#' @param simplify If TRUE and nreps=1, return an array draws by n_new by K.
#' @param ... Ignored.
#' @export
predict.adaptive_jkhaos <- function(object, newdata = NULL, mcmc.use = NULL,
                                    nugget = FALSE, nreps = 1, simplify = TRUE, ...) {
  if (is.null(newdata)) newdata <- object$X
  if (!is.matrix(newdata)) {
    if (is.data.frame(newdata)) newdata <- as.matrix(newdata) else newdata <- matrix(newdata, ncol = object$p)
  }
  if (is.null(mcmc.use)) mcmc.use <- seq_along(object$nbasis)
  if (!nugget) nreps <- 1L

  nnew <- nrow(newdata)
  K <- object$K
  ndraw <- length(mcmc.use) * nreps
  pred <- array(NA_real_, dim = c(ndraw, nnew, K))

  cnt <- 1L
  for (ii in mcmc.use) {
    M <- object$nbasis[ii]
    Hnew <- if (M > 0L) matrix(NA_real_, nrow = nnew, ncol = M) else matrix(nrow = nnew, ncol = 0L)
    if (M > 0L) {
      for (mm in seq_len(M)) {
        vv <- object$vars[ii, mm, seq_len(object$nint[ii, mm])]
        dd <- object$degs[ii, mm, seq_len(object$nint[ii, mm])]
        Hnew[, mm] <- .jk_make_basis(vv, dd, newdata)
      }
    }

    beta.curr <- object$beta[ii, seq_len(K + M)]
    mu <- matrix(NA_real_, nrow = nnew, ncol = K)
    for (kk in seq_len(K)) {
      mu[, kk] <- beta.curr[kk]
    }
    if (M > 0L) {
      for (mm in seq_len(M)) {
        rr <- .jk_vec_from_key(object$rkeys[ii, mm])
        if (length(rr) > 0L) {
          for (kk in rr) {
            mu[, kk] <- mu[, kk] + beta.curr[K + mm] * Hnew[, mm]
          }
        }
      }
    }

    for (rrr in seq_len(nreps)) {
      if (nugget) {
        noise <- matrix(0, nrow = nnew, ncol = K)
        for (kk in seq_len(K)) noise[, kk] <- stats::rnorm(nnew, 0, sqrt(object$s2[ii, kk]))
        pred[cnt, , ] <- mu + noise
      } else {
        pred[cnt, , ] <- mu
      }
      cnt <- cnt + 1L
    }
  }

  if (simplify) return(pred)
  pred
}

# -----------------------------------------------------------------------------
# Plot method.
# -----------------------------------------------------------------------------

#' Plot method for adaptive_jkhaos
#'
#' @param x Object returned by adaptive_jkhaos_ridge.
#' @param ... Ignored.
#' @export
plot.adaptive_jkhaos <- function(x, ...){

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))

  graphics::par(
    mfrow = c(2,2),
    mar = c(4,4,2,1),
    oma = c(0,0,0,0)
  )

  K <- ncol(x$Y)

  cols <- grDevices::rainbow(K)
  pchs <- seq_len(K)

  ##----------------------------------
  ## Predictions
  ##----------------------------------

  preds <- stats::predict(x, nugget = TRUE, nreps = 10)

  ## Assume predict returns:
  ## draws x n x K
  yhat <- apply(preds, c(2,3), mean)
  ci   <- 2 * apply(preds, c(2,3), stats::sd)

  ##----------------------------------
  ## Panel 1: Basis traces
  ##----------------------------------

  stats::ts.plot(
    x$nbasis,
    ylab = "nbasis",
    col = "black",
    lwd = 2
  )

  ## x$nbasis_output should be:
  ## iterations x K
  if(!is.null(x$nbasis_output)){

    for(k in seq_len(K)){
      graphics::lines(
        x$nbasis_output[,k],
        col = cols[k],
        lwd = 1
      )
    }

    graphics::legend(
      "topleft",
      legend = c("Total", paste0("Y",1:K)),
      col = c("black", cols),
      lwd = c(2, rep(1,K)),
      bty = "n",
      cex = 0.8
    )
  }

  ##----------------------------------
  ## Panel 2: sigma traces
  ##----------------------------------

  matplot(
    x$s2,
    type = "l",
    lty = 1,
    col = cols,
    xlab = "Iteration",
    ylab = expression(sigma^2)
  )

  graphics::legend(
    "topright",
    legend = paste0("Y",1:K),
    col = cols,
    lty = 1,
    bty = "n",
    cex = 0.8
  )

  ##----------------------------------
  ## Panel 3: Observed vs fitted
  ##----------------------------------

  rng <- range(c(x$y, yhat))

  graphics::plot(
    rng,
    rng,
    type = "n",
    xlab = "Observed",
    ylab = "Predicted"
  )

  graphics::abline(0,1,col="dodgerblue")

  for(k in seq_len(K)){

    graphics::segments(
      x0 = x$y[,k],
      y0 = yhat[,k]-ci[,k],
      y1 = yhat[,k]+ci[,k],
      col = adjustcolor(cols[k],0.4)
    )

    graphics::points(
      x$y[,k],
      yhat[,k],
      pch = pchs[k],
      col = cols[k]
    )
  }

  graphics::legend(
    "topleft",
    legend = paste0("Y",1:K),
    pch = pchs,
    col = cols,
    bty = "n",
    cex = 0.8
  )

  ##----------------------------------
  ## Panel 4: Residual KDEs
  ##----------------------------------

  residuals <- x$y - yhat

  rr <- as.vector(residuals)

  d <- stats::density(rr)

  graphics::plot(
    d,
    type = "n",
    xlab = "Residual",
    ylab = "Density"
  )

  for(k in seq_len(K)){

    dk <- stats::density(residuals[,k])

    graphics::lines(
      dk,
      col = cols[k],
      lwd = 2
    )
  }

  xx <- seq(min(rr), max(rr), length.out = 200)

  graphics::lines(
    xx,
    stats::dnorm(xx, mean(rr), stats::sd(rr)),
    col = "black",
    lwd = 2,
    lty = 2
  )

  graphics::legend(
    "topright",
    legend = c(paste0("Y",1:K),"Normal"),
    col = c(cols,"black"),
    lty = c(rep(1,K),2),
    lwd = 2,
    bty = "n",
    cex = 0.8
  )

  invisible(NULL)
}
