## ---------- unified signature --------------------------------
log_evidence_fast_core <- function(yty, v, g_diag, g0_sq,
                                   nu0, s0_sq, n)
{                 # <– fast orthogonal approximation
  logdet  <- 0.5 * sum(log(g0_sq * g_diag / (g0_sq * g_diag + 1)))
  quad    <- (yty - sum((g0_sq * g_diag * v^2) /
                          (n * (g0_sq * g_diag + 1)))) / 2
  out <- lgamma((nu0 + n)/2) - lgamma(nu0/2) -
    (nu0 + n)/2 * log(nu0*s0_sq/2 + quad) +
    logdet - n/2 * log(2*pi)
  out
}

log_evidence_full_core <- function(BtB, BtBi, Bty,
                                   g_diag, g0_sq,
                                   nu0, s0_sq, yty, n)
{                 # <– full non-orthogonal
  G      <- outer(rep(1, length(g_diag)), rep(1, length(g_diag))) +
    (1 / g0_sq) * tcrossprod(1 / sqrt(g_diag))
  Sigma_i <- G * BtB
  L       <- chol(Sigma_i)
  logdet  <- -2 * sum(log(diag(L)))

  a_hat   <- BtBi %*% Bty
  resid2  <- yty - crossprod(a_hat, Bty)

  quad <- (resid2 + nu0 * s0_sq) / 2
  out  <- lgamma((nu0 + n)/2) - lgamma(nu0/2) -
    (nu0 + n)/2 * log(quad) +
    0.5 * logdet - n/2 * log(2*pi)
  as.numeric(out)
}

## ---------- public wrapper -----------------------------------
log_evidence <- function(mode = c("fast","full"),
                         yty, v = NULL,              # fast needs v
                         BtB = NULL, BtBi = NULL, Bty = NULL,  # full needs these
                         g_diag, g0_sq,
                         nu0, s0_sq, n)
{
  mode <- match.arg(mode)
  if (mode == "fast") {
    log_evidence_fast_core(yty, v, g_diag, g0_sq,
                           nu0, s0_sq, n)
  } else {
    log_evidence_full_core(BtB, BtBi, Bty,
                           g_diag, g0_sq,
                           nu0, s0_sq, yty, n)
  }
}
