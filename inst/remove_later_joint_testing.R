library(BASS)
library(khaos)
library(coactivity)

X <- lhs::randomLHS(1000, 8)
y <- apply(X, 1, duqling::borehole)

fit_bass <- bass(X, y)
fit_khaos <- adaptive_khaos(X, y, prior_type = "ridge",
                            nmcmc=30000, nburn=28000, thin=2)

Cmat_bass <- C_bass(fit_bass)
Cmat_khaos <- C_khaos(fit_khaos)
Cmat_khaos_s <- Cmat_khaos[[1]] / 1000
for(i in 2:1000){
  Cmat_khaos_s <- Cmat_khaos_s + Cmat_khaos[[i]] / 1000
}

eb <- eigen(Cmat_bass)
ek <- eigen(Cmat_khaos_s)

par(mfrow=c(2,2))
plot(sqrt(eb$val))
points(sqrt(ek$val), col='blue')

plot(eb$vec[,1])
points(ek$vec[,1], col='blue', pch=2)

plot(eb$vec[,2])
points(-ek$vec[,2], col='blue', pch=2)

plot(eb$vec[,3])
points(-ek$vec[,3], col='blue', pch=2)


# A quick K=3 example
X <- lhs::randomLHS(1000, 3)
y1 <- apply(X,1, duqling::ishigami) / 3.8
y2 <- y1 + apply(X, 1, duqling::dms_additive) / 5
y3 <- y1 + apply(X, 1, duqling::dms_additive) / 5 + apply(X, 1, duqling::dms_complicated) / 10
y4 <- y2 - y1
Y <- cbind(y1, y2, y3, y4)

fit <- adaptive_jkhaos_ridge(X, Y)
plot_sharing_pyramid(fit, metric="coef_energy")
plot_sharing_pyramid(fit, metric="n_terms")
#plot_sharing_pyramid(fit, metric="trace_C")

class(fit) <- c("adaptive_khaos", "adaptive_jkhaos")
Cfull <- C_khaos(fit, mcmc.use=1000)

plot_shared_residuals(fit)
plot_shared_residuals_grid(fit)


plot_sharing_pyramid <- function(
    obj,
    mcmc.use = NULL,
    metric = c("n_terms", "coef_energy", "trace_C"),
    summarize = c("mean", "median"),
    min_strength = 0,
    stagger = 0.08,
    point_range = c(2, 9),
    line_width = 0.8,
    output_labels = NULL
) {
  metric <- match.arg(metric)
  summarize <- match.arg(summarize)

  if (!inherits(obj, "adaptive_jkhaos")) {
    stop("obj must be an object of class 'adaptive_jkhaos'.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  K <- obj$K
  if (is.null(K)) K <- ncol(obj$Y)

  if (is.null(output_labels)) {
    output_labels <- colnames(obj$Y)
    if (is.null(output_labels)) output_labels <- paste0("f", seq_len(K))
  }

  if (is.null(mcmc.use)) {
    mcmc.use <- seq_along(obj$nbasis)
  }

  subset_values <- list()

  # For trace_C metric, we need to precompute C matrices
  if (metric == "trace_C") {
    # Get C matrices for all mcmc.use samples
    C_list <- C_khaos(obj, mcmc.use = mcmc.use, verbose = FALSE)

    # Precompute integral tables once
    max_degree <- max(obj$degs, na.rm = TRUE)
    if (!is.finite(max_degree)) max_degree <- 0
    legendre_tabs <- make_khaos_integral_tables(max_degree)
    A <- legendre_tabs$A
    B <- legendre_tabs$B
    p <- ncol(obj$X)
  }

  for (ii in seq_along(mcmc.use)) {
    s <- mcmc.use[ii]
    M <- obj$nbasis[s]

    if (is.na(M) || M < 1) next

    # Build a map from subset keys to basis indices for this iteration
    key_to_bases <- list()
    for (m in seq_len(M)) {
      key <- obj$rkeys[s, m]
      if (is.na(key) || !nzchar(key)) next

      if (is.null(key_to_bases[[key]])) {
        key_to_bases[[key]] <- integer(0)
      }
      key_to_bases[[key]] <- c(key_to_bases[[key]], m)
    }

    # Compute values for each unique subset
    for (key in names(key_to_bases)) {
      bases_in_subset <- key_to_bases[[key]]

      val <- switch(
        metric,
        n_terms = length(bases_in_subset),
        coef_energy = {
          # Sum of squared coefficients for this subset
          sum(obj$beta[s, K + bases_in_subset]^2, na.rm = TRUE)
        },
        trace_C = {
          # Compute C matrix trace for this subset
          compute_subset_C_trace(
            obj = obj,
            s = s,
            basis_indices = bases_in_subset,
            p = p,
            A = A,
            B = B
          )
        }
      )

      if (is.null(subset_values[[key]])) {
        subset_values[[key]] <- numeric(length(mcmc.use))
      }

      subset_values[[key]][ii] <- val
    }
  }

  if (length(subset_values) == 0) {
    stop("No active subset structure found in obj$rkeys.")
  }

  reduce_fun <- switch(
    summarize,
    mean = mean,
    median = stats::median
  )

  subset_df <- do.call(
    rbind,
    lapply(names(subset_values), function(key) {
      members <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
      strength <- reduce_fun(subset_values[[key]], na.rm = TRUE)

      data.frame(
        subset_key = key,
        subset_label = paste(members, collapse = ":"),
        subset_size = length(members),
        xmin = min(members),
        xmax = max(members),
        strength = strength,
        stringsAsFactors = FALSE
      )
    })
  )

  subset_df <- subset_df[subset_df$strength > min_strength, , drop = FALSE]

  if (nrow(subset_df) == 0) {
    stop("No subsets remain after applying min_strength.")
  }

  subset_df$y_base <- K - subset_df$subset_size + 1
  subset_df$y <- NA_real_

  for (q in sort(unique(subset_df$subset_size))) {
    idx <- which(subset_df$subset_size == q)

    idx <- idx[
      order(
        subset_df$xmin[idx],
        subset_df$xmax[idx],
        subset_df$subset_label[idx]
      )
    ]

    offsets <- seq_along(idx) - mean(seq_along(idx))
    subset_df$y[idx] <- subset_df$y_base[idx] + stagger * offsets
  }

  point_df <- do.call(
    rbind,
    lapply(seq_len(nrow(subset_df)), function(i) {
      members <- as.integer(strsplit(subset_df$subset_key[i], ",", fixed = TRUE)[[1]])

      data.frame(
        subset_key = subset_df$subset_key[i],
        subset_label = subset_df$subset_label[i],
        output = members,
        y = subset_df$y[i],
        strength = subset_df$strength[i],
        subset_size = subset_df$subset_size[i],
        stringsAsFactors = FALSE
      )
    })
  )

  y_breaks <- K:1
  y_labels <- paste0("|r| = ", seq_len(K))

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = subset_df,
      ggplot2::aes(
        x = xmin,
        xend = xmax,
        y = y,
        yend = y,
        color = strength
      ),
      linewidth = line_width,
      alpha = 0.85
    ) +
    ggplot2::geom_point(
      data = point_df,
      ggplot2::aes(
        x = output,
        y = y,
        color = strength,
        size = strength
      ),
      alpha = 0.95
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(K),
      labels = output_labels,
      limits = c(0.5, K + 0.5)
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      limits = c(0.5, K + 0.5)
    ) +
    ggplot2::scale_size_continuous(range = point_range) +
    ggplot2::labs(
      x = "Output",
      y = "Sharing layer",
      color = metric,
      size = metric
    ) +
    ggplot2::theme_minimal()
}

# Helper function to compute C matrix trace for a subset of bases
compute_subset_C_trace <- function(obj, s, basis_indices, p, A, B) {
  M_subset <- length(basis_indices)

  if (M_subset == 0) return(0)

  # Build alpha matrix for this subset only
  alpha <- matrix(0L, nrow = M_subset, ncol = p)

  for (idx in seq_along(basis_indices)) {
    m <- basis_indices[idx]
    q <- obj$nint[s, m]

    if (is.na(q) || q == 0) next

    vars_m <- obj$vars[s, m, seq_len(q)]
    degs_m <- obj$degs[s, m, seq_len(q)]

    keep <- !is.na(vars_m) & !is.na(degs_m)

    vars_m <- vars_m[keep]
    degs_m <- degs_m[keep]

    alpha[idx, vars_m] <- degs_m
  }

  # Get beta values for this subset
  K <- obj$K
  beta <- obj$beta[s, K + basis_indices]

  # Compute C matrix (diagonal only for trace)
  C_trace <- 0

  for (i in seq_len(p)) {
    for (m1 in seq_len(M_subset)) {
      for (m2 in seq_len(M_subset)) {

        if (p > 1 && any(alpha[m1, -i] != alpha[m2, -i])) next

        a <- alpha[m1, i]
        b <- alpha[m2, i]

        C_trace <- C_trace + beta[m1] * beta[m2] * B[a + 1, b + 1]
      }
    }
  }

  return(as.numeric(C_trace))
}

plot_sharing_pyramid <- function(
    obj,
    mcmc.use = NULL,
    metric = c("n_terms", "coef_energy"),
    summarize = c("mean", "median"),
    min_strength = 0,
    stagger = 0.08,
    point_range = c(2, 9),
    line_width = 0.8,
    output_labels = NULL
) {
  metric <- match.arg(metric)
  summarize <- match.arg(summarize)

  if (!inherits(obj, "adaptive_jkhaos")) {
    stop("obj must be an object of class 'adaptive_jkhaos'.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  K <- obj$K
  if (is.null(K)) K <- ncol(obj$Y)

  if (is.null(output_labels)) {
    output_labels <- colnames(obj$Y)
    if (is.null(output_labels)) output_labels <- paste0("f", seq_len(K))
  }

  if (is.null(mcmc.use)) {
    mcmc.use <- seq_along(obj$nbasis)
  }

  subset_values <- list()

  for (ii in seq_along(mcmc.use)) {
    s <- mcmc.use[ii]
    M <- obj$nbasis[s]

    if (is.na(M) || M < 1) next

    for (m in seq_len(M)) {
      key <- obj$rkeys[s, m]

      if (is.na(key) || !nzchar(key)) next

      val <- switch(
        metric,
        n_terms = 1,
        coef_energy = obj$beta[s, K + m]^2
      )

      if (is.null(subset_values[[key]])) {
        subset_values[[key]] <- numeric(length(mcmc.use))
      }

      subset_values[[key]][ii] <- subset_values[[key]][ii] + val
    }
  }

  if (length(subset_values) == 0) {
    stop("No active subset structure found in obj$rkeys.")
  }

  reduce_fun <- switch(
    summarize,
    mean = mean,
    median = stats::median
  )

  subset_df <- do.call(
    rbind,
    lapply(names(subset_values), function(key) {
      members <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
      strength <- reduce_fun(subset_values[[key]], na.rm = TRUE)

      data.frame(
        subset_key = key,
        subset_label = paste(members, collapse = ":"),
        subset_size = length(members),
        xmin = min(members),
        xmax = max(members),
        strength = strength,
        stringsAsFactors = FALSE
      )
    })
  )

  subset_df <- subset_df[subset_df$strength > min_strength, , drop = FALSE]

  if (nrow(subset_df) == 0) {
    stop("No subsets remain after applying min_strength.")
  }

  subset_df$y_base <- K - subset_df$subset_size + 1
  subset_df$y <- NA_real_

  for (q in sort(unique(subset_df$subset_size))) {
    idx <- which(subset_df$subset_size == q)

    idx <- idx[
      order(
        subset_df$xmin[idx],
        subset_df$xmax[idx],
        subset_df$subset_label[idx]
      )
    ]

    offsets <- seq_along(idx) - mean(seq_along(idx))
    subset_df$y[idx] <- subset_df$y_base[idx] + stagger * offsets
  }

  point_df <- do.call(
    rbind,
    lapply(seq_len(nrow(subset_df)), function(i) {
      members <- as.integer(strsplit(subset_df$subset_key[i], ",", fixed = TRUE)[[1]])

      data.frame(
        subset_key = subset_df$subset_key[i],
        subset_label = subset_df$subset_label[i],
        output = members,
        y = subset_df$y[i],
        strength = subset_df$strength[i],
        subset_size = subset_df$subset_size[i],
        stringsAsFactors = FALSE
      )
    })
  )

  y_breaks <- K:1
  y_labels <- paste0("|r| = ", seq_len(K))

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = subset_df,
      ggplot2::aes(
        x = xmin,
        xend = xmax,
        y = y,
        yend = y,
        color = strength
      ),
      linewidth = line_width,
      alpha = 0.85
    ) +
    ggplot2::geom_point(
      data = point_df,
      ggplot2::aes(
        x = output,
        y = y,
        color = strength,
        size = strength
      ),
      alpha = 0.95
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(K),
      labels = output_labels,
      limits = c(0.5, K + 0.5)
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      limits = c(0.5, K + 0.5)
    ) +
    ggplot2::scale_size_continuous(range = point_range) +
    ggplot2::labs(
      x = "Output",
      y = "Sharing layer",
      color = metric,
      size = metric
    ) +
    ggplot2::theme_minimal()
}


extract_jkhaos_component <- function(
    obj,
    subset,
    exact = TRUE,
    keep_intercepts = FALSE
) {
  if (!inherits(obj, "adaptive_jkhaos")) {
    stop("obj must be an object of class 'adaptive_jkhaos'.")
  }

  K <- obj$K
  if (is.null(K)) K <- ncol(obj$Y)

  subset <- sort(unique(as.integer(subset)))
  if (any(is.na(subset)) || any(subset < 1) || any(subset > K)) {
    stop("subset must contain valid output indices.")
  }

  subset_key <- paste(subset, collapse = ",")

  out <- obj

  n_iter <- length(obj$nbasis)
  max_old <- dim(obj$vars)[2]

  keep_list <- vector("list", n_iter)
  nbasis_new <- integer(n_iter)
  max_new <- 0L

  for (s in seq_len(n_iter)) {
    M <- obj$nbasis[s]

    if (is.na(M) || M == 0) {
      keep_list[[s]] <- integer(0)
      next
    }

    keys_s <- obj$rkeys[s, seq_len(M)]

    if (exact) {
      keep <- which(keys_s == subset_key)
    } else {
      keep <- which(vapply(keys_s, function(key) {
        if (is.na(key) || !nzchar(key)) return(FALSE)
        members <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
        all(subset %in% members)
      }, logical(1)))
    }

    keep_list[[s]] <- keep
    nbasis_new[s] <- length(keep)
    max_new <- max(max_new, nbasis_new[s])
  }

  if (max_new == 0) {
    warning("No matching basis functions found.")
    max_new <- 1L
  }

  out$nbasis <- nbasis_new

  out$vars <- array(
    NA_integer_,
    dim = c(n_iter, max_new, dim(obj$vars)[3])
  )

  out$degs <- array(
    NA_integer_,
    dim = c(n_iter, max_new, dim(obj$degs)[3])
  )

  out$nint <- matrix(
    NA_integer_,
    nrow = n_iter,
    ncol = max_new
  )

  out$rkeys <- matrix(
    NA_character_,
    nrow = n_iter,
    ncol = max_new
  )

  out$beta <- matrix(
    0,
    nrow = n_iter,
    ncol = K + max_new
  )

  if (keep_intercepts) {
    out$beta[, seq_len(K)] <- obj$beta[, seq_len(K), drop = FALSE]
  }

  for (s in seq_len(n_iter)) {
    keep <- keep_list[[s]]
    M_new <- length(keep)

    if (M_new == 0) next

    out$vars[s, seq_len(M_new), ] <- obj$vars[s, keep, , drop = FALSE]
    out$degs[s, seq_len(M_new), ] <- obj$degs[s, keep, , drop = FALSE]
    out$nint[s, seq_len(M_new)] <- obj$nint[s, keep]
    out$rkeys[s, seq_len(M_new)] <- obj$rkeys[s, keep]

    out$beta[s, K + seq_len(M_new)] <- obj$beta[s, K + keep]
  }

  attr(out, "component_subset") <- subset
  attr(out, "component_subset_key") <- subset_key
  attr(out, "component_exact") <- exact

  class(out) <- unique(c("jkhaos_component", class(obj)))

  out
}


#' Visualize Residuals from Shared vs Full Model Components with FVU
#'
#' Creates a visualization comparing FVU (Fraction of Variance Unexplained)
#' when using only shared bases (|r| > 1) versus using all bases versus no model.
#'
#' @param obj An object of class "adaptive_jkhaos".
#' @param mcmc.use Vector of MCMC sample indices to use for prediction.
#' @param output_labels Optional vector of labels for outputs.
#' @param colors Optional vector of colors for each output.
#' @param show_legend Logical. Should a legend be displayed?
#'
#' @return Invisibly returns a list with residual statistics.
#' @export
plot_shared_residuals <- function(obj,
                                  mcmc.use = NULL,
                                  output_labels = NULL,
                                  colors = NULL,
                                  show_legend = TRUE) {

  if (!inherits(obj, "adaptive_jkhaos")) {
    stop("obj must be an object of class 'adaptive_jkhaos'.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  K <- obj$K
  if (is.null(K)) K <- ncol(obj$Y)

  if (is.null(output_labels)) {
    output_labels <- colnames(obj$Y)
    if (is.null(output_labels)) output_labels <- paste0("f", seq_len(K))
  }

  if (is.null(colors)) {
    colors <- grDevices::rainbow(K)
  }

  if (is.null(mcmc.use)) {
    mcmc.use <- seq_along(obj$nbasis)
  }

  # Helper function to predict with subset filter
  predict_filtered <- function(obj, subset_filter, mcmc.use) {
    nnew <- nrow(obj$X)
    K <- obj$K
    ndraw <- length(mcmc.use)
    pred <- array(0, dim = c(ndraw, nnew, K))

    for (idx in seq_along(mcmc.use)) {
      ii <- mcmc.use[idx]
      M <- obj$nbasis[ii]

      # Start with intercepts
      for (kk in seq_len(K)) {
        pred[idx, , kk] <- obj$beta[ii, kk]
      }

      if (M > 0) {
        # Add basis contributions
        for (mm in seq_len(M)) {
          key <- obj$rkeys[ii, mm]
          if (is.na(key) || !nzchar(key)) next

          # Check if this basis passes the filter
          r <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
          if (!subset_filter(r)) next

          # Build basis function
          vv <- obj$vars[ii, mm, seq_len(obj$nint[ii, mm])]
          dd <- obj$degs[ii, mm, seq_len(obj$nint[ii, mm])]
          basis_vals <- .jk_make_basis(vv, dd, obj$X)

          # Add to relevant outputs
          for (kk in r) {
            pred[idx, , kk] <- pred[idx, , kk] + obj$beta[ii, K + mm] * basis_vals
          }
        }
      }
    }

    apply(pred, c(2, 3), mean)
  }

  # Get predictions from shared bases only (|r| > 1)
  pred_shared <- predict_filtered(obj, function(r) length(r) > 1, mcmc.use)

  # Get predictions from all bases
  pred_all <- predict_filtered(obj, function(r) TRUE, mcmc.use)

  # Baseline: just use mean
  pred_baseline <- matrix(colMeans(obj$Y), nrow = nrow(obj$Y), ncol = K, byrow = TRUE)

  # Compute residuals
  resid_baseline <- obj$Y - pred_baseline
  resid_shared <- obj$Y - pred_shared
  resid_all <- obj$Y - pred_all

  # Compute total variance for each output
  total_var <- apply(obj$Y, 2, var)

  # Compute FVU (Fraction of Variance Unexplained)
  fvu_baseline <- colMeans(resid_baseline^2) / total_var
  fvu_shared <- colMeans(resid_shared^2) / total_var
  fvu_all <- colMeans(resid_all^2) / total_var

  # R² = 1 - FVU
  r2_baseline <- 1 - fvu_baseline  # Should be ~0
  r2_shared <- 1 - fvu_shared
  r2_all <- 1 - fvu_all

  # Prepare data for plotting
  n <- nrow(obj$Y)
  plot_data <- data.frame(
    output = rep(rep(output_labels, each = n), 3),
    output_num = rep(rep(seq_len(K), each = n), 3),
    residual = c(as.vector(resid_baseline), as.vector(resid_shared), as.vector(resid_all)),
    type = rep(c("Baseline (mean only)", "Shared bases only", "All bases"), each = n * K),
    stringsAsFactors = FALSE
  )

  plot_data$type <- factor(plot_data$type,
                           levels = c("Baseline (mean only)", "Shared bases only", "All bases"))

  stats_data <- data.frame(
    output = rep(output_labels, 3),
    output_num = rep(seq_len(K), 3),
    type = rep(c("Baseline (mean only)", "Shared bases only", "All bases"), each = K),
    fvu = c(fvu_baseline, fvu_shared, fvu_all),
    r2 = c(r2_baseline, r2_shared, r2_all),
    pct_var_explained = c(r2_baseline, r2_shared, r2_all) * 100,
    stringsAsFactors = FALSE
  )
  stats_data$type <- factor(stats_data$type,
                            levels = c("Baseline (mean only)", "Shared bases only", "All bases"))

  # Create the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = output, y = residual, fill = type)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_violin(alpha = 0.6, position = ggplot2::position_dodge(width = 0.9)) +
    ggplot2::geom_boxplot(width = 0.15, alpha = 0.8,
                          position = ggplot2::position_dodge(width = 0.9),
                          outlier.alpha = 0.3) +
    ggplot2::scale_fill_manual(values = c("Baseline (mean only)" = "gray70",
                                          "Shared bases only" = "coral",
                                          "All bases" = "steelblue")) +
    ggplot2::labs(
      x = "Output",
      y = "Residual",
      fill = "Model",
      title = "Residual Comparison: Baseline vs Shared vs Full Model",
      subtitle = "Lower residuals indicate better model fit"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = if (show_legend) "top" else "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  # Add FVU annotations
  fvu_labels <- data.frame(
    output = rep(output_labels, 3),
    type = rep(c("Baseline (mean only)", "Shared bases only", "All bases"), each = K),
    label = sprintf("FVU: %.3f", c(fvu_baseline, fvu_shared, fvu_all)),
    y = rep(max(abs(plot_data$residual)) * 0.95, K * 3),
    stringsAsFactors = FALSE
  )
  fvu_labels$type <- factor(fvu_labels$type,
                            levels = c("Baseline (mean only)", "Shared bases only", "All bases"))

  p <- p + ggplot2::geom_text(
    data = fvu_labels,
    ggplot2::aes(label = label, y = y),
    position = ggplot2::position_dodge(width = 0.9),
    size = 2.5,
    vjust = 1
  )

  print(p)

  # Print summary table
  cat("\n=== Variance Explained Summary ===\n")
  for (kk in seq_len(K)) {
    cat(sprintf("\n%s:\n", output_labels[kk]))
    cat(sprintf("  Baseline       - FVU: %.4f (%.1f%% var explained)\n",
                fvu_baseline[kk], r2_baseline[kk] * 100))
    cat(sprintf("  Shared only    - FVU: %.4f (%.1f%% var explained)\n",
                fvu_shared[kk], r2_shared[kk] * 100))
    cat(sprintf("  All bases      - FVU: %.4f (%.1f%% var explained)\n",
                fvu_all[kk], r2_all[kk] * 100))

    # Show how much shared structure explains
    if (fvu_baseline[kk] > 0) {
      shared_contrib <- (fvu_baseline[kk] - fvu_shared[kk]) / fvu_baseline[kk] * 100
      total_contrib <- (fvu_baseline[kk] - fvu_all[kk]) / fvu_baseline[kk] * 100
      specific_contrib <- (fvu_shared[kk] - fvu_all[kk]) / fvu_baseline[kk] * 100

      cat(sprintf("  → Shared structure explains %.1f%% of total variance\n", shared_contrib))
      cat(sprintf("  → Output-specific structure explains %.1f%% of total variance\n", specific_contrib))
      cat(sprintf("  → Total model explains %.1f%% of total variance\n", total_contrib))
    }
  }

  # Return statistics invisibly
  invisible(list(
    residuals_baseline = resid_baseline,
    residuals_shared = resid_shared,
    residuals_all = resid_all,
    fvu_baseline = fvu_baseline,
    fvu_shared = fvu_shared,
    fvu_all = fvu_all,
    r2_baseline = r2_baseline,
    r2_shared = r2_shared,
    r2_all = r2_all,
    stats_data = stats_data
  ))
}


#' Create a Multi-Panel Residual Comparison Plot with FVU
#'
#' Creates a comprehensive visualization with multiple panels showing:
#' 1. Observed vs predicted for shared-only model
#' 2. Observed vs predicted for full model
#' 3. Residual distributions
#' 4. FVU comparison bar plot (with variance decomposition)
#'
#' @param obj An object of class "adaptive_jkhaos".
#' @param mcmc.use Vector of MCMC sample indices to use for prediction.
#' @param output_labels Optional vector of labels for outputs.
#'
#' @export
plot_shared_residuals_grid <- function(obj,
                                       mcmc.use = NULL,
                                       output_labels = NULL) {

  if (!inherits(obj, "adaptive_jkhaos")) {
    stop("obj must be an object of class 'adaptive_jkhaos'.")
  }

  K <- obj$K
  if (is.null(K)) K <- ncol(obj$Y)

  if (is.null(output_labels)) {
    output_labels <- colnames(obj$Y)
    if (is.null(output_labels)) output_labels <- paste0("f", seq_len(K))
  }

  if (is.null(mcmc.use)) {
    mcmc.use <- seq_along(obj$nbasis)
  }

  # Set up the plotting area
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))

  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  colors <- grDevices::rainbow(K)

  # Get predictions (reuse logic from previous function)
  predict_filtered <- function(obj, subset_filter, mcmc.use) {
    nnew <- nrow(obj$X)
    K <- obj$K
    ndraw <- length(mcmc.use)
    pred <- array(0, dim = c(ndraw, nnew, K))

    for (idx in seq_along(mcmc.use)) {
      ii <- mcmc.use[idx]
      M <- obj$nbasis[ii]

      for (kk in seq_len(K)) {
        pred[idx, , kk] <- obj$beta[ii, kk]
      }

      if (M > 0) {
        for (mm in seq_len(M)) {
          key <- obj$rkeys[ii, mm]
          if (is.na(key) || !nzchar(key)) next

          r <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
          if (!subset_filter(r)) next

          vv <- obj$vars[ii, mm, seq_len(obj$nint[ii, mm])]
          dd <- obj$degs[ii, mm, seq_len(obj$nint[ii, mm])]
          basis_vals <- .jk_make_basis(vv, dd, obj$X)

          for (kk in r) {
            pred[idx, , kk] <- pred[idx, , kk] + obj$beta[ii, K + mm] * basis_vals
          }
        }
      }
    }

    apply(pred, c(2, 3), mean)
  }

  pred_shared <- predict_filtered(obj, function(r) length(r) > 1, mcmc.use)
  pred_all <- predict_filtered(obj, function(r) TRUE, mcmc.use)
  pred_baseline <- matrix(colMeans(obj$Y), nrow = nrow(obj$Y), ncol = K, byrow = TRUE)

  resid_baseline <- obj$Y - pred_baseline
  resid_shared <- obj$Y - pred_shared
  resid_all <- obj$Y - pred_all

  total_var <- apply(obj$Y, 2, var)
  fvu_baseline <- colMeans(resid_baseline^2) / total_var
  fvu_shared <- colMeans(resid_shared^2) / total_var
  fvu_all <- colMeans(resid_all^2) / total_var

  # Panel 1: Observed vs Predicted (Shared only)
  plot(range(obj$Y), range(obj$Y), type = "n",
       xlab = "Observed", ylab = "Predicted",
       main = "Shared Bases Only")
  graphics::abline(0, 1, col = "gray", lty = 2, lwd = 2)
  for (kk in seq_len(K)) {
    graphics::points(obj$Y[, kk], pred_shared[, kk],
                     col = colors[kk], pch = 16, cex = 0.8)
  }
  if (K <= 6) {
    graphics::legend("topleft", legend = output_labels,
                     col = colors, pch = 16, bty = "n", cex = 0.8)
  }
  # Add R² to title
  graphics::mtext(sprintf("Mean R² = %.3f", mean(1 - fvu_shared)),
                  side = 3, line = 0.5, cex = 0.8)

  # Panel 2: Observed vs Predicted (All bases)
  plot(range(obj$Y), range(obj$Y), type = "n",
       xlab = "Observed", ylab = "Predicted",
       main = "All Bases")
  graphics::abline(0, 1, col = "gray", lty = 2, lwd = 2)
  for (kk in seq_len(K)) {
    graphics::points(obj$Y[, kk], pred_all[, kk],
                     col = colors[kk], pch = 16, cex = 0.8)
  }
  graphics::mtext(sprintf("Mean R² = %.3f", mean(1 - fvu_all)),
                  side = 3, line = 0.5, cex = 0.8)

  # Panel 3: Residual distributions
  all_resids <- c(resid_baseline, resid_shared, resid_all)
  plot(c(0, K + 1), range(all_resids), type = "n",
       ylab = "Residual", xlab = "",
       main = "Residual Distributions", yaxt = "n", xaxt="n")
  graphics::axis(1, at = seq_len(K), labels = output_labels, las = 2)
  graphics::axis(2, at=-5:5, labels = -5:5)
  graphics::abline(v = 0, col = "gray", lty = 2)

  for (kk in seq_len(K)) {
    # Baseline
    vioplot_simple(resid_baseline[, kk], at = kk - 0.3,
                   col = adjustcolor("gray70", 0.5), width = 0.2)
    # Shared only
    vioplot_simple(resid_shared[, kk], at = kk,
                   col = adjustcolor(colors[kk], 0.5), width = 0.2)
    # All bases
    vioplot_simple(resid_all[, kk], at = kk + 0.3,
                   col = adjustcolor(colors[kk], 0.8), width = 0.2)
  }

  # Add mini legend
  graphics::legend("topright",
                   legend = c("Baseline", "Shared", "All"),
                   fill = c("gray70", adjustcolor(colors[1], 0.5), adjustcolor(colors[1], 0.8)),
                   bty = "n", cex = 0.7)

  # Panel 4: Stacked bar chart showing variance decomposition
  # For each output: [Shared] [Output-specific] [Unexplained]
  shared_pct <- pmax(0, (fvu_baseline - fvu_shared) / fvu_baseline * 100)
  specific_pct <- pmax(0, (fvu_shared - fvu_all) / fvu_baseline * 100)
  unexplained_pct <- pmax(0, fvu_all / fvu_baseline * 100)

  barplot_data <- rbind(shared_pct, specific_pct, unexplained_pct)
  colnames(barplot_data) <- output_labels

  bp <- graphics::barplot(barplot_data,
                          col = c("coral", "lightblue", "gray80"),
                          ylab = "Percentage of Total Variance",
                          main = "Variance Decomposition",
                          legend.text = c("Shared structure",
                                          "Output-specific structure",
                                          "Unexplained"),
                          args.legend = list(x = "topright", bty = "n", cex = 0.7),
                          ylim = c(0, 100))

  # Add total R² labels
  for (kk in seq_len(K)) {
    r2_total <- 1 - fvu_all[kk]
    graphics::text(bp[kk], 102,
                   sprintf("R²=%.2f", r2_total),
                   cex = 0.7, font = 2)
  }

  graphics::mtext("Residual Analysis: Shared vs Full Model",
                  outer = TRUE, cex = 1.2, font = 2)

  invisible(list(
    pred_baseline = pred_baseline,
    pred_shared = pred_shared,
    pred_all = pred_all,
    resid_baseline = resid_baseline,
    resid_shared = resid_shared,
    resid_all = resid_all,
    fvu_baseline = fvu_baseline,
    fvu_shared = fvu_shared,
    fvu_all = fvu_all,
    r2_baseline = 1 - fvu_baseline,
    r2_shared = 1 - fvu_shared,
    r2_all = 1 - fvu_all
  ))
}

# Simple violin plot helper (unchanged)
vioplot_simple <- function(x, at, col, width = 0.3) {
  d <- stats::density(x, adjust = 1)
  d$y <- d$y / max(d$y) * width
  graphics::polygon(c(at - d$y, rev(at + d$y)),
                    c(d$x, rev(d$x)),
                    col = col, border = NA)
  graphics::segments(at, stats::quantile(x, 0.25),
                     at, stats::quantile(x, 0.75), lwd = 2)
  graphics::points(at, median(x), pch = 16, cex = 0.8)
}



#' Visualize Basis Sharing Patterns
#'
#' Creates visualizations showing which outputs share which basis functions.
#'
#' @param obj An adaptive_jkhaos object
#' @param mcmc.use Which MCMC samples to use
#' @param output_labels Labels for outputs
#' @export
plot_sharing_patterns <- function(obj, mcmc.use = NULL, output_labels = NULL) {

  if (is.null(mcmc.use)) {
    mcmc.use <- seq_along(obj$nbasis)
  }

  K <- obj$K
  if (is.null(output_labels)) {
    output_labels <- if (!is.null(colnames(obj$Y))) colnames(obj$Y) else paste0("Y", 1:K)
  }

  # Build sharing matrix: entry (k1, k2) = how often outputs k1 and k2 share a basis
  share_matrix <- matrix(0, K, K)
  rownames(share_matrix) <- colnames(share_matrix) <- output_labels

  # Track which specific outputs use each basis
  output_usage <- matrix(0, K, K)

  for (ii in mcmc.use) {
    M <- obj$nbasis[ii]
    if (M > 0) {
      for (mm in 1:M) {
        key <- obj$rkeys[ii, mm]
        if (is.na(key) || !nzchar(key)) next

        r <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])

        # Update co-occurrence matrix
        for (i in r) {
          for (j in r) {
            share_matrix[i, j] <- share_matrix[i, j] + 1
          }
          output_usage[i, length(r)] <- output_usage[i, length(r)] + 1
        }
      }
    }
  }

  # Normalize
  share_matrix <- share_matrix / length(mcmc.use)
  output_usage <- output_usage / length(mcmc.use)

  # Create visualization
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  graphics::par(mfrow = c(1, 2), mar = c(5, 5, 3, 2))

  # Panel 1: Heatmap of sharing
  image(1:K, 1:K, share_matrix,
        col = grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
        xlab = "Output", ylab = "Output",
        main = "Basis Sharing Frequency",
        axes = FALSE)
  axis(1, at = 1:K, labels = output_labels, las = 2)
  axis(2, at = 1:K, labels = output_labels, las = 2)

  # Add text values
  for (i in 1:K) {
    for (j in 1:K) {
      if (share_matrix[i,j] > 0) {
        graphics::text(i, j, sprintf("%.1f", share_matrix[i,j]),
                       col = if(share_matrix[i,j] > max(share_matrix)/2) "white" else "black",
                       cex = 0.8)
      }
    }
  }

  # Panel 2: Distribution of subset sizes per output
  barplot(output_usage, beside = TRUE,
          col = grDevices::rainbow(K, alpha = 0.7),
          xlab = "Subset Size", ylab = "Average Count",
          main = "Basis Functions by Subset Size",
          legend.text = output_labels,
          args.legend = list(x = "topright", bty = "n", cex = 0.7))

  # Return the matrices for further analysis
  invisible(list(
    share_matrix = share_matrix,
    output_usage = output_usage
  ))
}

#' Network Visualization of Output Relationships
#'
#' Creates a network graph showing which outputs are connected via shared bases.
#' Requires igraph package.
#'
#' @param obj An adaptive_jkhaos object
#' @param mcmc.use Which MCMC samples to use
#' @param threshold Minimum sharing frequency to draw an edge
#' @param output_labels Labels for outputs
#' @export
plot_output_network <- function(obj, mcmc.use = NULL, threshold = 0.1,
                                output_labels = NULL) {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for network visualization.")
  }

  sharing <- plot_sharing_patterns(obj, mcmc.use, output_labels)
  share_mat <- sharing$share_matrix

  # Zero out diagonal and apply threshold
  diag(share_mat) <- 0
  share_mat[share_mat < threshold] <- 0

  # Create graph
  g <- igraph::graph_from_adjacency_matrix(share_mat,
                                           mode = "undirected",
                                           weighted = TRUE,
                                           diag = FALSE)

  # Plot
  igraph::plot.igraph(g,
                      vertex.size = 30,
                      vertex.label.cex = 0.8,
                      edge.width = igraph::E(g)$weight * 5,
                      edge.color = "gray50",
                      main = "Output Relationship Network\n(edge width = sharing frequency)")

  invisible(g)
}
