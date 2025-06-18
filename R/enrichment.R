

enrichment_0 <- function(p, d, q, A_curr){
  A_new <- list(rep(0, p)) # Initialize with intercept

  # Check to see if we need to filter based on q (Yes, if q_curr >= q)
  q_curr <- max(apply(A_curr, 1, function(alpha) sum(alpha > 0)))

  # Loop over alpha vectors
  for(i in seq_len(nrow(A_curr))){
    alpha <- A_curr[i,]

    # Start tracking the enriched versions
    enriched_list <- list(alpha)

    for(j in seq_len(p)){
      # Do deletion and bind to A_new (set to 0 all values >= 1)
      if(alpha[j] >= 1){
        tmp <- alpha
        tmp[j] <- 0
        enriched_list <- append(enriched_list, list(tmp))
      }

      # Do simplification and bind to A_new (-1 to all values >= 2)
      if(alpha[j] >= 2){
        tmp <- alpha
        tmp[j] <- tmp[j] - 1
        enriched_list <- append(enriched_list, list(tmp))
      }

      # Do promotion and bind to A_new (set to d, all values < d)
      if(alpha[j] <= d-1){
        tmp <- alpha
        tmp[j] <- d
        enriched_list <- append(enriched_list, list(tmp))
      }

      # Do complication and bind to A_new (+1 to all values < d-1)
      if(alpha[j] <= d-2){
        tmp <- alpha
        tmp[j] <- tmp[j] + 1
        enriched_list <- append(enriched_list, list(tmp))
      }
    }
    A_new <- append(A_new, enriched_list)
  }

  A_new_mat <- unique(do.call(rbind, A_new))   # Remove duplicates
  A_new_mat <- A_new_mat[-1,,drop=FALSE]       # Remove intercept

  # Screen for max order only if necesarry
  # This works because effective q can be at most q_curr + 1.
  if(q_curr >= q){
    valid_idx <- rowSums(A_new_mat > 0) <= q
    A_new_mat <- A_new_mat[valid_idx, , drop = FALSE]
  }
  return(A_new_mat)
}

# Helper to generate the candidates only using certain active vars
generate_A_active <- function(p, d, q, active_vars){
  p_active <- length(active_vars)

  # Generate basis functions in reduced active variable space
  A_reduced <- generate_A(p_active, d, q)
  M <- nrow(A_reduced)

  # Embed into full p-dimensional space
  A_full <- matrix(0L, nrow = M, ncol = p)
  A_full[, active_vars] <- A_reduced

  return(A_full)
}

# Shao version
enrichment_1 <- function(p, d, q, A_curr){
  # Find the active vars
  activity <- apply(A_curr, 2, sum)
  index_active <- which(activity > 0)

  A_new <- generate_A_active(p, d, q, index_active)
  return(A_new)
}

enrichment_2 <- function(p, d, q, A_curr){
  # Find the active vars
  activity <- apply(A_curr, 2, sum)
  index_active <- which(activity > 0)
  index_inactive <- (1:p)[-index_active]

  # Check to see if we need to filter based on q (Yes, if q_curr >= q)
  q_curr <- max(apply(A_curr, 1, function(alpha) sum(alpha > 0)))

  # Build full-active
  A_new <- list(generate_A_active(p, d, q, index_active))

  # Return early if everything is active
  if(length(index_inactive) == 0)
    return(A_new[[1]])

  # Enrich inactives around the previously accepted terms
  for(i in seq_len(nrow(A_curr))){
    alpha <- A_curr[i,]
    for(j in index_inactive){
      # Do promotion and bind to A_new (set to d, all values < d)
      if(alpha[j] <= d-1){
        tmp <- alpha
        tmp[j] <- d
        A_new <- append(A_new, list(tmp))
      }

      # Do complication and bind to A_new (+1 to all values < d-1)
      if(alpha[j] <= d-2){
        tmp <- alpha
        tmp[j] <- tmp[j] + 1
        A_new <- append(A_new, list(tmp))
      }
    }
  }

  # Combine original + enriched
  A_new_mat <- do.call(rbind, A_new)

  # Screen for max order only if necessary
  # This works because effective-q can be at most q_curr + 1.
  if(q_curr >= q){
    valid_idx <- rowSums(A_new_mat > 0) <= q
    A_new_mat <- A_new_mat[valid_idx, , drop = FALSE]
  }
  return(A_new_mat)
}

enrichment_3 <- function(p, d, q, A_curr){
  # Find the active vars
  activity <- apply(A_curr, 2, sum)
  index_active <- which(activity > 0)
  index_inactive <- (1:p)[-index_active]

  # Check to see if we need to filter based on q (Yes, if q_curr >= q)
  q_curr <- max(apply(A_curr, 1, function(alpha) sum(alpha > 0)))

  # Build full-active
  A_base <- generate_A_active(p, d, q, index_active)
  A_new <- list(A_base)

  # Return early if everything is active
  if(length(index_inactive) == 0)
    return(A_base)

  # Enrich inactives around the previously accepted terms
  for(i in seq_len(nrow(A_base))){
    alpha <- A_base[i,]
    for(j in index_inactive){
      # Do promotion and bind to A_new (set to d, all values < d)
      if(alpha[j] <= d-1){
        tmp <- alpha
        tmp[j] <- d
        A_new <- append(A_new, list(tmp))
      }

      # Do complication and bind to A_new (+1 to all values < d-1)
      if(alpha[j] <= d-2){
        tmp <- alpha
        tmp[j] <- tmp[j] + 1
        A_new <- append(A_new, list(tmp))
      }
    }
  }

  # Combine original + enriched
  A_new_mat <- do.call(rbind, A_new)

  # Screen for max order only if necessary
  # This works because effective-q can be at most q_curr + 1.
  if(q_curr >= q){
    valid_idx <- rowSums(A_new_mat > 0) <= q
    A_new_mat <- A_new_mat[valid_idx, , drop = FALSE]
  }
  return(A_new_mat)
}

# Generate everything; A_curr is ignored
enrichment_4 <- function(p, d, q, A_curr){
  generate_A(p, d, q)
}

