

normalize_tensor_LD <- function(D_rs, r_order, s_order, Phi_A, Phi_B, QA, QB) {
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    return(list(D_norm = NA_real_, D_plus = NA_real_, D_minus = NA_real_))
  }
  
  c_ploidy <- length(QA) - 1
  cost_mat <- outer(Phi_A[, r_order + 1], Phi_B[, s_order + 1])
  
  row_signs <- rep("=", c_ploidy + 1)
  col_signs <- rep("=", c_ploidy + 1)
  
  max_res <- lpSolve::lp.transport(cost.mat = cost_mat, direction = "max", 
                                   row.signs = row_signs, row.rhs = QA, 
                                   col.signs = col_signs, col.rhs = QB, 
                                   integers = NULL)
  min_res <- lpSolve::lp.transport(cost.mat = cost_mat, direction = "min", 
                                   row.signs = row_signs, row.rhs = QA, 
                                   col.signs = col_signs, col.rhs = QB, 
                                   integers = NULL)
  
  if (max_res$status != 0 || min_res$status != 0) {
    return(list(D_norm = NA_real_, D_plus = NA_real_, D_minus = NA_real_))
  }
  
  D_plus <- max_res$objval
  D_minus <- min_res$objval
  
  if (D_rs >= 0) {
    D_norm <- ifelse(abs(D_plus) < 1e-12, NA_real_, D_rs / D_plus)
  } else {
    D_norm <- ifelse(abs(D_minus) < 1e-12, NA_real_, D_rs / abs(D_minus))
  }
  
  return(list(D_norm = D_norm, D_plus = D_plus, D_minus = D_minus))
}


AutoLD_Permutation_Test <- function(
    geno,
    c_ploidy,
    type = c("full", "partial"),
    B = 1000L,
    strata = NULL,
    seed = NULL,
    sd_tol = 1e-8,
    min_convergence_rate = 0.95,
    max_attempt_factor = 2,
    obs_res = NULL  # <--- 核心新增：接收已计算的观测值，避免重复 EM
) {
  
  type <- match.arg(type)
  
  if (!is.null(seed)) set.seed(seed)
  
  if (!is.data.frame(geno) || !"class" %in% names(geno)) {
    stop("geno must be a data.frame containing a 'class' column.")
  }
  
  if (c_ploidy < 2L || c_ploidy != as.integer(c_ploidy)) {
    stop("c_ploidy must be an integer >= 2.")
  }
  
  if (B < 99L) {
    stop("B should be at least 99; B >= 999 is recommended.")
  }
  
  # Full 6x or higher tensor is not identifiable from 3-class data.
  if (type == "partial" && c_ploidy >= 3L) {
    stop(
      "The full higher-order tensor is not identifiable from ",
      "three-class coarsened dosage data at hexaploid or higher ploidy."
    )
  }
  
  # Fix the analysis sample across all permutations.
  valid <- !is.na(geno$class)
  geno_use <- geno[valid, , drop = FALSE]
  
  if (nrow(geno_use) < 20L) {
    stop("Too few complete observations.")
  }
  
  if (!is.null(strata)) {
    strata <- strata[valid]
    
    if (length(strata) != nrow(geno_use)) {
      stop("strata must have one value per original observation.")
    }
  }
  
  v <- 2L * c_ploidy
  S_B <- if (type == "full") v + 1L else 3L
  
  class_vec <- as.integer(geno_use$class)
  
  max_class <- S_B^2 - 1L
  if (any(class_vec < 0L | class_vec > max_class)) {
    stop("geno$class contains values outside the valid class range.")
  }
  
  locus_A <- class_vec %/% S_B
  locus_B <- class_vec %% S_B
  
  # Helper: extract D matrix.
  extract_D <- function(res) {
    D <- matrix(NA_real_, c_ploidy, c_ploidy)
    
    for (r in seq_len(c_ploidy)) {
      for (s in seq_len(c_ploidy)) {
        nm <- paste0("D_", r, s)
        
        if (is.null(res[[nm]])) {
          stop("Missing estimator output: ", nm)
        }
        
        D[r, s] <- as.numeric(res[[nm]])
      }
    }
    
    D
  }

  permute_B <- function(x, strata = NULL) {
    if (is.null(strata)) {
      return(x[sample.int(length(x))])
    }
    
    out <- x
    groups <- split(seq_along(x), strata)
    
    for (idx in groups) {
      out[idx] <- x[idx][sample.int(length(idx))]
    }
    
    out
  }
  

  if (is.null(obs_res)) {
    obs_res <- estimate_LD_decoupled(
      geno = geno_use,
      c_ploidy = c_ploidy,
      type = type,
      calc_norm = FALSE
    )
  }
  
  if (!isTRUE(obs_res$em_converged)) {
    stop("EM failed on the observed data.")
  }
  
  obs_D <- extract_D(obs_res)
  
  # Generate B successfully converged permutations.
  null_D_cube <- array(
    NA_real_,
    dim = c(c_ploidy, c_ploidy, B)
  )
  
  valid_b <- 0L
  attempts <- 0L
  max_attempts <- ceiling(max_attempt_factor * B)
  
  while (valid_b < B && attempts < max_attempts) {
    attempts <- attempts + 1L
    
    permuted_B <- permute_B(locus_B, strata)
    
    perm_class <- locus_A * S_B + permuted_B
    perm_geno <- data.frame(class = perm_class)
    
    perm_out <- tryCatch(
      estimate_LD_decoupled(
        geno = perm_geno,
        c_ploidy = c_ploidy,
        type = type,
        calc_norm = FALSE 
      ),
      error = function(e) NULL
    )
    
    if (!is.null(perm_out) && isTRUE(perm_out$em_converged)) {
      valid_b <- valid_b + 1L
      null_D_cube[, , valid_b] <- extract_D(perm_out)
    }
  }
  
  convergence_rate <- valid_b / attempts
  
  if (valid_b < ceiling(min_convergence_rate * B)) {
    stop(
      sprintf(
        "Only %d valid permutations were obtained after %d attempts.",
        valid_b, attempts
      )
    )
  }
  
  null_D_cube <- null_D_cube[, , seq_len(valid_b), drop = FALSE]
  
  # Component-specific empirical null moments and raw P values.
  null_mean <- apply(null_D_cube, c(1, 2), mean, na.rm = TRUE)
  null_sd <- apply(null_D_cube, c(1, 2), sd, na.rm = TRUE)
  
  estimable <- is.finite(null_sd) & null_sd >= sd_tol
  
  Z_obs <- matrix(
    NA_real_,
    c_ploidy,
    c_ploidy
  )
  
  P_raw <- matrix(
    NA_real_,
    c_ploidy,
    c_ploidy
  )
  
  for (r in seq_len(c_ploidy)) {
    for (s in seq_len(c_ploidy)) {
      
      if (!estimable[r, s]) next
      
      null_dist <- null_D_cube[r, s, ]
      
      obs_dev <- abs(obs_D[r, s] - null_mean[r, s])
      null_dev <- abs(null_dist - null_mean[r, s])
      
      Z_obs[r, s] <-
        (obs_D[r, s] - null_mean[r, s]) / null_sd[r, s]
      
      P_raw[r, s] <-
        (1 + sum(null_dev >= obs_dev)) / (valid_b + 1)
    }
  }
  
  # Standardized null cube.
  Z_null_cube <- array(
    NA_real_,
    dim = dim(null_D_cube)
  )
  
  for (r in seq_len(c_ploidy)) {
    for (s in seq_len(c_ploidy)) {
      if (!estimable[r, s]) next
      
      Z_null_cube[r, s, ] <-
        (null_D_cube[r, s, ] - null_mean[r, s]) /
        null_sd[r, s]
    }
  }
  
  order_index <- seq_len(c_ploidy)
  
  higher_mask <- outer(
    order_index,
    order_index,
    FUN = function(r, s) r > 1L | s > 1L
  )
  
  higher_mask <- higher_mask & estimable
  all_mask <- estimable
  
  if (!any(higher_mask)) {
    stop("No higher-order component is estimable.")
  }
  
  obs_T_HO <- max(abs(Z_obs[higher_mask]), na.rm = TRUE)
  obs_T_all <- max(abs(Z_obs[all_mask]), na.rm = TRUE)
  
  null_T_HO <- numeric(valid_b)
  null_T_all <- numeric(valid_b)
  
  for (b in seq_len(valid_b)) {
    Z_b <- Z_null_cube[, , b]
    
    null_T_HO[b] <-
      max(abs(Z_b[higher_mask]), na.rm = TRUE)
    
    null_T_all[b] <-
      max(abs(Z_b[all_mask]), na.rm = TRUE)
  }
  
  P_HigherOrder <-
    (1 + sum(null_T_HO >= obs_T_HO)) /
    (valid_b + 1)
  
  P_AllTensor <-
    (1 + sum(null_T_all >= obs_T_all)) /
    (valid_b + 1)
  
  # MaxT-adjusted component P values.
  P_MaxT <- matrix(
    NA_real_,
    c_ploidy,
    c_ploidy
  )
  
  for (r in seq_len(c_ploidy)) {
    for (s in seq_len(c_ploidy)) {
      if (!estimable[r, s]) next
      
      reference_T <-
        if (higher_mask[r, s]) null_T_HO else null_T_all
      
      P_MaxT[r, s] <-
        (1 + sum(reference_T >= abs(Z_obs[r, s]))) /
        (valid_b + 1)
    }
  }
  
  list(
    D_observed = obs_D,
    Z_calibrated = Z_obs,
    P_D11 = P_raw[1, 1],
    P_HigherOrder = P_HigherOrder,
    P_AllTensor = P_AllTensor,
    P_rs_raw = P_raw,
    P_rs_MaxT = P_MaxT,
    null_mean = null_mean,
    null_sd = null_sd,
    estimable_components = estimable,
    valid_permutations = valid_b,
    permutation_attempts = attempts,
    convergence_rate = convergence_rate
  )
}