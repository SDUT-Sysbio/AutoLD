

simulate_unified_zygotes <- function(QA, QB, D_matrix, c_ploidy,
                                     n = 500,
                                     type = "full",
                                     hwd_f = 0,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!type %in% c("full", "partial")) stop("type must be either 'full' or 'partial'.")
  QA <- as.numeric(QA); QB <- as.numeric(QB)
  QA[QA < 0] <- 0; QB[QB < 0] <- 0
  QA <- QA / sum(QA); QB <- QB / sum(QB)
  
  Phi_A <- build_orthogonal_basis(QA, c_ploidy)
  Phi_B <- build_orthogonal_basis(QB, c_ploidy)
  P_AB <- matrix(0, nrow = c_ploidy + 1, ncol = c_ploidy + 1)
  
  for (i in 0:c_ploidy) {
    for (j in 0:c_ploidy) {
      interaction <- 0
      for (r in 1:c_ploidy) {
        for (s in 1:c_ploidy) {
          interaction <- interaction + D_matrix[r, s] * Phi_A[i + 1, r + 1] * Phi_B[j + 1, s + 1]
        }
      }
      P_AB[i + 1, j + 1] <- QA[i + 1] * QB[j + 1] * (1 + interaction)
    }
  }
  
  if (min(P_AB) < -1e-12) {
    stop(sprintf(
      "Infeasible D matrix: minimum P_AB = %.6g. The injected higher-order tensors violate the simplex bounds for this allele frequency.",
      min(P_AB)
    ))
  }
  P_AB[P_AB < 0] <- 0
  P_AB <- P_AB / sum(P_AB)
  
  # 4. 终极闭环验收 (Mathematical Loop-back Test)
  D_ext <- extract_Drs_tensor(P_AB, c_ploidy)
  stopifnot(
    max(abs(rowSums(P_AB) - QA)) < 1e-10,
    max(abs(colSums(P_AB) - QB)) < 1e-10,
    max(abs(D_ext - D_matrix)) < 1e-8
  )
  
  gamete_types <- expand.grid(LocusA = 0:c_ploidy, LocusB = 0:c_ploidy)
  gamete_probs <- as.vector(P_AB)
  
  
  idx1 <- sample(seq_len(nrow(gamete_types)), n, replace = TRUE, prob = gamete_probs)
  
  is_inbred <- runif(n) < hwd_f 
  idx2_random <- sample(seq_len(nrow(gamete_types)), n, replace = TRUE, prob = gamete_probs)
  idx2 <- ifelse(is_inbred, idx1, idx2_random)
  
  g1 <- gamete_types[idx1, ]
  g2 <- gamete_types[idx2, ]
  
  zygote_A <- g1$LocusA + g2$LocusA
  zygote_B <- g1$LocusB + g2$LocusB
  
  v_ploidy <- c_ploidy * 2
  full_class <- as.integer(zygote_A * (v_ploidy + 1) + zygote_B)
  
  collapse_dosage <- function(x, v) { ifelse(x == 0, 0, ifelse(x == v, 2, 1)) }
  
  if (type == "partial") {
    obs_A <- collapse_dosage(zygote_A, v_ploidy)
    obs_B <- collapse_dosage(zygote_B, v_ploidy)
    obs_class <- as.integer(obs_A * 3 + obs_B)
  } else {
    obs_class <- full_class
  }
  
  return(list(
    class = obs_class,
    full_class = full_class,
    type = type,
    zygote_A = zygote_A,
    zygote_B = zygote_B,
    P_AB_true = P_AB,
    QA_true = QA,
    QB_true = QB,
    D_tensor_input = D_matrix,
    D_tensor_true = extract_Drs_tensor(P_AB, c_ploidy),
    D_eab_true = calc_Deab_comp(P_AB)
  ))
}

simulate_from_eta <- function(v_ploidy, pA, pB, eta1_A, eta1_B, D_matrix, 
                              n = 500, type = "full", 
                              hwd_f = 0,
                              seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  c_ploidy <- v_ploidy / 2
  
  QA_true <- get_equilibrium_Q_one_eta(v = v_ploidy, p = pA, eta1 = eta1_A)
  QB_true <- get_equilibrium_Q_one_eta(v = v_ploidy, p = pB, eta1 = eta1_B)
  
  sim_res <- simulate_unified_zygotes(
    QA = QA_true, 
    QB = QB_true, 
    D_matrix = D_matrix, 
    c_ploidy = c_ploidy, 
    n = n,
    type = type,
    hwd_f = hwd_f,
    seed = seed
  )
  
  sim_res$v_ploidy <- v_ploidy
  sim_res$pA_true <- pA
  sim_res$pB_true <- pB
  sim_res$eta1_A_true <- eta1_A
  sim_res$eta1_B_true <- eta1_B
  sim_res$hwd_f_true <- hwd_f
  
  return(sim_res)
}


evaluate_AutoLD_pipeline <- function(
    reps = 100, n = 500, c_ploidy = 2, type = "full",
    pA = 0.5, pB = 0.5, eta_A = 0, eta_B = 0, D_true_matrix = NULL, hwd_f = 0,
    calc_norm = FALSE, calc_pvalue = FALSE, B_perm = 200, alpha = 0.05, n_cores = 1, seed = 123
) {
  
  set.seed(seed)
  v_ploidy <- 2 * c_ploidy
  if(is.null(D_true_matrix)) D_true_matrix <- matrix(0, c_ploidy, c_ploidy)
  
  has_ldsep <- requireNamespace("ldsep", quietly = TRUE)
  
  cat(sprintf("\n[AutoLD Simulator] Starting %d reps | %dx %s | n=%d | Cores=%d\n", 
              reps, v_ploidy, type, n, n_cores))
  if (calc_pvalue) cat(sprintf(">>> Permutation Test Enabled (B=%d). This will take time...\n", B_perm))
  if (calc_norm) cat(">>> Tensor Normalization (D') Enabled. Using Simplex constraints...\n")
  if (has_ldsep) cat(">>> Baseline Comparison: ldsep::ldest() Enabled.\n")
  

  run_single_sim <- function(i) {
    sim_data <- tryCatch({
      simulate_from_eta(v_ploidy, pA, pB, eta_A, eta_B, D_true_matrix, n, type, hwd_f)
    }, error = function(e) {
      message(sprintf("\n[Rep %d] Simulation Error: %s", i, e$message))
      return(NULL)
    })
    
    if (is.null(sim_data)) return(NULL)
    
    r_ldsep <- NA_real_; r2_ldsep <- NA_real_
    if (has_ldsep) {
      ldsep_res <- tryCatch({
        ldsep::ldest(sim_data$zygote_A, sim_data$zygote_B, K = v_ploidy, type = "comp")
      }, error = function(e) return(NULL))
      
      if (!is.null(ldsep_res)) {
        if ("r" %in% names(ldsep_res)) r_ldsep <- unname(ldsep_res["r"])
        if ("r2" %in% names(ldsep_res)) r2_ldsep <- unname(ldsep_res["r2"])
      }
    }
    # -------------------------------------------------------------
    
    geno_df <- data.frame(class = sim_data$class)
    
    if (calc_pvalue) {
      perm_res <- tryCatch({
        AutoLD_Permutation_Test(geno = geno_df, c_ploidy = c_ploidy, type = type, B = B_perm, sd_tol = 1e-8)
      }, error = function(e) {
        message(sprintf("\n[Rep %d] Permutation Error: %s", i, e$message))
        return(NULL)
      })
      
      if (is.null(perm_res)) return(NULL)
      
      res_flat <- list(em_converged = TRUE)
      for(r in 1:c_ploidy) {
        for(s in 1:c_ploidy) {
          res_flat[[paste0("D_", r, s)]] <- perm_res$D_observed[r,s]
          res_flat[[paste0("P_raw_", r, s)]] <- perm_res$P_rs_raw[r,s]
          res_flat[[paste0("P_MaxT_", r, s)]] <- perm_res$P_rs_MaxT[r,s]
        }
      }
      res_flat$P_D11 <- perm_res$P_D11
      res_flat$P_HO <- perm_res$P_HigherOrder
      res_flat$P_All <- perm_res$P_AllTensor
      
      base_res <- tryCatch({
        estimate_LD_decoupled(geno_df, c_ploidy, type, calc_norm = calc_norm)
      }, error = function(e) return(NULL))
      if(is.null(base_res)) return(NULL)
      
      res_flat$eta_A_est <- base_res$eta_A_est
      res_flat$eta_B_est <- base_res$eta_B_est
      res_flat$D_comp <- base_res$D_comp
      res_flat$r_comp <- base_res$r_comp   # 提取 r
      res_flat$r2_comp <- base_res$r2_comp # 提取 r2
      
      # 挂载 ldsep 结果
      res_flat$r_ldsep <- r_ldsep
      res_flat$r2_ldsep <- r2_ldsep
      
      if (calc_norm) {
        for(r in 1:c_ploidy) for(s in 1:c_ploidy) {
          res_flat[[paste0("D_norm_", r, s)]] <- base_res[[paste0("D_norm_", r, s)]]
        }
      }
      return(res_flat)
      
    } else {
      res <- tryCatch({
        estimate_LD_decoupled(geno_df, c_ploidy, type, calc_norm = calc_norm)
      }, error = function(e) {
        message(sprintf("\n[Rep %d] Estimation Error: %s", i, e$message))
        return(NULL)
      })
      
      if (!is.null(res)) {
        # 挂载 ldsep 结果
        res$r_ldsep <- r_ldsep
        res$r2_ldsep <- r2_ldsep
      }
      return(res)
    }
  }
  
  results <- if (.Platform$OS.type == "windows") {
    if (n_cores > 1L) warning("On Windows, evaluate_AutoLD_pipeline currently uses one core.")
    lapply(1:reps, run_single_sim)
  } else {
    parallel::mclapply(1:reps, run_single_sim, mc.cores = n_cores)
  }
  
  valid_res <- results[!sapply(results, is.null)]
  converged_mask <- vapply(valid_res, function(x) isTRUE(x$em_converged), logical(1))
  converged_res <- valid_res[converged_mask]
  
  if (length(converged_res) == 0) stop("All simulations failed to converge.")
  success_rate <- length(converged_res) / reps
  
  # 防御性提取（防止列表有 NULL 元素导致均值报错）
  extract_vec <- function(var_name) {
    vec <- sapply(converged_res, function(x) x[[var_name]])
    vec <- unlist(vec[!sapply(vec, is.null)])
    return(as.numeric(vec))
  }
  
  cat("\n================ EVALUATION SUMMARY ================\n")
  cat(sprintf("Successful convergence: %d / %d (%.1f%%)\n\n", length(converged_res), reps, success_rate*100))
  
  etaA_est <- extract_vec("eta_A_est"); etaB_est <- extract_vec("eta_B_est")
  cat("[1. Marginal Segregation Distortion]\n")
  cat(sprintf("  Eta_A | True: %.2f | Mean: %.4f | Bias: % .4f | MSE: %.6f\n", 
              eta_A, mean(etaA_est), mean(etaA_est) - eta_A, mean((etaA_est - eta_A)^2)))
  cat(sprintf("  Eta_B | True: %.2f | Mean: %.4f | Bias: % .4f | MSE: %.6f\n\n", 
              eta_B, mean(etaB_est), mean(etaB_est) - eta_B, mean((etaB_est - eta_B)^2)))
  
  # -------------------------------------------------------------
  cat("[2. Classic Composite LD Comparison]\n")
  cat(sprintf("  [AutoLD] D_comp Mean : % .4f\n", mean(extract_vec("D_comp"), na.rm=TRUE)))
  cat(sprintf("  [AutoLD] r_comp Mean : % .4f\n", mean(extract_vec("r_comp"), na.rm=TRUE)))
  cat(sprintf("  [AutoLD] r2_comp Mean: % .4f\n", mean(extract_vec("r2_comp"), na.rm=TRUE)))
  
  if (has_ldsep) {
    cat(sprintf("  [ldsep ] r Mean      : % .4f\n", mean(extract_vec("r_ldsep"), na.rm=TRUE)))
    cat(sprintf("  [ldsep ] r2 Mean     : % .4f\n\n", mean(extract_vec("r2_ldsep"), na.rm=TRUE)))
  } else {
    cat("\n")
  }
  # -------------------------------------------------------------
  
  cat("[3. Tensor Components D_rs]\n")
  out_df <- data.frame()
  for (r in 1:c_ploidy) {
    for (s in 1:c_ploidy) {
      d_est <- extract_vec(paste0("D_", r, s))
      d_true <- D_true_matrix[r, s]
      bias <- mean(d_est, na.rm=TRUE) - d_true
      mse <- mean((d_est - d_true)^2, na.rm=TRUE)
      
      row_info <- data.frame(
        Component = paste0("D_", r, s), True = d_true, Mean_Est = mean(d_est, na.rm=TRUE),
        Bias = bias, MSE = mse
      )
      
      if (calc_norm) {
        d_norm <- extract_vec(paste0("D_norm_", r, s))
        row_info$Mean_Norm <- mean(d_norm, na.rm = TRUE)
      }
      
      if (calc_pvalue) {
        p_raw <- extract_vec(paste0("P_raw_", r, s))
        p_maxt <- extract_vec(paste0("P_MaxT_", r, s))
        row_info$Power_Raw <- mean(p_raw < alpha, na.rm=TRUE)
        row_info$Power_MaxT <- mean(p_maxt < alpha, na.rm=TRUE)
      }
      out_df <- rbind(out_df, row_info)
    }
  }
  print(out_df, row.names = FALSE, digits = 4)
  
  if (calc_pvalue) {
    cat("\n[4. Omnibus Tests (Power / FPR at alpha = 0.05)]\n")
    cat(sprintf("  P_D11 (Linear)    : %.3f\n", mean(extract_vec("P_D11") < alpha, na.rm=TRUE)))
    cat(sprintf("  P_HigherOrder (HO): %.3f\n", mean(extract_vec("P_HO") < alpha, na.rm=TRUE)))
    cat(sprintf("  P_AllTensor (All) : %.3f\n", mean(extract_vec("P_All") < alpha, na.rm=TRUE)))
  }
  cat("====================================================\n")
  
  invisible(list(converged_res = converged_res, summary = out_df))
}