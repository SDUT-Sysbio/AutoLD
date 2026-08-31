






AutoLD_GenomeScan <- function(chr_name, geno_df, info_df, 
                              window_kb = 500, ploidy = 4, type = "full",
                              calc_norm = TRUE, calc_pvalue = FALSE, 
                              pvalue_trigger_r2 = NULL, pvalue_trigger_HO = NULL, 
                              B_perm = 200, n_cores = 1, min_N_eff = 30) {
  
  cat(sprintf("\n[%s] Starting AutoLD Genomic Scan...\n", chr_name))
  cat(sprintf(" | Ploidy: %dx | Type: %s | Window: %d kb\n", ploidy, type, window_kb))
  cat(sprintf(" | Calc Norm (D'): %s | Calc P-value: %s\n", calc_norm, calc_pvalue))
  
  if (calc_pvalue) {
    if (!is.null(pvalue_trigger_r2)) cat(sprintf(" | [Trigger 1] Permutation if r2_comp >= %.3f\n", pvalue_trigger_r2))
    if (!is.null(pvalue_trigger_HO)) cat(sprintf(" | [Trigger 2] Permutation if Max|Higher-Order D| >= %.3f\n", pvalue_trigger_HO))
  }
  
  c_ploidy <- ploidy / 2
  v_ploidy <- ploidy
  
  idx_chr <- which(info_df$CHROM == chr_name)
  if (length(idx_chr) == 0) stop(sprintf("No markers found on chromosome %s", chr_name))
  
  sub_info <- info_df[idx_chr, ]
  sub_geno <- as.matrix(geno_df[idx_chr, ])
  
  ord <- order(sub_info$POS)
  sub_info <- sub_info[ord, ]
  sub_geno <- sub_geno[ord, ]
  positions <- sub_info$POS
  
  end_idxs <- findInterval(positions + (window_kb * 1000), positions)
  valid_mask <- (end_idxs - seq_along(positions)) > 0
  
  tasks_list <- lapply(which(valid_mask), function(i) {
    j_range <- (i + 1):end_idxs[i]
    cbind(i, j_range)
  })
  tasks_mat <- do.call(rbind, tasks_list)
  tasks <- split(tasks_mat, row(tasks_mat))
  
  cat(sprintf("[%s] Total SNP pairs to evaluate: %d\n", chr_name, length(tasks)))
  

  calc_one_pair <- function(pair_idx) {
    i <- pair_idx[1]; j <- pair_idx[2]
    gA <- sub_geno[i, ]
    gB <- sub_geno[j, ]
    

    valid <- !is.na(gA) & !is.na(gB)
    N_eff <- sum(valid)
    if (N_eff < min_N_eff) return(NULL) 
    

    if (type == "partial") {
      gA_use <- ifelse(gA[valid] == 0, 0, ifelse(gA[valid] == v_ploidy, 2, 1))
      gB_use <- ifelse(gB[valid] == 0, 0, ifelse(gB[valid] == v_ploidy, 2, 1))
      geno_input <- data.frame(class = as.integer(gA_use * 3 + gB_use))
    } else {
      geno_input <- data.frame(class = as.integer(gA[valid] * (v_ploidy + 1) + gB[valid]))
    }
    

    base_res <- tryCatch({
      estimate_LD_decoupled(geno_input, c_ploidy, type = type, calc_norm = calc_norm)
    }, error = function(e) return(NULL))
    
    if (is.null(base_res) || !isTRUE(base_res$em_converged)) return(NULL)
    

    out_dict <- list(
      CHR = chr_name, POS_A = positions[i], POS_B = positions[j], Dist_bp = positions[j] - positions[i],
      N_eff = N_eff,
      Eta_A = base_res$eta_A_est, Eta_B = base_res$eta_B_est,
      D_comp = base_res$D_comp, r_comp = base_res$r_comp, r2_comp = base_res$r2_comp
    )
    

    for(r in 1:c_ploidy) {
      for(s in 1:c_ploidy) {
        out_dict[[paste0("D_raw_", r, s)]] <- base_res[[paste0("D_", r, s)]]
        if (calc_norm) {
          out_dict[[paste0("D_norm_", r, s)]] <- base_res[[paste0("D_norm_", r, s)]]
          out_dict[[paste0("D_plus_", r, s)]] <- base_res[[paste0("D_plus_", r, s)]]
          out_dict[[paste0("D_minus_", r, s)]] <- base_res[[paste0("D_minus_", r, s)]]
        }
      }
    }
    
    out_dict$Permutation_tested <- FALSE
    out_dict$Permutation_status <- NA_character_
    out_dict$P_D11 <- NA_real_; out_dict$P_HO_Omnibus <- NA_real_; out_dict$P_All_Tensor <- NA_real_
    for(r in 1:c_ploidy) for(s in 1:c_ploidy) out_dict[[paste0("P_MaxT_", r, s)]] <- NA_real_
    
    if (calc_pvalue) {
      trigger_D11 <- FALSE
      if (!is.null(pvalue_trigger_r2) && !is.na(base_res$r2_comp)) {
        trigger_D11 <- base_res$r2_comp >= pvalue_trigger_r2
      }
      
      trigger_HO <- FALSE
      if (!is.null(pvalue_trigger_HO)) {
        ho_values <- numeric(0)
        for (r in seq_len(c_ploidy)) {
          for (s in seq_len(c_ploidy)) {
            if (r == 1L && s == 1L) next
            nm <- if (calc_norm) paste0("D_norm_", r, s) else paste0("D_", r, s)
            ho_values <- c(ho_values, base_res[[nm]])
          }
        }
        if (length(ho_values) > 0) {
          trigger_HO <- max(abs(ho_values), na.rm = TRUE) >= pvalue_trigger_HO
        }
      }
      
      run_permutation <- FALSE
      if (is.null(pvalue_trigger_r2) && is.null(pvalue_trigger_HO)) {
        run_permutation <- TRUE # 未设阈值则无条件强制触发
      } else {
        run_permutation <- trigger_D11 || trigger_HO
      }
      
      if (run_permutation) {
        perm_res <- tryCatch({
          AutoLD_Permutation_Test(geno = geno_input, c_ploidy = c_ploidy, type = type, 
                                  B = B_perm, obs_res = base_res)
        }, error = function(e) return(NULL))
        
        if (!is.null(perm_res)) {
          out_dict$Permutation_tested <- TRUE
          out_dict$Permutation_status <- "success"
          out_dict$P_D11 <- perm_res$P_D11
          out_dict$P_HO_Omnibus <- perm_res$P_HigherOrder
          out_dict$P_All_Tensor <- perm_res$P_AllTensor
          for(r in 1:c_ploidy) for(s in 1:c_ploidy) out_dict[[paste0("P_MaxT_", r, s)]] <- perm_res$P_rs_MaxT[r, s]
        } else {
          out_dict$Permutation_tested <- TRUE
          out_dict$Permutation_status <- "failed"
        }
      } else {
        out_dict$Permutation_tested <- FALSE
        out_dict$Permutation_status <- "below_screening_threshold"
      }
    }
    
    return(out_dict)
  }
  
  cat(sprintf("[%s] Executing parallel scan on %d cores...\n", chr_name, n_cores))
  if (Sys.info()[['sysname']] == "Windows") {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, varlist=ls(envir=globalenv()), envir=globalenv())
    res_list <- parallel::parLapply(cl, tasks, calc_one_pair)
    parallel::stopCluster(cl)
  } else {
    res_list <- parallel::mclapply(tasks, calc_one_pair, mc.cores = n_cores)
  }
  
  res_list <- res_list[!sapply(res_list, is.null)]
  if (length(res_list) == 0) {
    warning(sprintf("[%s] No valid pairs successfully calculated.", chr_name))
    return(data.frame())
  }
  
  final_df <- as.data.frame(data.table::rbindlist(res_list, fill = TRUE))
  cat(sprintf("[%s] Scan complete. Retained %d valid pairs.\n", chr_name, nrow(final_df)))
  
  return(final_df)
}


AutoLD_GenomeScan_can <- function(chr_name, geno_df, info_df, 
                              window_kb = 500, ploidy = 4, type = "full",
                              calc_norm = TRUE, calc_pvalue = FALSE, 
                              pvalue_trigger_r2 = NULL, pvalue_trigger_HO = NULL, 
                              B_perm = 200, n_cores = 1, min_N_eff = 30,
                              candidate_pairs = NULL,
                              force_pvalue = FALSE) {
  
  cat(sprintf("\n[%s] Starting AutoLD Genomic Scan...\n", chr_name))
  cat(sprintf(" | Ploidy: %dx | Type: %s | Window: %d kb\n", ploidy, type, window_kb))
  cat(sprintf(" | Calc Norm (D'): %s | Calc P-value: %s\n", calc_norm, calc_pvalue))
  cat(sprintf(" | Candidate-only mode: %s\n", !is.null(candidate_pairs)))
  cat(sprintf(" | Force permutation: %s\n", force_pvalue))
  
  if (calc_pvalue && !force_pvalue) {
    if (!is.null(pvalue_trigger_r2)) {
      cat(sprintf(" | [Trigger 1] Permutation if r2_comp >= %.3f\n", pvalue_trigger_r2))
    }
    if (!is.null(pvalue_trigger_HO)) {
      cat(sprintf(" | [Trigger 2] Permutation if Max|Higher-Order D| >= %.3f\n", pvalue_trigger_HO))
    }
  }
  
  c_ploidy <- ploidy / 2
  v_ploidy <- ploidy
  
  if (!"CHROM" %in% names(info_df)) {
    stop("info_df must contain column 'CHROM'.")
  }
  if (!"POS" %in% names(info_df)) {
    stop("info_df must contain column 'POS'.")
  }
  
  idx_chr <- which(info_df$CHROM == chr_name)
  if (length(idx_chr) == 0) {
    stop(sprintf("No markers found on chromosome %s", chr_name))
  }
  
  sub_info <- as.data.frame(info_df[idx_chr, ])
  sub_geno <- as.matrix(geno_df[idx_chr, ])
  
  ord <- order(sub_info$POS)
  sub_info <- sub_info[ord, , drop = FALSE]
  sub_geno <- sub_geno[ord, , drop = FALSE]
  positions <- sub_info$POS
  
  if (!is.null(candidate_pairs)) {
    
    cand <- as.data.table(candidate_pairs)
    
    if (!"CHROM" %in% names(cand)) {
      if ("CHR" %in% names(cand)) {
        setnames(cand, "CHR", "CHROM")
      } else {
        cand[, CHROM := chr_name]
      }
    }
    
    cand <- cand[CHROM == chr_name]
    
    if (nrow(cand) == 0) {
      cat(sprintf("[%s] No candidate pairs for this chromosome.\n", chr_name))
      return(data.frame())
    }
    
    if (!all(c("POS_A", "POS_B") %in% names(cand))) {
      stop("candidate_pairs must contain POS_A and POS_B columns.")
    }
    
    pos_map <- data.table(
      POS = positions,
      local_idx = seq_along(positions)
    )
    
    pos_map <- unique(pos_map, by = "POS")
    setkey(pos_map, POS)
    
    cand[, idx_A := pos_map[.(POS_A), local_idx]]
    cand[, idx_B := pos_map[.(POS_B), local_idx]]
    
    n_miss_A <- sum(is.na(cand$idx_A))
    n_miss_B <- sum(is.na(cand$idx_B))
    
    if (n_miss_A > 0 || n_miss_B > 0) {
      warning(sprintf(
        "[%s] Missing position mapping: A=%d, B=%d. These pairs are removed.",
        chr_name, n_miss_A, n_miss_B
      ))
      cand <- cand[!is.na(idx_A) & !is.na(idx_B)]
    }
    
    if (nrow(cand) == 0) {
      cat(sprintf("[%s] No candidate pairs after POS mapping.\n", chr_name))
      return(data.frame())
    }
    
    cand <- cand[idx_A != idx_B]
    
    cand[, i := pmin(idx_A, idx_B)]
    cand[, j := pmax(idx_A, idx_B)]
    
    cand <- unique(cand, by = c("i", "j"))
    
    tasks_mat <- as.matrix(cand[, .(i, j)])
    tasks <- split(tasks_mat, row(tasks_mat))
    
    cat(sprintf("[%s] Candidate pairs to evaluate: %d\n", chr_name, length(tasks)))
    
  } else {
    
    end_idxs <- findInterval(positions + (window_kb * 1000), positions)
    valid_mask <- (end_idxs - seq_along(positions)) > 0
    
    tasks_list <- lapply(which(valid_mask), function(i) {
      j_range <- (i + 1):end_idxs[i]
      cbind(i, j_range)
    })
    
    tasks_mat <- do.call(rbind, tasks_list)
    tasks <- split(tasks_mat, row(tasks_mat))
    
    cat(sprintf("[%s] Total SNP pairs to evaluate: %d\n", chr_name, length(tasks)))
  }
  
  if (length(tasks) == 0) {
    warning(sprintf("[%s] No SNP pairs to evaluate.", chr_name))
    return(data.frame())
  }
  
  ## ------------------------------------------------------------
  ## 3. Single-pair worker
  ## ------------------------------------------------------------
  calc_one_pair <- function(pair_idx) {
    
    i <- pair_idx[1]
    j <- pair_idx[2]
    
    gA <- sub_geno[i, ]
    gB <- sub_geno[j, ]
    
    valid <- !is.na(gA) & !is.na(gB)
    N_eff <- sum(valid)
    
    if (N_eff < min_N_eff) return(NULL)
    
    if (type == "partial") {
      gA_use <- ifelse(gA[valid] == 0, 0, ifelse(gA[valid] == v_ploidy, 2, 1))
      gB_use <- ifelse(gB[valid] == 0, 0, ifelse(gB[valid] == v_ploidy, 2, 1))
      geno_input <- data.frame(class = as.integer(gA_use * 3 + gB_use))
    } else {
      geno_input <- data.frame(
        class = as.integer(gA[valid] * (v_ploidy + 1) + gB[valid])
      )
    }
    
    base_res <- tryCatch({
      estimate_LD_decoupled(
        geno_input,
        c_ploidy,
        type = type,
        calc_norm = calc_norm
      )
    }, error = function(e) return(NULL))
    
    if (is.null(base_res) || !isTRUE(base_res$em_converged)) return(NULL)
    
    out_dict <- list(
      CHR = chr_name,
      POS_A = positions[i],
      POS_B = positions[j],
      Dist_bp = positions[j] - positions[i],
      N_eff = N_eff,
      Eta_A = base_res$eta_A_est,
      Eta_B = base_res$eta_B_est,
      D_comp = base_res$D_comp,
      r_comp = base_res$r_comp,
      r2_comp = base_res$r2_comp
    )
    
    for (r in 1:c_ploidy) {
      for (s in 1:c_ploidy) {
        out_dict[[paste0("D_raw_", r, s)]] <- base_res[[paste0("D_", r, s)]]
        
        if (calc_norm) {
          out_dict[[paste0("D_norm_", r, s)]] <- base_res[[paste0("D_norm_", r, s)]]
          out_dict[[paste0("D_plus_", r, s)]] <- base_res[[paste0("D_plus_", r, s)]]
          out_dict[[paste0("D_minus_", r, s)]] <- base_res[[paste0("D_minus_", r, s)]]
        }
      }
    }
    
    out_dict$Permutation_tested <- FALSE
    out_dict$Permutation_status <- NA_character_
    out_dict$P_D11 <- NA_real_
    out_dict$P_HO_Omnibus <- NA_real_
    out_dict$P_All_Tensor <- NA_real_
    
    for (r in 1:c_ploidy) {
      for (s in 1:c_ploidy) {
        out_dict[[paste0("P_MaxT_", r, s)]] <- NA_real_
      }
    }
    
    if (calc_pvalue) {
      
      if (isTRUE(force_pvalue)) {
        
        run_permutation <- TRUE
        
      } else {
        
        trigger_D11 <- FALSE
        
        if (!is.null(pvalue_trigger_r2) && !is.na(base_res$r2_comp)) {
          trigger_D11 <- base_res$r2_comp >= pvalue_trigger_r2
        }
        
        trigger_HO <- FALSE
        
        if (!is.null(pvalue_trigger_HO)) {
          ho_values <- numeric(0)
          
          for (r in seq_len(c_ploidy)) {
            for (s in seq_len(c_ploidy)) {
              if (r == 1L && s == 1L) next
              
              nm <- if (calc_norm) {
                paste0("D_norm_", r, s)
              } else {
                paste0("D_", r, s)
              }
              
              ho_values <- c(ho_values, base_res[[nm]])
            }
          }
          
          if (length(ho_values) > 0) {
            trigger_HO <- max(abs(ho_values), na.rm = TRUE) >= pvalue_trigger_HO
          }
        }
        
        if (is.null(pvalue_trigger_r2) && is.null(pvalue_trigger_HO)) {
          run_permutation <- TRUE
        } else {
          run_permutation <- trigger_D11 || trigger_HO
        }
      }
      
      if (run_permutation) {
        
        perm_res <- tryCatch({
          AutoLD_Permutation_Test(
            geno = geno_input,
            c_ploidy = c_ploidy,
            type = type,
            B = B_perm,
            obs_res = base_res
          )
        }, error = function(e) return(NULL))
        
        if (!is.null(perm_res)) {
          out_dict$Permutation_tested <- TRUE
          out_dict$Permutation_status <- "success"
          out_dict$P_D11 <- perm_res$P_D11
          out_dict$P_HO_Omnibus <- perm_res$P_HigherOrder
          out_dict$P_All_Tensor <- perm_res$P_AllTensor
          
          for (r in 1:c_ploidy) {
            for (s in 1:c_ploidy) {
              out_dict[[paste0("P_MaxT_", r, s)]] <- perm_res$P_rs_MaxT[r, s]
            }
          }
          
        } else {
          out_dict$Permutation_tested <- TRUE
          out_dict$Permutation_status <- "failed"
        }
        
      } else {
        out_dict$Permutation_tested <- FALSE
        out_dict$Permutation_status <- "below_screening_threshold"
      }
    }
    
    return(out_dict)
  }
  
  cat(sprintf("[%s] Executing parallel scan on %d cores...\n", chr_name, n_cores))
  
  if (Sys.info()[["sysname"]] == "Windows") {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterExport(
      cl,
      varlist = ls(envir = globalenv()),
      envir = globalenv()
    )
    
    parallel::clusterExport(
      cl,
      varlist = c(
        "sub_geno", "positions", "tasks", "calc_one_pair",
        "ploidy", "type", "calc_norm", "calc_pvalue",
        "pvalue_trigger_r2", "pvalue_trigger_HO",
        "B_perm", "min_N_eff", "force_pvalue",
        "c_ploidy", "v_ploidy"
      ),
      envir = environment()
    )
    
    res_list <- parallel::parLapply(cl, tasks, calc_one_pair)
    
  } else {
    res_list <- parallel::mclapply(
      tasks,
      calc_one_pair,
      mc.cores = n_cores
    )
  }
  
  res_list <- res_list[!sapply(res_list, is.null)]
  
  if (length(res_list) == 0) {
    warning(sprintf("[%s] No valid pairs successfully calculated.", chr_name))
    return(data.frame())
  }
  
  final_df <- as.data.frame(data.table::rbindlist(res_list, fill = TRUE))
  
  cat(sprintf("[%s] Scan complete. Retained %d valid pairs.\n", chr_name, nrow(final_df)))
  
  return(final_df)
}
