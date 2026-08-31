
build_orthogonal_basis <- function(Q, c_ploidy) {
  
  Q <- as.numeric(Q)
  if (any(Q < -1e-12)) stop("Q contains negative probabilities.")
  Q[Q < 0] <- 0
  Q <- Q / sum(Q)
  
  support_size <- sum(Q > 1e-10)
  max_order_estimable <- support_size - 1
  if (max_order_estimable < c_ploidy) {
    stop(sprintf(
      "Only orders up to %d are estimable from the marginal support. Requested %d. (MAF too low or severe bottleneck).",
      max_order_estimable, c_ploidy
    ))
  }
  
  x <- 0:c_ploidy
  basis <- matrix(0, nrow = c_ploidy + 1, ncol = c_ploidy + 1)
  
  basis[, 1] <- 1
  
  for (r in 1:c_ploidy) {
    v <- x^r 
  
    for (k in 1:r) {
      norm_k <- sum(Q * basis[, k]^2)
      if (norm_k > 1e-12) {
        proj_coef <- sum(Q * v * basis[, k]) / norm_k
        v <- v - proj_coef * basis[, k]
      }
    }
    basis[, r + 1] <- v
  }
  
  for (r in 2:(c_ploidy + 1)) {
    norm_factor <- sqrt(sum(Q * basis[, r]^2))
    if(norm_factor > 1e-12) {
      basis[, r] <- basis[, r] / norm_factor
    }
  }
  
  return(basis)
}

extract_Drs_tensor <- function(P_AB, c_ploidy) {
  
  if (!is.matrix(P_AB)) stop("P_AB must be a matrix.")
  if (nrow(P_AB) != c_ploidy + 1 || ncol(P_AB) != c_ploidy + 1) {
    stop("P_AB dimensions must be (c_ploidy + 1) x (c_ploidy + 1).")
  }
  if (any(P_AB < -1e-12)) stop("P_AB contains negative probabilities.")
  P_AB[P_AB < 0] <- 0
  P_AB <- P_AB / sum(P_AB)
  
  Q_A <- rowSums(P_AB)
  Q_B <- colSums(P_AB)
  
  Phi_A <- build_orthogonal_basis(Q_A, c_ploidy)
  Phi_B <- build_orthogonal_basis(Q_B, c_ploidy)
  
  D_matrix <- matrix(0, nrow = c_ploidy, ncol = c_ploidy)
  
  for (r in 1:c_ploidy) {
    for (s in 1:c_ploidy) {
      D_rs <- 0
      for (i in 0:c_ploidy) {
        for (j in 0:c_ploidy) {
          D_rs <- D_rs + P_AB[i+1, j+1] * Phi_A[i+1, r+1] * Phi_B[j+1, s+1]
        }
      }
      D_matrix[r, s] <- D_rs
    }
  }
  
  rownames(D_matrix) <- paste0("Order_", 1:c_ploidy, "_A")
  colnames(D_matrix) <- paste0("Order_", 1:c_ploidy, "_B")
  
  return(D_matrix)
}

convert_b_to_a_4x <- function(b) { return((2 * b) / (3 - b)) }

reshape_EM_to_PAB <- function(em_vector, c_ploidy, byrow = TRUE) {
  P_AB <- matrix(em_vector, nrow = c_ploidy + 1, ncol = c_ploidy + 1, byrow = byrow)
  return(P_AB)
}


calc_Deab_comp <- function(P_AB) {
  c_ploidy <- nrow(P_AB) - 1
  QA <- rowSums(P_AB)
  QB <- colSums(P_AB)
  pA <- sum((0:c_ploidy) * QA) / c_ploidy
  pB <- sum((0:c_ploidy) * QB) / c_ploidy
  
  pAB <- sum(outer(0:c_ploidy, 0:c_ploidy, "*") * P_AB) / (c_ploidy^2)
  return(pAB - pA * pB)
}