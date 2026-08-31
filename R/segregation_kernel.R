
choose0 <- function(n, k) {
  if (!is.finite(n) || !is.finite(k)) return(0)
  if (n < 0 || k < 0 || k > n) return(0)
  return(choose(n, k))
}


transition_kernel_biallelic <- function(v, eta) {
  c_val <- v / 2
  L <- floor(v / 4)
  
  if (length(eta) != L + 1) {
    stop(paste("For ploidy v =", v, ", the segregation vector eta must be of length", L + 1))
  }
  if (abs(sum(eta) - 1.0) > 1e-8 || any(eta < 0)) {
    stop("The segregation vector eta must lie within the standard probability simplex.")
  }
  
  T_mat <- matrix(0, nrow = c_val + 1, ncol = v + 1)
  rownames(T_mat) <- paste0("i=", 0:c_val)
  colnames(T_mat) <- paste0("k=", 0:v)
  
  for (i in 0:c_val) {          
    for (k in 0:v) {            
      prob_ik <- 0
      for (l in 0:L) {          
        eta_l <- eta[l + 1]
        if (eta_l == 0) next
        
        sum_j <- 0
        for (j in 0:l) {
          num1 <- choose0(k, j)
          num2 <- choose0(k - j, i - 2 * j)
          num3 <- choose0(v - k, l - j)
          num4 <- choose0(v - k - l + j, c_val - i - 2 * l + 2 * j)
          
          den1 <- choose0(v, l)
          den2 <- choose0(v - l, c_val - 2 * l)
          
          if (den1 > 0 && den2 > 0) {
            term <- (num1 * num2 * num3 * num4) / (den1 * den2)
            sum_j <- sum_j + term
          }
        }
        prob_ik <- prob_ik + eta_l * sum_j
      }
      T_mat[i + 1, k + 1] <- prob_ik
    }
  }
  

  if (any(T_mat < -1e-10)) {
    stop("Transition kernel contains negative probabilities. Eta may be outside the simplex.")
  }
  T_mat[T_mat < 0] <- 0
  
  col_sums <- colSums(T_mat)
  if (max(abs(col_sums - 1)) > 1e-8) {
    stop("Transition kernel columns do not sum to 1. Check combinatorial formula or eta.")
  }
  T_mat <- sweep(T_mat, 2, col_sums, "/")
  
  expected_i <- colSums(T_mat * matrix(0:c_val, nrow = c_val + 1, ncol = v + 1))
  target_i <- (0:v) / 2
  if (max(abs(expected_i - target_i)) > 1e-8) {
    warning("Transition kernel does not conserve allele dosage expectation.")
  }
  
  return(T_mat)
}

convolve_gametes <- function(Q) {
  c_val <- length(Q) - 1
  v <- 2 * c_val
  Z <- numeric(v + 1)
  for(i in 0:c_val) {
    for(j in 0:c_val) {
      Z[i + j + 1] <- Z[i + j + 1] + Q[i + 1] * Q[j + 1]
    }
  }
  return(Z)
}

solve_equilibrium_Q <- function(v, p, trans_matrix, tol=1e-12, max_iter=2000) {
  c_val <- v / 2
  Q_current <- dbinom(0:c_val, size = c_val, prob = p)
  
  for (iter in 1:max_iter) {
    Z_current <- convolve_gametes(Q_current) 
    Q_next <- as.vector(trans_matrix %*% Z_current)
    
    if (any(Q_next < -1e-10)) {
      stop("Equilibrium iteration produced negative gametic probabilities.")
    }
    Q_next[abs(Q_next) < 1e-15] <- 0
    Q_next <- Q_next / sum(Q_next)
    
    if (max(abs(Q_next - Q_current)) < tol) {
      return(Q_next)
    }
    Q_current <- Q_next
  }
  warning(paste("Equilibrium solver failed to converge for p =", p))
  return(Q_current)
}


get_equilibrium_Q_one_eta <- function(v, p, eta1) {
  L <- floor(v / 4)
  eta <- numeric(L + 1)
  
  eta[1] <- 1 - eta1
  eta[2] <- eta1
  
  T_mat <- transition_kernel_biallelic(v, eta)
  Q_eq <- solve_equilibrium_Q(v, p, T_mat)
  return(Q_eq)
}


project_ESD_general <- function(Q_emp) {
  c_val <- length(Q_emp) - 1
  v <- c_val * 2
  p_emp <- sum(Q_emp * (0:c_val)) / c_val
  
  if (p_emp < 1e-4 || p_emp > 1 - 1e-4) {
    return(list(eta_ESD = NA_real_, p = p_emp, loss = 0))
  }
  
  loss_func <- function(eta1) {
    Q_theo <- get_equilibrium_Q_one_eta(v, p_emp, eta1)
    return(sum((Q_emp - Q_theo)^2))
  }
  
  opt_res <- optim(par = 0.05, fn = loss_func, method = "L-BFGS-B", lower = 0, upper = 1)
  
  return(list(eta_ESD = opt_res$par, p = p_emp, loss = opt_res$value))
}


project_ESD_4x <- function(Q_emp) project_ESD_general(Q_emp)
project_ESD_6x <- function(Q_emp) project_ESD_general(Q_emp)



# The numerical fixed-point solver is validated against the closed-form 
# autotetraploid equilibrium solution.

get_equilibrium_Q_4x_closed <- function(p, a) {
  Q <- numeric(3)
  Q[3] <- ((2 - 2*a)/(2 + a)) * p^2 + (3*a)/(2 + a) * p
  Q[2] <- 2 * ((2 - 2*a)/(2 + a)) * p * (1 - p)
  Q[1] <- ((2 - 2*a)/(2 + a)) * (1 - p)^2 + (3*a)/(2 + a) * (1 - p)
  return(Q)
}

validate_4x_solver <- function() {
  p_test <- 0.3
  eta1_test <- 1/6
  
  Q_num <- get_equilibrium_Q_one_eta(v = 4, p = p_test, eta1 = eta1_test)
  Q_closed <- get_equilibrium_Q_4x_closed(p = p_test, a = eta1_test)
  
  max_err <- max(abs(Q_num - Q_closed))
  if (max_err < 1e-10) {
    cat(sprintf("Unit Test Passed: Numeric vs Closed-form max error = %g\n", max_err))
  } else {
    warning(sprintf("Unit Test Failed: Numeric vs Closed-form max error = %g\n", max_err))
  }
}


# validate_4x_solver()
