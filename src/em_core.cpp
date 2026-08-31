
#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
List EM_full_cpp(NumericVector counts, int c_ploidy, double tol = 1e-8, int max_iter = 10000) {
  int v_ploidy = 2 * c_ploidy;
  int n_zyg = (v_ploidy + 1) * (v_ploidy + 1);
  int n_gam = (c_ploidy + 1) * (c_ploidy + 1);
  
  NumericVector pi(n_gam, 1.0 / n_gam);
  double max_diff = 1.0;
  int iter = 0;
  
  double prev_logLik = -1e30;
  bool monotone = true;
  
  while(max_diff > tol && iter < max_iter) {
    NumericVector pi_new(n_gam, 0.0);
    double total_w = 0.0; 
    double current_logLik = 0.0;
    
    for(int Z = 0; Z < n_zyg; ++Z) {
      if(counts[Z] == 0) continue;
      double w_z = counts[Z]; // 权重已撤销，全为 1
      total_w += w_z;
      
      int x = Z / (v_ploidy + 1);
      int y = Z % (v_ploidy + 1);
      double prob_Z = 0.0;
      
      for(int i1 = std::max(0, x - c_ploidy); i1 <= std::min(c_ploidy, x); ++i1) {
        int i2 = x - i1;
        for(int j1 = std::max(0, y - c_ploidy); j1 <= std::min(c_ploidy, y); ++j1) {
          int j2 = y - j1;
          prob_Z += pi[i1 * (c_ploidy + 1) + j1] * pi[i2 * (c_ploidy + 1) + j2];
        }
      }
      
      if(prob_Z > 0) {
        current_logLik += counts[Z] * std::log(prob_Z);
        
        for(int i1 = std::max(0, x - c_ploidy); i1 <= std::min(c_ploidy, x); ++i1) {
          int i2 = x - i1;
          for(int j1 = std::max(0, y - c_ploidy); j1 <= std::min(c_ploidy, y); ++j1) {
            int j2 = y - j1;
            int s = i1 * (c_ploidy + 1) + j1;
            int t = i2 * (c_ploidy + 1) + j2;
            double posterior = (pi[s] * pi[t]) / prob_Z; 
            pi_new[s] += w_z * posterior;
            pi_new[t] += w_z * posterior;
          }
        }
      }
    }
    
    // 检查似然单调性 (允许微小数值误差)
    if (iter > 0 && current_logLik < prev_logLik - 1e-8) {
      monotone = false;
    }
    prev_logLik = current_logLik;
    
    max_diff = 0.0;
    for(int s = 0; s < n_gam; ++s) {
      pi_new[s] /= (2.0 * total_w);
      double diff = std::abs(pi_new[s] - pi[s]);
      if(diff > max_diff) max_diff = diff;
      pi[s] = pi_new[s];
    }
    iter++;
  }
  
  return List::create(
    _["pi"] = pi,
    _["converged"] = (max_diff <= tol),
                       _["iterations"] = iter,
                       _["max_diff"] = max_diff,
                       _["logLik"] = prev_logLik,
                       _["likelihood_monotone"] = monotone
  );
}


// ==============================================================================
// 2. 部分信息纯净 EM 算法 (Partial / Coarsened Information EM)
// 适用于仅能观测到 0, mid, v 状态退化的数据 (9种组合)
// 自带收敛监控与似然单调性检查
// ==============================================================================
// [[Rcpp::export]]
List EM_partial_cpp(NumericVector counts, int c_ploidy, double tol = 1e-8, int max_iter = 10000) {
  int v_ploidy = 2 * c_ploidy;
  int n_gam = (c_ploidy + 1) * (c_ploidy + 1);
  NumericVector pi(n_gam, 1.0 / n_gam);
  double max_diff = 1.0;
  int iter = 0;
  
  double prev_logLik = -1e30;
  bool monotone = true;
  
  while(max_diff > tol && iter < max_iter) {
    NumericVector pi_new(n_gam, 0.0);
    double total_w = 0.0; 
    NumericVector prob_Y(9, 0.0); 
    double current_logLik = 0.0;
    
    // Step 1: E-step 前向计算粗粒化概率 prob_Y
    for(int i1 = 0; i1 <= c_ploidy; ++i1) {
      for(int j1 = 0; j1 <= c_ploidy; ++j1) {
        for(int i2 = 0; i2 <= c_ploidy; ++i2) {
          for(int j2 = 0; j2 <= c_ploidy; ++j2) {
            int s = i1 * (c_ploidy + 1) + j1;
            int t = i2 * (c_ploidy + 1) + j2;
            int x = i1 + i2; 
            int y = j1 + j2; 
            int yA = (x == 0) ? 0 : ((x == v_ploidy) ? 2 : 1);
            int yB = (y == 0) ? 0 : ((y == v_ploidy) ? 2 : 1);
            prob_Y[yA * 3 + yB] += pi[s] * pi[t];
          }
        }
      }
    }
    
    // Step 2: 累加当前对数似然 (Log-Likelihood)
    for(int Y = 0; Y < 9; ++Y) {
      if(counts[Y] > 0) {
        total_w += counts[Y]; // 纯净的样本量，无权重
        if(prob_Y[Y] > 0) {
          current_logLik += counts[Y] * std::log(prob_Y[Y]);
        }
      }
    }
    
    // Step 3: 后验概率重分配 (E-step -> M-step)
    for(int i1 = 0; i1 <= c_ploidy; ++i1) {
      for(int j1 = 0; j1 <= c_ploidy; ++j1) {
        for(int i2 = 0; i2 <= c_ploidy; ++i2) {
          for(int j2 = 0; j2 <= c_ploidy; ++j2) {
            int x = i1 + i2;
            int y = j1 + j2;
            int yA = (x == 0) ? 0 : ((x == v_ploidy) ? 2 : 1);
            int yB = (y == 0) ? 0 : ((y == v_ploidy) ? 2 : 1);
            int Y = yA * 3 + yB;
            
            if(counts[Y] > 0 && prob_Y[Y] > 0) {
              int s = i1 * (c_ploidy + 1) + j1;
              int t = i2 * (c_ploidy + 1) + j2;
              double posterior = (pi[s] * pi[t]) / prob_Y[Y];
              pi_new[s] += counts[Y] * posterior;
              pi_new[t] += counts[Y] * posterior;
            }
          }
        }
      }
    }
    
    // 收敛单调性防线 (允许极微小的机器浮点误差)
    if (iter > 0 && current_logLik < prev_logLik - 1e-8) {
      monotone = false;
    }
    prev_logLik = current_logLik;
    
    // Step 4: 更新并计算最大偏离
    max_diff = 0.0;
    for(int s = 0; s < n_gam; ++s) {
      pi_new[s] /= (2.0 * total_w);
      double diff = std::abs(pi_new[s] - pi[s]);
      if(diff > max_diff) max_diff = diff;
      pi[s] = pi_new[s];
    }
    iter++;
  }
  
  // 返回极其规范的诊断 List
  return List::create(
    _["pi"] = pi,
    _["converged"] = (max_diff <= tol),
                       _["iterations"] = iter,
                       _["max_diff"] = max_diff,
                       _["logLik"] = prev_logLik,
                       _["likelihood_monotone"] = monotone
  );
}