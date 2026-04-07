#' Utility Functions for AutoLD Package
#' @noRd
solve_a_from_QAA_4x <- function(p, QAA){
  if(p <= 0 || p >= 1) stop("p must be in (0,1)")
  pAA <- p^2; pAa <- 2*p*(1-p); paa <- (1-p)^2
  S1 <- 2*pAA*pAa; S2 <- 2*pAA*paa + pAa^2; S3 <- 2*pAa*paa
  C0 <- pAA^2 + 0.5*S1 + (1/6)*S2
  C1 <- 0.25*S1 + (1/3)*S2 + 0.25*S3
  a <- (QAA - C0) / C1
  return(a)
}

solve_a_from_T_6x <- function(p, q_obs){
  q0 <- gamete_margin_DR_6x(p, 0)
  q1 <- gamete_margin_DR_6x(p, 1)
  T_obs <- (q_obs["AAA"] + q_obs["aaa"]) - (q_obs["AAa"] + q_obs["Aaa"])
  T0    <- (q0["AAA"] + q0["aaa"]) - (q0["AAa"] + q0["Aaa"])
  T1    <- (q1["AAA"] + q1["aaa"]) - (q1["AAa"] + q1["Aaa"])
  denom <- (T1 - T0)
  if(abs(denom) < 1e-12) return(NA_real_)
  a <- (T_obs - T0) / denom
  max(0, min(1, a))
}

estimate_a_from_gamete6x_freq <- function(p, q_obs){
  q0 <- gamete_margin_DR_6x(p, a = 0)
  q1 <- gamete_margin_DR_6x(p, a = 1)
  ord <- c("AAA","AAa","Aaa","aaa")
  q_obs <- q_obs[ord]; q0 <- q0[ord]; q1 <- q1[ord]
  d <- q1 - q0
  if(all(abs(d) < 1e-12)) return(list(a = NA_real_, reason = "Not identifiable"))
  a_hat <- sum((q_obs - q0) * d) / sum(d^2)
  a_hat <- max(0, min(1, a_hat))
  list(a = a_hat, q0 = q0, q1 = q1, d = d)
}

solveD <- function(pA, pB, pv){
  HM <- H(pA=pA, pB=pB)  # 需确保 H() 函数已在你的 auto6_util.R 中定义
  HM1 <- HM[,-1]
  pv1 <- c(pv-HM[,1])
  rr <- MASS::ginv(HM1) %*% pv1
  return(as.numeric(rr))
}

gamete_margin_DR_4x <- function(p, a, D){
  p <- as.numeric(p); a <- as.numeric(a); D <- as.numeric(D)
  q0 <- c(AA = p^2, Aa = 2*p*(1-p), aa = (1-p)^2)
  qE <- c(AA = p, Aa = 0, aa = 1-p)
  q <- (1-a)*q0 + a*qE
  q_final <- c(AA = q["AA"] + D, Aa = q["Aa"] - 2*D, aa = q["aa"] + D)
  q_final[q_final < 0] <- 0
  return(q_final / sum(q_final))
}

gamete_margin_DR_4x_full <- function(p, a){
  pAA <- p^2; pAa <- 2*p*(1-p); paa <- (1-p)^2
  QAA <- pAA^2 + (3/4*a + 1/2*(1-a))*(2*pAA*pAa) + (1/2*a + 1/6*(1-a))*(2*pAA*paa + pAa^2) + 1/4*a*(2*pAa*paa)
  QAa <- 1/2*(1-a)*(2*pAA*pAa) + 2/3*(1-a)*(2*pAA*paa + pAa^2) + 1/2*(1-a)*(2*pAa*paa)
  Qaa <- 1/4*a*(2*pAA*pAa) + (1/2*a + 1/6*(1-a))*(2*pAA*paa + pAa^2) + (3/4*a + 1/2*(1-a))*(2*pAa*paa) + paa^2
  q <- c(QAA, QAa, Qaa); names(q) <- c("AA","Aa","aa")
  q
}
comDA2_B2 <- function(x,ii,p){

  xx <- x[ii]
  xx1 <- xx[1:2];xx2 <- xx[3:4];xx3 <- xx[5:6];
  signs <- (-1)^c(sum(xx1),sum(xx2),sum(xx3))
  pii <- c(p[x[-ii[1:2]]],p[x[-ii[3:4]]],p[x[-ii[5:6]]])
  sum(pii*signs*x[4]*x[5]*x[6])

}

comDab1 <- function(x,p,q){

  xx1 <- x[c(1,4)];xx2 <- x[c(2,5)];xx3 <- x[c(3,6)];
  signs <- (-1)^c(sum(xx1),sum(xx2),sum(xx3))

  pii <- c(prod(c(p[x[1:3][-1]],q[x[4:6][-1]])),prod(c(p[x[1:3][-2]],q[x[4:6][-2]])),
           prod(c(p[x[1:3][-3]],q[x[4:6][-3]])))

  sum(pii*signs)

}

comDab2 <- function(x,p,q){

  xx1 <- x[c(1,5)];xx2 <- x[c(1,6)];xx3 <- x[c(2,4)];
  xx4 <- x[c(2,6)];xx5 <- x[c(3,4)];xx6 <- x[c(3,5)];
  signs <- (-1)^c(sum(xx1),sum(xx2),sum(xx3),sum(xx4),sum(xx5),sum(xx6))

  pii <- c(prod(c(p[x[1:3][-1]],q[x[4:6][-2]])),prod(c(p[x[1:3][-1]],q[x[4:6][-3]])),
           prod(c(p[x[1:3][-2]],q[x[4:6][-1]])),prod(c(p[x[1:3][-2]],q[x[4:6][-3]])),
           prod(c(p[x[1:3][-3]],q[x[4:6][-1]])),prod(c(p[x[1:3][-3]],q[x[4:6][-2]])))

  sum(pii*signs)

}



comDAb <- function(x,p,q){

  x1 <- x[1:3];x2 <- x[4:6]
  nc1 <- rbind(c(1,2),c(1,2),c(1,2),c(1,3),c(1,3),c(1,3),c(2,3),c(2,3),c(2,3))
  nc2 <- rep(c(1,2,3),3)

  signs <- -(-1)^rowSums(cbind(matrix(x1[nc1],ncol=2,byrow=F),x2[nc2]))
  pii <- c(p[x1[3]]*prod(q[x2[c(2,3)]]),p[x1[3]]*prod(q[x2[c(1,3)]]),p[x1[3]]*prod(q[x2[c(1,2)]]),
           p[x1[2]]*prod(q[x2[c(2,3)]]),p[x1[2]]*prod(q[x2[c(1,3)]]),p[x1[2]]*prod(q[x2[c(1,2)]]),
           p[x1[1]]*prod(q[x2[c(2,3)]]),p[x1[1]]*prod(q[x2[c(1,3)]]),p[x1[1]]*prod(q[x2[c(1,2)]]))
  sum(pii*signs)

}


comDAAb <- function(x,p){

  xx1 <- x[c(1:3,4)];xx2 <- x[c(1:3,5)];xx3 <- x[c(1:3,6)];
  signs <- (-1)^c(sum(xx1),sum(xx2),sum(xx3))

  pii <- c(prod(p[x[c(5,6)]]),prod(p[x[c(4,6)]]),
           prod(p[x[c(4,5)]]))

  sum(pii*signs)

}


comDAB <- function(x,p,q){

  xx1 <- x[c(1,2,4,5)];xx2 <- x[c(1,2,4,6)];xx3 <- x[c(1,2,5,6)];
  xx4 <- x[c(1,3,4,5)];xx5 <- x[c(1,3,4,6)];xx6 <- x[c(1,3,5,6)];
  xx7 <- x[c(2,3,4,5)];xx8 <- x[c(2,3,4,6)];xx9 <- x[c(2,3,5,6)];
  signs <- (-1)^c(sum(xx1),sum(xx2),sum(xx3),sum(xx4),sum(xx5),sum(xx6),sum(xx7),sum(xx8),sum(xx9))

  pii <- c(prod(c(p[x[3]],q[x[6]])),prod(c(p[x[3]],q[x[5]])),prod(c(p[x[3]],q[x[4]])),
           prod(c(p[x[2]],q[x[6]])),prod(c(p[x[2]],q[x[5]])),prod(c(p[x[2]],q[x[4]])),
           prod(c(p[x[1]],q[x[6]])),prod(c(p[x[1]],q[x[5]])),prod(c(p[x[1]],q[x[4]])))

  sum(pii*signs)

}


comDAAB <- function(x,p){

  xx1 <- x[c(1:3,4,5)];xx2 <- x[c(1:3,4,6)];xx3 <- x[c(1:3,5,6)]
  signs <- -(-1)^c(sum(xx1),sum(xx2),sum(xx3))

  pii <- c(p[x[6]],p[x[5]],p[x[4]])

  sum(pii*signs)

}

DF4_cLDp <- function(para,hwdA,hwdB){

  gfA <- gamete_margin_DR_4x(p = para[1], a = para[3], D=hwdA)
  gfB <- gamete_margin_DR_4x(p = para[2], a = para[4], D=hwdB)

  gem <- DF4(gfA=gfA,gfB=gfB,pall=para[-c(1:4)])

  gp0 <- auto_DF4_from_gp_partial(gp=gem)
  gp0
}


auto_DF4_from_gp_partial <- function(gp){

  PAABB <- gp[1]; PAABb <- gp[2]; PAAbb <- gp[3]
  PAaBB <- gp[4]; PAaBb <- gp[5]; PAabb <- gp[6]
  PaaBB <- gp[7]; PaaBb <- gp[8]; Paabb <- gp[9]

  zy1 <- (PAABB)^2
  zy2 <- 2*PAABB*PAABb
  zy3 <- 2*PAABB*PAAbb+(PAABb)^2
  zy4 <- 2*PAABb*PAAbb
  zy5 <- (PAAbb)^2

  zy6 <- 2*PAABB*PAaBB
  zy7 <- 2*PAABB*PAaBb+2*PAaBB*PAABb
  zy8 <- 2*PAABb*PAaBb+2*PAaBB*PAAbb+2*PAabb*PAABB
  zy9 <- 2*PAAbb*PAaBb+2*PAabb*PAABb
  zy10 <- 2*PAabb*PAAbb

  zy11 <- (PAaBB)^2+2*PAABB*PaaBB
  zy12 <- 2*PAaBB*PAaBb+2*PAABB*PaaBb+2*PAABb*PaaBB
  zy13 <- (PAaBb)^2+2*PAABb*PaaBb+2*PAABB*Paabb+2*PAAbb*PaaBB+2*PAabb*PAaBB
  zy14 <- 2*PAabb*PAaBb+2*PAABb*Paabb+2*PAAbb*PaaBb
  zy15 <- (PAabb)^2+2*PAAbb*Paabb

  zy16 <- 2*PAaBB*PaaBB
  zy17 <- 2*PaaBB*PAaBb+2*PaaBb*PAaBB
  zy18 <- 2*PaaBb*PAaBb+2*Paabb*PAaBB+2*PaaBB*PAabb
  zy19 <- 2*Paabb*PAaBb+2*PaaBb*PAabb
  zy20 <- 2*Paabb*PAabb

  zy21 <- (PaaBB)^2
  zy22 <- 2*PaaBB*PaaBb
  zy23 <- (PaaBb)^2+2*pmax(PaaBB,0)*pmax(Paabb,0) # 注意这里的 pmax 是原代码保留的逻辑
  zy24 <- 2*PaaBb*Paabb
  zy25 <- (Paabb)^2

  prob_mat_25 <- matrix(c(zy1, zy2, zy3, zy4, zy5,
                          zy6, zy7, zy8, zy9, zy10,
                          zy11,zy12,zy13,zy14,zy15,
                          zy16,zy17,zy18,zy19,zy20,
                          zy21,zy22,zy23,zy24,zy25),
                        nrow=5, byrow=TRUE)


  idx_A_hom1 <- 1
  idx_A_het  <- 2:4
  idx_A_hom2 <- 5

  idx_B_hom1 <- 1
  idx_B_het  <- 2:4
  idx_B_hom2 <- 5

  p11 <- sum(prob_mat_25[idx_A_hom1, idx_B_hom1])
  p12 <- sum(prob_mat_25[idx_A_hom1, idx_B_het])
  p13 <- sum(prob_mat_25[idx_A_hom1, idx_B_hom2])

  p21 <- sum(prob_mat_25[idx_A_het, idx_B_hom1])
  p22 <- sum(prob_mat_25[idx_A_het, idx_B_het])
  p23 <- sum(prob_mat_25[idx_A_het, idx_B_hom2])

  p31 <- sum(prob_mat_25[idx_A_hom2, idx_B_hom1])
  p32 <- sum(prob_mat_25[idx_A_hom2, idx_B_het])
  p33 <- sum(prob_mat_25[idx_A_hom2, idx_B_hom2])

  tmp <- c(p11, p12, p13, p21, p22, p23, p31, p32, p33)
  return(tmp)
}


DF6_cLDp <- function(para,hwdA,hwdB){

  gfA <- gamete_margin_DR_6xD(p = para[1], a = para[3],D=hwdA)  # A位点配子边际
  gfB <- gamete_margin_DR_6xD(p = para[2], a = para[4],D=hwdB)  # B位点配子边际
  names(gfA) <- names(gfB) <- c("aaa","Aaa","AAa","AAA")
  Da <- DA2_DA3_from_gamete6x(q=gfA)
  Db <- DA2_DA3_from_gamete6x(q=gfB)

  gem <- DF6(pall=c(para[1],para[2],Da[2],Db[2],para[5],Da[3],Db[3],para[-c(1:5)]))
  gp0 <- auto_DF6_from_gp_partial(gp=gem)
  gp0
}





auto_DF6_from_gp_partial <- function(gp){

  PAAABBB <- gp[1];  PAAABBb <- gp[2];  PAAABbb <- gp[3];  PAAAbbb <- gp[4]
  PAAaBBB <- gp[5];  PAAaBBb <- gp[6];  PAAaBbb <- gp[7];  PAAabbb <- gp[8]
  PAaaBBB <- gp[9];  PAaaBBb <- gp[10]; PAaaBbb <- gp[11]; PAaabbb <- gp[12]
  PaaaBBB <- gp[13]; PaaaBBb <- gp[14]; PaaaBbb <- gp[15]; Paaabbb <- gp[16]

  # Row 1: A Dosage 6 (AAAAAA)
  zy1 <- (PAAABBB)^2
  zy2 <- 2*PAAABBB*PAAABBb
  zy3 <- 2*PAAABBB*PAAABbb+(PAAABBb)^2
  zy4 <- 2*PAAABBB*PAAAbbb+2*PAAABBb*PAAABbb
  zy5 <- 2*PAAABBb*PAAAbbb+(PAAABbb)^2
  zy6 <- 2*PAAABbb*PAAAbbb
  zy7 <- (PAAAbbb)^2

  zy8 <- 2*PAAABBB*PAAaBBB
  zy9 <- 2*PAAABBB*PAAaBBb+2*PAAaBBB*PAAABBb
  zy10 <- 2*PAAABBB*PAAaBbb+2*PAAaBBB*PAAABbb+2*PAAABBb*PAAaBBb
  zy11 <- 2*PAAABBB*PAAabbb+2*PAAaBBB*PAAAbbb+2*PAAABBb*PAAaBbb+2*PAAaBBb*PAAABbb
  zy12 <- 2*PAAABBb*PAAabbb+2*PAAaBBb*PAAAbbb+2*PAAABbb*PAAaBbb
  zy13 <- 2*PAAABbb*PAAabbb+2*PAAaBbb*PAAAbbb
  zy14 <- 2*PAAAbbb*PAAabbb

  zy15 <- 2*PAAABBB*PAaaBBB+(PAAaBBB)^2
  zy16 <- 2*PAAABBB*PAaaBBb+2*PAAaBBB*PAAaBBb+2*PAaaBBB*PAAABBb
  zy17 <- (PAAaBBb)^2+2*PAAABBB*PAaaBbb+2*PAAaBBB*PAAaBbb+2*PAaaBBB*PAAABbb+2*PAAABBb*PAaaBBb
  zy18 <- 2*PAAABBB*PAaabbb+2*PAAaBBB*PAAabbb+2*PAaaBBB*PAAAbbb+2*PAAABBb*PAaaBbb+2*PAAaBBb*PAAaBbb+2*PAaaBBb*PAAABbb
  zy19 <- (PAAaBbb)^2+2*PAAABBb*PAaabbb+2*PAAaBBb*PAAabbb+2*PAaaBBb*PAAAbbb+2*PAAABbb*PAaaBbb
  zy20 <- 2*PAAABbb*PAaabbb+2*PAaaBbb*PAAAbbb+2*PAAaBbb*PAAabbb
  zy21 <- (PAAabbb)^2+2*PAAAbbb*PAaabbb

  zy22 <- 2*PAAABBB*PaaaBBB+2*PAAaBBB*PAaaBBB
  zy23 <- 2*PAAABBB*PaaaBBb+2*PAAaBBB*PAaaBBb+2*PAaaBBB*PAAaBBb+2*PaaaBBB*PAAABBb
  zy24 <- 2*PAAABBB*PaaaBbb+2*PAAaBBB*PAaaBbb+2*PAaaBBB*PAAaBbb+2*PaaaBBB*PAAABbb+2*PAAABBb*PaaaBBb+2*PAAaBBb*PAaaBBb
  zy25 <- 2*PAAABBB*Paaabbb+2*PAAaBBB*PAaabbb+2*PAaaBBB*PAAabbb+2*PaaaBBB*PAAAbbb+2*PAAABBb*PaaaBbb+2*PAAaBBb*PAaaBbb+2*PAaaBBb*PAAaBbb+2*PaaaBBb*PAAABbb
  zy26 <- 2*PAAABBb*Paaabbb+2*PAAaBBb*PAaabbb+2*PAaaBBb*PAAabbb+2*PaaaBBb*PAAAbbb+2*PAAABbb*PaaaBbb+2*PAAaBbb*PAaaBbb
  zy27 <- 2*PAAABbb*Paaabbb+2*PAAaBbb*PAaabbb+2*PAaaBbb*PAAabbb+2*PaaaBbb*PAAAbbb
  zy28 <- 2*PAAAbbb*Paaabbb+2*PAAabbb*PAaabbb


  zy29 <- (PAaaBBB)^2+2*PAAaBBB*PaaaBBB
  zy30 <- 2*PAAaBBB*PaaaBBb+2*PAaaBBB*PAaaBBb+2*PaaaBBB*PAAaBBb
  zy31 <- (PAaaBBb)^2+2*PAAaBBB*PaaaBbb+2*PAaaBBB*PAaaBbb+2*PaaaBBB*PAAaBbb+2*PAAaBBb*PaaaBBb
  zy32 <- 2*PAAaBBB*Paaabbb+2*PAaaBBB*PAaabbb+2*PaaaBBB*PAAabbb+2*PAAaBBb*PaaaBbb+2*PAaaBBb*PAaaBbb+2*PaaaBBb*PAAaBbb
  zy33 <- (PAaaBbb)^2+2*PAAaBBb*Paaabbb+2*PAaaBBb*PAaabbb+2*PaaaBBb*PAAabbb+2*PAAaBbb*PaaaBbb
  zy34 <- 2*PAAaBbb*Paaabbb+2*PAaaBbb*PAaabbb+2*PaaaBbb*PAAabbb
  zy35 <- (PAaabbb)^2+2*PAAabbb*Paaabbb


  zy36 <- 2*PAaaBBB*PaaaBBB
  zy37 <- 2*PAaaBBB*PaaaBBb+2*PaaaBBB*PAaaBBb
  zy38 <- 2*PAaaBBB*PaaaBbb+2*PaaaBBB*PAaaBbb+2*PAaaBBb*PaaaBBb
  zy39 <- 2*PAaaBBB*Paaabbb+2*PaaaBBB*PAaabbb+2*PAaaBBb*PaaaBbb+2*PaaaBBb*PAaaBbb
  zy40 <- 2*PAaaBBb*Paaabbb+2*PaaaBBb*PAaabbb+2*PAaaBbb*PaaaBbb
  zy41 <- 2*PAaaBbb*Paaabbb+2*PaaaBbb*PAaabbb
  zy42 <- 2*PAaabbb*Paaabbb

  zy43 <- (PaaaBBB)^2
  zy44 <- 2*PaaaBBB*PaaaBBb
  zy45 <- (PaaaBBb)^2+2*PaaaBBB*PaaaBbb
  zy46 <- 2*PaaaBBB*Paaabbb+2*PaaaBBb*PaaaBbb
  zy47 <- (PaaaBbb)^2+2*PaaaBBb*Paaabbb
  zy48 <- 2*PaaaBbb*Paaabbb
  zy49 <- (Paaabbb)^2


  prob_mat_49 <- matrix(c(zy1,zy2,zy3,zy4,zy5,zy6,zy7,
                          zy8,zy9,zy10,zy11,zy12,zy13,zy14,
                          zy15,zy16,zy17,zy18,zy19,zy20,zy21,
                          zy22,zy23,zy24,zy25,zy26,zy27,zy28,
                          zy29,zy30,zy31,zy32,zy33,zy34,zy35,
                          zy36,zy37,zy38,zy39,zy40,zy41,zy42,
                          zy43,zy44,zy45,zy46,zy47,zy48,zy49),
                        nrow=7, byrow=TRUE)


  idx_hom1 <- 1
  idx_het  <- 2:6
  idx_hom2 <- 7


  # Row 1: A is Hom1 (Dosage 6)
  p11 <- sum(prob_mat_49[idx_hom1, idx_hom1]) # A=Hom1, B=Hom1
  p12 <- sum(prob_mat_49[idx_hom1, idx_het])  # A=Hom1, B=Het
  p13 <- sum(prob_mat_49[idx_hom1, idx_hom2]) # A=Hom1, B=Hom2

  # Row 2: A is Het (Dosage 1-5)
  p21 <- sum(prob_mat_49[idx_het, idx_hom1])  # A=Het, B=Hom1
  p22 <- sum(prob_mat_49[idx_het, idx_het])   # A=Het, B=Het
  p23 <- sum(prob_mat_49[idx_het, idx_hom2])  # A=Het, B=Hom2

  # Row 3: A is Hom2 (Dosage 0)
  p31 <- sum(prob_mat_49[idx_hom2, idx_hom1]) # A=Hom2, B=Hom1
  p32 <- sum(prob_mat_49[idx_hom2, idx_het])  # A=Hom2, B=Het
  p33 <- sum(prob_mat_49[idx_hom2, idx_hom2]) # A=Hom2, B=Hom2

  tmp <- c(p11, p12, p13, p21, p22, p23, p31, p32, p33)

  return(tmp)
}
