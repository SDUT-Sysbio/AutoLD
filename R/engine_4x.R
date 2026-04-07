# === Engine functions for Autotetraploid LD (Full & Partial) ===

auto4_est_LD <- function(dose){
  nts <- rep(0,25)
  nt <- table(dose)
  nts[as.numeric(names(nt))] <- as.numeric(nt)
  nn <- sum(nts)
  res4 <- EM_C4(nts=nts)

  pAA <- sum(res4[1:3]); pAa <- sum(res4[4:6]); paa <- sum(res4[7:9])
  pBB <- sum(res4[c(1,4,7)]); pBb <- sum(res4[c(2,5,8)]); pbb <- sum(res4[c(3,6,9)])
  pA <- pAA+pAa/2; pa <- paa+pAa/2; pB <- pBB+pBb/2; pb <- pbb+pBb/2
  DA <- pAA-pA^2; DB <- pBB-pB^2

  pAB <- (2*res4[1] + sum(res4[c(2,4)])+res4[5]/2)/2
  pAb <- (2*res4[3] + sum(res4[c(2,6)])+res4[5]/2)/2
  paB <- (2*res4[7] + sum(res4[c(4,8)])+res4[5]/2)/2
  pab <- (2*res4[9] + sum(res4[c(6,8)])+res4[5]/2)/2
  Deab <- 2*(pAB-pA*pB)
  DAb <- res4[1]+res4[2]/2-pA*Deab-pB*DA-pA^2*pB
  DaB <- res4[1]+res4[4]/2-pB*Deab-pA*DB-pA*pB^2
  DAB <- res4[1]-(pA^2*pB^2+pB^2*DA+pA^2*DB+2*pB*DAb+2*pA*DaB+2*pA*pB*Deab)-DA*DB-Deab^2

  estp <- as.numeric(c(pA,pB,Deab,DAb,DaB,DAB,pAA,pBB))
  return(estp)
}

auto4_est_LDp <- function(dose){
  nts <- rep(0,9)
  nt <- table(dose)
  nts[as.numeric(names(nt))] <- as.numeric(nt)
  nn <- sum(nts)
  res4 <- EM_C4P(nts=nts)

  pAA <- sum(res4[1:3]); pAa <- sum(res4[4:6]); paa <- sum(res4[7:9])
  pBB <- sum(res4[c(1,4,7)]); pBb <- sum(res4[c(2,5,8)]); pbb <- sum(res4[c(3,6,9)])
  pA <- pAA+pAa/2; pa <- paa+pAa/2; pB <- pBB+pBb/2; pb <- pbb+pBb/2
  DA <- pAA-pA^2; DB <- pBB-pB^2

  pAB <- (2*res4[1] + sum(res4[c(2,4)])+res4[5]/2)/2
  pAb <- (2*res4[3] + sum(res4[c(2,6)])+res4[5]/2)/2
  paB <- (2*res4[7] + sum(res4[c(4,8)])+res4[5]/2)/2
  pab <- (2*res4[9] + sum(res4[c(6,8)])+res4[5]/2)/2
  Deab <- 2*(pAB-pA*pB)
  DAb <- res4[1]+res4[2]/2-pA*Deab-pB*DA-pA^2*pB
  DaB <- res4[1]+res4[4]/2-pB*Deab-pA*DB-pA*pB^2
  DAB <- res4[1]-(pA^2*pB^2+pB^2*DA+pA^2*DB+2*pB*DAb+2*pA*DaB+2*pA*pB*Deab)-DA*DB-Deab^2

  estp <- as.numeric(c(pA,pB,Deab,DAb,DaB,DAB,pAA,pBB))
  return(estp)
}

EM_C4 <- function(nts){
  p4 <- rep(1/9,9)
  nn <- sum(nts)
  iter <- 0
  nts11 <- matrix(rep(nts,9),nrow=9,byrow=T)
  while(1){
    p41 <- p4
    phi <- CEM_E4(p=p4)
    emm <- CEM_M4(phi=phi)
    p4 <- rowSums(emm*nts11)/(2*nn)
    iter <- iter + 1
    if(max(abs(p4-p41))<1e-5) break
  }
  return(p4)
}

EM_C4P <- function(nts){
  p4 <- rep(1/9,9)
  nn <- sum(nts)
  iter <- 0
  nts11 <- matrix(rep(nts,9),nrow=9,byrow=T)
  while(1){
    p41 <- p4
    emm <- CEM4P(p=p4)
    p4 <- rowSums(emm*nts11)/(2*nn)
    iter <- iter + 1
    if(max(abs(p4-p41))<1e-5) break
  }
  return(p4)
}

CEM_E4 <- function(p){
  PAABB = p[1]; PAABb = p[2]; PAAbb = p[3]
  PAaBB = p[4]; PAaBb = p[5]; PAabb = p[6]
  PaaBB = p[7]; PaaBb = p[8]; Paabb = p[9]

  phi_1 <- (2*PAABB*PAAbb)/(2*PAABB*PAAbb+(PAABb)^2)
  phi_2 <- (2*PAABB*PAaBb)/(2*PAABB*PAaBb+2*PAaBB*PAABb)
  phi_31 <- (2*PAabb*PAABB)/(2*PAABb*PAaBb+2*PAaBB*PAAbb+2*PAabb*PAABB)
  phi_32 <- (2*PAABb*PAaBb)/(2*PAABb*PAaBb+2*PAaBB*PAAbb+2*PAabb*PAABB)
  phi_4 <- (2*PAABB*PaaBB)/((PAaBB)^2+2*PAABB*PaaBB)
  phi_51 <- (2*PAABB*PaaBb)/(2*PAaBB*PAaBb+2*PAABB*PaaBb+2*PAABb*PaaBB)
  phi_52 <- (2*PAABb*PaaBB)/(2*PAaBB*PAaBb+2*PAABB*PaaBb+2*PAABb*PaaBB)
  phi_61 <- (2*PAABB*Paabb)/((PAaBb)^2+2*PAABb*PaaBb+2*PAABB*Paabb+2*PAAbb*PaaBB+2*PAabb*PAaBB)
  phi_62 <- (2*PAABb*PaaBb)/((PAaBb)^2+2*PAABb*PaaBb+2*PAABB*Paabb+2*PAAbb*PaaBB+2*PAabb*PAaBB)
  phi_63 <- (2*PAAbb*PaaBB)/((PAaBb)^2+2*PAABb*PaaBb+2*PAABB*Paabb+2*PAAbb*PaaBB+2*PAabb*PAaBB)
  phi_64 <- (2*PAabb*PAaBB)/((PAaBb)^2+2*PAABb*PaaBb+2*PAABB*Paabb+2*PAAbb*PaaBB+2*PAabb*PAaBB)
  phi_7 <- (2*PAabb*PAABb)/(2*PAAbb*PAaBb+2*PAabb*PAABb)
  phi_81 <- (2*PAABb*Paabb)/(2*PAabb*PAaBb+2*PAABb*Paabb+2*PAAbb*PaaBb)
  phi_82 <- (2*PAAbb*PaaBb)/(2*PAabb*PAaBb+2*PAABb*Paabb+2*PAAbb*PaaBb)
  phi_9 <- (2*PAAbb*Paabb)/((PAabb)^2+2*PAAbb*Paabb)
  phi_10 <- (2*PaaBb*PAaBB)/(2*PaaBB*PAaBb+2*PaaBb*PAaBB)
  phi_111 <- (2*Paabb*PAaBB)/(2*PaaBb*PAaBb+2*Paabb*PAaBB+2*PaaBB*PAabb)
  phi_112 <- (2*PaaBb*PAaBb)/(2*PaaBb*PAaBb+2*Paabb*PAaBB+2*PaaBB*PAabb)
  phi_12 <- (2*Paabb*PAaBb)/(2*Paabb*PAaBb+2*PaaBb*PAabb)
  phi_13 <- (2*PaaBB*Paabb)/((PaaBb)^2+2*PaaBB*Paabb)

  return(c(phi_1,phi_2,phi_31,phi_32,phi_4,phi_51,phi_52,phi_61,phi_62,phi_63,phi_64,phi_7,phi_81,phi_82,phi_9,phi_10,phi_111,phi_112,phi_12,phi_13))
}

CEM_M4 <- function(phi){
  phi1 <- phi[1]; phi2 <- phi[2]; phi31 <- phi[3]; phi32 <- phi[4]
  phi4 <- phi[5]; phi51 <- phi[6]; phi52 <- phi[7]; phi61 <- phi[8]
  phi62 <- phi[9]; phi63 <- phi[10]; phi64 <- phi[11]; phi7 <- phi[12]
  phi81 <- phi[13]; phi82 <- phi[14]; phi9 <- phi[15]; phi10 <- phi[16]
  phi111 <- phi[17]; phi112 <- phi[18]; phi12 <- phi[19]; phi13 <- phi[20]

  epAABB <- c(2,1,phi1,0,0,1,phi2,phi31,0,0,phi4,phi51,phi61,0,0,0,0,0,0,0,0,0,0,0,0)
  epAABb <- c(0,1,2*(1-phi1),1,0,0,1-phi2,phi32,phi7,0,0,phi52,phi62,phi81,0,0,0,0,0,0,0,0,0,0,0)
  epAAbb <- c(0,0,phi1,1,2,0,0,1-phi31-phi32,1-phi7,1,0,0,phi63,phi82,phi9,0,0,0,0,0,0,0,0,0,0)
  epAaBB <- c(0,0,0,0,0,1,1-phi2,1-phi31-phi32,0,0,2*(1-phi4),1-phi51-phi52,phi64,0,0,1,phi10,phi111,0,0,0,0,0,0,0)
  epAaBb <- c(0,0,0,0,0,0,phi2,phi32,1-phi7,0,0,1-phi51-phi52,2*(1-phi61-phi62-phi63-phi64),1-phi81-phi82,0,0,1-phi10,phi112,phi12,0,0,0,0,0,0)
  epAabb <- c(0,0,0,0,0,0,0,phi31,phi7,1,0,0,phi64,1-phi81-phi82,2*(1-phi9),0,0,1-phi111-phi112,1-phi12,1,0,0,0,0,0)
  epaaBB <- c(0,0,0,0,0,0,0,0,0,0,phi4,phi52,phi63,0,0,1,1-phi10,1-phi111-phi112,0,0,2,1,phi13,0,0)
  epaaBb <- c(0,0,0,0,0,0,0,0,0,0,0,phi51,phi62,phi82,0,0,phi10,phi112,1-phi12,0,0,1,2*(1-phi13),1,0)
  epaabb <- c(0,0,0,0,0,0,0,0,0,0,0,0,phi61,phi81,phi9,0,0,phi111,phi12,1,0,0,phi13,1,2)
  return(rbind(epAABB,epAABb,epAAbb,epAaBB,epAaBb,epAabb,epaaBB,epaaBb,epaabb))
}

CEM4P <- function(p){
  PAABB=p[1]; PAABb=p[2]; PAAbb=p[3]; PAaBB=p[4]; PAaBb=p[5]; PAabb=p[6]; PaaBB=p[7]; PaaBb=p[8]; Paabb=p[9]

  T1 <- c(2*PAABB*PAABb,2*PAABB*PAAbb,(PAABb)^2,2*PAABb*PAAbb)
  phi_11 <- (2*PAABB*PAABb)/sum(T1); phi_12 <- (2*PAABB*PAAbb)/sum(T1); phi_13 <- ((PAABb)^2)/sum(T1)

  T2 <- c(2*PAABB*PAaBB , (PAaBB)^2,2*PAABB*PaaBB , 2*PAaBB*PaaBB)
  phi_21 <- (2*PAABB*PAaBB)/sum(T2); phi_22 <- ((PAaBB)^2)/sum(T2); phi_23 <- (2*PAABB*PaaBB)/sum(T2)

  T3 <- c(2*PAABB*PAaBb,2*PAaBB*PAABb,2*PAABb*PAaBb,2*PAaBB*PAAbb,2*PAabb*PAABB,2*PAAbb*PAaBb,2*PAabb*PAABb,2*PAaBB*PAaBb,2*PAABB*PaaBb,2*PAABb*PaaBB,(PAaBb)^2,2*PAABb*PaaBb,2*PAABB*Paabb,2*PAAbb*PaaBB,2*PAabb*PAaBB,2*PAabb*PAaBb,2*PAABb*Paabb,2*PAAbb*PaaBb,2*PaaBB*PAaBb,2*PaaBb*PAaBB,2*PaaBb*PAaBb,2*Paabb*PAaBB,2*PaaBB*PAabb,2*Paabb*PAaBb,2*PaaBb*PAabb)
  phi3 <- T3/sum(T3)

  phi3_AABB <- sum(phi3[c(1,5,9,13)]); phi3_AABb <- sum(phi3[c(2,3,7,10,12,17)]); phi3_AAbb <- sum(phi3[c(4,6,14,18)])
  phi3_AaBB <- sum(phi3[c(2,4,8,15,20,22)]); phi3_AaBb <- sum(phi3[c(1,3,6,8,16,19,21,24)])+2*phi3[11]; phi3_Aabb <- sum(phi3[c(5,7,15,16,23,25)])
  phi3_aaBB <- sum(phi3[c(10,14,19,23)]); phi3_aaBb <- sum(phi3[c(9,12,18,20,21,25)]); phi3_aabb <- sum(phi3[c(13,17,22,24)])

  T4 <- c(2*PAabb*PAAbb,(PAabb)^2,2*PAAbb*Paabb,2*Paabb*PAabb)
  phi_41 <- (2*PAabb*PAAbb)/sum(T4); phi_42 <- ((PAabb)^2)/sum(T4); phi_43 <- (2*PAAbb*Paabb)/sum(T4)

  T5 <- c(2*PaaBB*PaaBb + (PaaBb)^2+2*PaaBB*Paabb + 2*PaaBb*Paabb)
  phi_51 <- (2*PaaBB*PaaBb)/sum(T5); phi_52 <- ((PaaBb)^2)/sum(T5); phi_53 <- (2*PaaBB*Paabb)/sum(T5)

  epAABB <- c(2,phi_11+phi_12,0,phi_21+phi_23,phi3_AABB,0,0,0,0)
  epAABb <- c(0,1-phi_12+phi_13,0,0,phi3_AABb,0,0,0,0)
  epAAbb <- c(0,1-phi_11-phi_13,2,0,phi3_AAbb,phi_41+phi_43,0,0,0)
  epAaBB <- c(0,0,0,1+phi_22-phi_23,phi3_AaBB,0,0,0,0)
  epAaBb <- c(0,0,0,0,phi3_AaBb,0,0,0,0)
  epAabb <- c(0,0,0,0,phi3_Aabb,1+phi_42-phi_43,0,0,0)
  epaaBB <- c(0,0,0,1-phi_21-phi_22,phi3_aaBB,0,2,phi_51+phi_53,0)
  epaaBb <- c(0,0,0,0,phi3_aaBb,0,0,1+phi_52-phi_53,0)
  epaabb <- c(0,0,0,0,phi3_aabb,1-phi_41-phi_42,0,1-phi_51-phi_52,2)
  return(rbind(epAABB,epAABb,epAAbb,epAaBB,epAaBb,epAabb,epaaBB,epaaBb,epaabb))
}

run_EM_engine_hap <- function(nts, mode) {
  p4 <- rep(0.25, 4)
  nn <- sum(nts)
  iter <- 0
  while(iter < 1000) {
    p4_prev <- p4
    if(mode == "full") { emm <- get_EM_matrix_25_hap(p4) } else { emm <- get_EM_matrix_9_hap(p4) }
    p1 <- sum(emm[1,] * nts) / (4 * nn)
    p2 <- sum(emm[2,] * nts) / (4 * nn)
    p3 <- sum(emm[3,] * nts) / (4 * nn)
    p4_new <- sum(emm[4,] * nts) / (4 * nn)
    p4 <- c(p1, p2, p3, p4_new)
    iter <- iter + 1
    if(max(abs(p4 - p4_prev)) < 1e-6) break
  }
  return(p4)
}

get_EM_matrix_25_hap <- function(p) {
  PAB <- p[1]; PAb <- p[2]; PaB <- p[3]; Pab <- p[4]
  term1 <- 4*PAB^3*Pab + 12*PAB^2*PAb*PaB; phi_1 <- ifelse(term1==0, 0, (4*PAB^3*Pab)/term1)
  term2 <- 12*PAB^2*PAb*Pab + 12*PAB*PAb^2*PaB; phi_2 <- ifelse(term2==0, 0, (12*PAB^2*PAb*Pab)/term2)
  term3 <- 12*PAB*PAb^2*Pab + 4*PAb^3*PaB; phi_3 <- ifelse(term3==0, 0, (12*PAB*PAb^2*Pab)/term3)
  term4 <- 12*PAB^2*PaB*Pab + 12*PAB*PAb*PaB^2; phi_4 <- ifelse(term4==0, 0, (12*PAB^2*PaB*Pab)/term4)
  term5 <- 6*PAB^2*Pab^2 + 24*PAB*PAb*PaB*Pab + 6*PAb^2*PaB^2
  phi_51 <- ifelse(term5==0, 0, (6*PAB^2*Pab^2)/term5); phi_52 <- ifelse(term5==0, 0, (6*PAb^2*PaB^2)/term5)
  term6 <- 12*PAB*PAb*Pab^2 + 12*PAb^2*PaB*Pab; phi_6 <- ifelse(term6==0, 0, (12*PAB*PAb*Pab^2)/term6)
  term7 <- 12*PAB*PaB^2*Pab + 4*PAb*PaB^3; phi_7 <- ifelse(term7==0, 0, (12*PAB*PaB^2*Pab)/term7)
  term8 <- 12*PAB*PaB*Pab^2 + 12*PAb*PaB^2*Pab; phi_8 <- ifelse(term8==0, 0, (12*PAB*PaB*Pab^2)/term8)
  term9 <- 4*PAB*Pab^3 + 12*PAb*PaB*Pab^2; phi_9 <- ifelse(term9==0, 0, (4*PAB*Pab^3)/term9)

  epAB <- c(4,3,2,1,0, 3,phi_1+2,1+phi_2,phi_3,0, 2,1+phi_4,phi_51-phi_52+1,phi_6,0, 1,phi_7,phi_8,phi_9,0, 0,0,0,0,0)
  epAb <- c(0,1,2,3,4, 0,1-phi_1,2-phi_2,3-phi_3,3, 0,1-phi_4,phi_52-phi_51+1,2-phi_6,2, 0,1-phi_7,1-phi_8,1-phi_9,1, 0,0,0,0,0)
  epaB <- c(0,0,0,0,0, 1,1-phi_1,1-phi_2,1-phi_3,0, 2,2-phi_4,phi_52-phi_51+1,1-phi_6,0, 3,3-phi_7,2-phi_8,1-phi_9,0, 4,3,2,1,0)
  epab <- c(0,0,0,0,0, 0,phi_1,phi_2,phi_3,1,     0,phi_4,phi_51-phi_52+1,1+phi_6,2, 0,phi_7,1+phi_8,2+phi_9,3, 0,1,2,3,4)
  return(rbind(epAB, epAb, epaB, epab))
}

get_EM_matrix_9_hap <- function(p) {
  PAB <- p[1]; PAb <- p[2]; PaB <- p[3]; Pab <- p[4]
  T1_sum <- 4*PAB^3*PAb + 6*PAB^2*PAb^2 + 4*PAB*PAb^3
  phi_11 <- ifelse(T1_sum==0,0, (4*PAB^3*PAb)/T1_sum); phi_12 <- ifelse(T1_sum==0,0, (6*PAB^2*PAb^2)/T1_sum)
  T2_sum <- 4*PAB^3*PaB + 6*PAB^2*PaB^2 + 4*PAB*PaB^3
  phi_21 <- ifelse(T2_sum==0,0, (4*PAB^3*PaB)/T2_sum); phi_22 <- ifelse(T2_sum==0,0, (6*PAB^2*PaB^2)/T2_sum)
  T4_sum <- 4*PAb^3*Pab + 6*PAb^2*Pab^2 + 4*PAb*Pab^3
  phi_41 <- ifelse(T4_sum==0,0, (4*PAb^3*Pab)/T4_sum); phi_42 <- ifelse(T4_sum==0,0, (6*PAb^2*Pab^2)/T4_sum)
  T5_sum <- 4*PaB^3*Pab + 6*PaB^2*Pab^2 + 4*PaB*Pab^3
  phi_51 <- ifelse(T5_sum==0,0, (4*PaB^3*Pab)/T5_sum); phi_52 <- ifelse(T5_sum==0,0, (6*PaB^2*Pab^2)/T5_sum)

  T3 <- c(4*PAB^3*Pab, 12*PAB^2*PAb*PaB, 12*PAB^2*PAb*Pab, 12*PAB*PAb^2*PaB, 12*PAB*PAb^2*Pab, 4*PAb^3*PaB, 12*PAB^2*PaB*Pab, 12*PAB*PAb*PaB^2, 6*PAB^2*Pab^2, 24*PAB*PAb*PaB*Pab, 6*PAb^2*PaB^2, 12*PAB*PAb*Pab^2, 12*PAb^2*PaB*Pab, 12*PAB*PaB^2*Pab, 4*PAb*PaB^3, 12*PAB*PaB*Pab^2, 12*PAb*PaB^2*Pab, 4*PAB*Pab^3, 12*PAb*PaB*Pab^2)
  sT3 <- sum(T3); if(sT3==0) sT3<-1e-9; phi3 <- T3/sT3

  i11 <- c(3,2,2,1,1,2,1,2,1,1,1,1,1); i12 <- c(1,2,3,4,5,7,8,9,10,12,14,16,18); phi3_AB <- sum(i11*phi3[i12])
  i21 <- c(1,1,2,2,3,1,1,2,1,2,1,1,1); i22 <- c(2,3,4,5,6,8,10,11,12,13,15,17,19); phi3_Ab <- sum(i21*phi3[i22])
  i31 <- c(1,1,1,1,2,1,2,1,2,3,1,2,1); i32 <- c(2,4,6,7,8,10,11,13,14,15,16,17,19); phi3_aB <- sum(i31*phi3[i32])
  i41 <- c(1,1,1,1,2,1,2,1,1,2,1,3,2); i42 <- c(1,3,5,7,9,10,12,13,14,16,17,18,19); phi3_ab <- sum(i41*phi3[i42])

  epAB <- c(4, 1+2*phi_11+phi_12, 0, 1+2*phi_21+phi_22, phi3_AB, 0, 0, 0, 0)
  epAb <- c(0, 3-2*phi_11-phi_12, 4, 0, phi3_Ab, 1+2*phi_41+phi_42, 0, 0, 0)
  epaB <- c(0, 0, 0, 3-2*phi_21-phi_22, phi3_aB, 0, 4, 1+2*phi_51+phi_52, 0)
  epab <- c(0, 0, 0, 0, phi3_ab, 3-2*phi_41-phi_42, 0, 3-2*phi_51-phi_52, 4)
  return(rbind(epAB, epAb, epaB, epab))
}

EME_p <- function(p,a){
  mt <- MM_p(p=p,a=a); n <- length(mt); EL <- c()
  for(i in 1:n){
    tmp <- mt[[i]]; QT <- tmp[,1]*tmp[,2]; QN <- tmp[,1]*tmp[,3]; QA <- tmp[,2]*tmp[,4]
    EQ <- colSums(cbind(QN,QA))/sum(QT)
    EL <- rbind(EL,EQ)
  }
  return(EL)
}

EMM_p1 <- function(ELL_p,nl){
  nl <- as.matrix(nl); nlt <- sum(nl)
  p_e <- sum(nl*ELL_p[,1])/(8*nlt)
  return(c(p_e))
}

LL2  <- function(x,p,nl){
  a <- x; pAA <- p^2; pAa <- 2*p*(1-p); paa <- (1-p)^2
  QAA <- pAA^2 + (3/4*a+1/2*(1-a))*(2*pAA*pAa) + (1/2*a+1/6*(1-a))*(2*pAA*paa + pAa^2) + 1/4*a*(2*pAa*paa)
  QAa <- 1/2*(1-a)*(2*pAA*pAa) + 2/3*(1-a)*(2*pAA*paa + pAa^2) + 1/2*(1-a)*(2*pAa*paa)
  Qaa <- 1/4*a*(2*pAA*pAa) + (1/2*a+1/6*(1-a))*(2*pAA*paa + pAa^2) + (3/4*a+1/2*(1-a))*(2*pAa*paa) + paa^2
  QAAAA <- QAA^2; QAAAa <- 2*QAA*QAa; QAAaa <- 2*QAA*Qaa + QAa^2; QAaaa <- 2*QAa*Qaa; Qaaaa <- Qaa^2
  LL <- (nl[1]*log(QAAAA)+nl[2]*log(QAAAa)+nl[3]*log(QAAaa)+nl[4]*log(QAaaa)+nl[5]*log(Qaaaa))
  -LL
}

AAAAM_p <- function(p,a){ matrix(c(1,(1/2+1/4*a)^2,(1/6+1/3*a)^2,1/16*a^2,1+1/2*a,1/3+2/3*a,1/2*a,(1+1/2*a)*(1/6+1/3*a),1/2*a*(1/6+1/3*a),1/4*a*(1+1/2*a), p^8,16*p^6*(1-p)^2,36*p^4*(1-p)^4,16*p^2*(1-p)^6,4*p^7*(1-p),6*p^6*(1-p)^2,4*p^5*(1-p)^3,24*p^5*(1-p)^3,24*p^3*(1-p)^5,16*p^4*(1-p)^4, 8*p^8,96*p^6*(1-p)^2,144*p^4*(1-p)^4,32*p^2*(1-p)^6,28*p^7*(1-p),36*p^6*(1-p)^2,20*p^5*(1-p)^3,120*p^5*(1-p)^3,72*p^3*(1-p)^5,64*p^4*(1-p)^4, 0,3/8*a*(2+a),1/3*a*(1/2+a),1/8*a^2,3/2*a,a,1/2*a,3/4*a*(1+a),1/12*a*(1+5*a),1/4*a*(1+2*a)), nrow=10) }
AAAaM_p <- function(p,a){ matrix(c(1-a,(1-a)*(1/2+1/4*a),(1-a)*(5/6+2/3*a),4/3*(1-a),2/9*(1-a)*(1+2*a),1-a,(1-a)*(1/2+1/2*a),(1-a)*(1/6+2/3*a),1/4*a*(1-a), 4*p^7*(1-p),16*p^6*(1-p)^2,24*p^5*(1-p)^3,6*p^6*(1-p)^2,36*p^4*(1-p)^4,4*p^5*(1-p)^3,16*p^4*(1-p)^4,24*p^3*(1-p)^5,16*p^2*(1-p)^6, 28*p^7*(1-p),96*p^6*(1-p)^2,120*p^5*(1-p)^3,36*p^6*(1-p)^2,144*p^4*(1-p)^4,20*p^5*(1-p)^3,64*p^4*(1-p)^4,72*p^3*(1-p)^5,32*p^2*(1-p)^6, 0,3/4*a*(1-a),3/2*a*(1-a),0,2/3*a*(1-a),0,a*(1-a),5/6*a*(1-a),1/4*a*(1-a)), nrow=9) }
AAaaM_p <- function(p,a){ matrix(c(1/2*a,1/3*(1+2*a),1+1/2*a,1/4*(1-a+3/2*a^2),5/6*(1-a+6/5*a^2),3/4*a^2-a/2+1,2,2/3*a^2-2/3*a+1/2,5/6*(1-a+6/5*a^2),1/4*(1-a+3/2*a^2),1+a/2,1/3*(1+2*a),1/2*a, 4*p^7*(1-p),6*p^6*(1-p)^2,4*p^5*(1-p)^3,16*p^6*(1-p)^2,24*p^5*(1-p)^3,16*p^4*(1-p)^4,p^4*(1-p)^4,36*p^4*(1-p)^4,24*p^3*(1-p)^5,16*p^2*(1-p)^6,4*p^3*(1-p)^5,6*p^2*(1-p)^6,4*p*(1-p)^7, 28*p^7*(1-p),36*p^6*(1-p)^2,20*p^5*(1-p)^3,96*p^6*(1-p)^2,120*p^5*(1-p)^3,64*p^4*(1-p)^4,4*p^4*(1-p)^4,144*p^4*(1-p)^4,72*p^3*(1-p)^5,32*p^2*(1-p)^6,12*p^3*(1-p)^5,12*p^2*(1-p)^6,4*p*(1-p)^7, 1/2*a,a,3/2*a,1/2*a*(1/2+a),1/6*a*(5+7*a),1/2*a*(3+2*a),0,1/3*a*(1+2*a),1/6*a*(5+7*a),1/2*a*(1/2+a),3/2*a,a,1/2*a), nrow=13) }
AaaaM_p <- function(p,a){ matrix(c(1/4*a*(1-a),(1-a)*(1/6+2/3*a),1/2*(1-a)*(1+a),1-a,4/9*(1-a)*(1/2+a),4/3*(1-a),1/3*(1-a)*(5/2+2*a),1/2*(1-a)*(1+1/2*a),1-a, 16*p^6*(1-p)^2,24*p^5*(1-p)^3,16*p^4*(1-p)^4,4*p^3*(1-p)^5,36*p^4*(1-p)^4,6*p^2*(1-p)^6,24*p^3*(1-p)^5,16*p^2*(1-p)^6,4*p*(1-p)^7, 96*p^6*(1-p)^2,120*p^5*(1-p)^3,64*p^4*(1-p)^4,12*p^3*(1-p)^5,144*p^4*(1-p)^4,12*p^2*(1-p)^6,72*p^3*(1-p)^5,32*p^2*(1-p)^6,4*p*(1-p)^7, 1/4*a*(1-a),5/6*a*(1-a),a*(1-a),0,2/3*a*(1-a),0,3/2*a*(1-a),3/4*a*(1-a),0), nrow=9) }
aaaaM_p <- function(p,a){ matrix(c(1/4*a*(1+1/2*a),1/2*a*(1/6+1/3*a),(1+1/2*a)*(1/6+1/3*a),1/2*a,1/3+2/3*a,1+1/2*a,1/16*a^2,(1/6+1/3*a)^2,(1/2+1/4*a)^2,1, 16*p^4*(1-p)^4,24*p^5*(1-p)^3,24*p^3*(1-p)^5,4*p^3*(1-p)^5,6*p^2*(1-p)^6,4*p*(1-p)^7,16*p^6*(1-p)^2,36*p^4*(1-p)^4,16*p^2*(1-p)^6,(1-p)^8, 64*p^4*(1-p)^4,120*p^5*(1-p)^3,72*p^3*(1-p)^5,12*p^3*(1-p)^5,12*p^2*(1-p)^6,4*p*(1-p)^7,96*p^6*(1-p)^2,144*p^4*(1-p)^4,32*p^2*(1-p)^6,0, 1/4*a*(1+2*a),1/12*a*(1+5*a),3/4*a*(1+a),1/2*a,a,3/2*a,1/8*a^2,1/3*a*(1/2+a),3/8*a*(2+a),0), nrow=10) }
MM_p <- function(p,a){ ML <- list(); ML[[1]] <- AAAAM_p(p=p,a=a); ML[[2]] <- AAAaM_p(p=p,a=a); ML[[3]] <- AAaaM_p(p=p,a=a); ML[[4]] <- AaaaM_p(p=p,a=a); ML[[5]] <- aaaaM_p(p=p,a=a); return(ML) }

auto4_est_c <- function(dose,cpar){
  nl1 <- rep(0,5); st <- table(dose); nl1[as.numeric(as.character(names(st)))+1] <- as.numeric(st)
  s1_pari <- cpar; iter1 <- 0
  while(1){
    pari1 <- s1_pari; r1 <- EME_p(p=pari1[1],a=pari1[2]); parii <- EMM_p1(ELL_p=r1,nl=nl1)
    tmpa <-  optimize(LL2,tol = 0.0000001,p=parii,nl=(nl1),lower=0,upper=0.25)
    s1_pari <- c(parii,tmpa$minimum); iter1 <- iter1 + 1
    if(max(abs(s1_pari-pari1))<1e-5||iter1>1000) break
  }
  s1_pari
}
