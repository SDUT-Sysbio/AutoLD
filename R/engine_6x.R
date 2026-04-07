# === Engine functions for Autohexaploid LD (Full & Partial) ===

auto6_est_LD <- function(dose){
  nts <- rep(0,49)
  nt <- table(dose)
  nts[as.numeric(names(nt))] <- as.numeric(nt)
  nn <- sum(nts)
  res4 <- EM_C6(nts=nts)

  pAAA <- sum(res4[1:4]); pAAa <- sum(res4[5:8]); pAaa <- sum(res4[9:12]); paaa <- sum(res4[13:16])
  pBBB <- sum(res4[c(1,5,9,13)]); pBBb <- sum(res4[c(2,6,10,14)]); pBbb <- sum(res4[c(3,7,11,15)]); pbbb <- sum(res4[c(4,8,12,16)])
  pA <- pAAA+pAAa*2/3+pAaa/3; pa <- paaa+pAaa*2/3+pAAa/3
  pB <- pBBB+pBBb*2/3+pBbb/3; pb <- pbbb+pBbb*2/3+pBBb/3

  rr1 <- solveD(pA=pA, pB=pB, pv=res4)
  estp <- c(pA,pB,rr1)
  gamA <- c(pAAA,pAAa,pAaa,paaa); names(gamA) <- c("AAA","AAa","Aaa","aaa")
  gamB <- c(pBBB,pBBb,pBbb,pbbb); names(gamB) <- c("AAA","AAa","Aaa","aaa")
  return(list(estp=estp,gamA=gamA,gamB=gamB))
}

auto6_est_LDp <- function(dose){
  nts <- rep(0,9)
  nt <- table(dose)
  nts[as.numeric(names(nt))] <- as.numeric(nt)
  nn <- sum(nts)
  res4 <- EM_C6p(nts=nts)

  pAAA <- sum(res4[1:4]); pAAa <- sum(res4[5:8]); pAaa <- sum(res4[9:12]); paaa <- sum(res4[13:16])
  pBBB <- sum(res4[c(1,5,9,13)]); pBBb <- sum(res4[c(2,6,10,14)]); pBbb <- sum(res4[c(3,7,11,15)]); pbbb <- sum(res4[c(4,8,12,16)])
  pA <- pAAA+pAAa*2/3+pAaa/3; pa <- paaa+pAaa*2/3+pAAa/3
  pB <- pBBB+pBBb*2/3+pBbb/3; pb <- pbbb+pBbb*2/3+pBBb/3

  rr1 <- solveD(pA=pA, pB=pB, pv=res4)
  estp <- c(pA,pB,rr1)
  gamA <- c(pAAA,pAAa,pAaa,paaa); names(gamA) <- c("AAA","AAa","Aaa","aaa")
  gamB <- c(pBBB,pBBb,pBbb,pbbb); names(gamB) <- c("AAA","AAa","Aaa","aaa")
  return(list(estp=estp,gamA=gamA,gamB=gamB))
}

EM_C6 <- function(nts){
  p4 <- rep(1/16,16); nn <- sum(nts); iter <- 0
  nts11 <- matrix(rep(nts,16),nrow=16,byrow=T)
  H <- exM6()
  while(1){
    p41 <- p4; phi <- ZIf(p=p4)
    emm <- cm6(HH=H, PHI=phi)
    p4 <- rowSums(emm*nts11)/(2*nn)
    iter <- iter + 1
    if(max(abs(p4-p41))<1e-8) break
  }
  return(p4)
}

EM_C6p <- function(nts){
  p4 <- rep(1/16,16); nn <- sum(nts); iter <- 0
  nts11 <- matrix(rep(nts,16),nrow=16,byrow=T)
  H1 <- exM6p()
  while(1){
    p41 <- p4; phi <- ZIfp(p=p4)
    emm <- cm6(HH=H1, PHI=phi)
    p4 <- rowSums(emm*nts11)/(2*nn)
    iter <- iter + 1
    if(max(abs(p4-p41))<1e-8) break
  }
  return(p4)
}

cm6 <- function(HH, PHI){
  ll <- length(HH); pc <- c()
  for(i in 1:ll){
    pc <- rbind(pc, rowSums(HH[[i]]*PHI))
  }
  pc
}

run_EM_engine_6x_hap <- function(nts, mode) {
  p4 <- rep(0.25, 4); nn <- sum(nts); iter <- 0
  while(iter < 2000) {
    p4_prev <- p4
    if(mode == "full") {
      phi <- EM_Em6_hap(p4); emm <- EM_Mm6_hap(phi)
    } else { emm <- EM_Emp6_hap(p4) }
    p1 <- sum(emm[1,] * nts) / (6 * nn); p2 <- sum(emm[2,] * nts) / (6 * nn)
    p3 <- sum(emm[3,] * nts) / (6 * nn); p4_new <- sum(emm[4,] * nts) / (6 * nn)
    p4 <- c(p1, p2, p3, p4_new)
    iter <- iter + 1
    if(max(abs(p4 - p4_prev)) < 1e-8) break
  }
  return(p4)
}



require(matlib)


auto6_est <- function(geno,type="full"){

  dose <- geno$class
  if(type=="full"){
    dp <- auto6_est_LD(dose=dose)
    a1 <- estimate_a_from_gamete6x_freq(p=dp$estp[1], dp$gamA)$a
    a2 <- estimate_a_from_gamete6x_freq(p=dp$estp[2], dp$gamB)$a
  }else{
    dp <- auto6_est_LDp(dose=dose)
    a1 <- estimate_a_from_gamete6x_freq(p=dp$estp[1], dp$gamA)$a
    a2 <- estimate_a_from_gamete6x_freq(p=dp$estp[2], dp$gamB)$a
  }

  init.par <- c(dp$estp[1:2],a1,a2,dp$estp[c(-c(1:4,6:7))])
  names(init.par) <- c("pA","pB","aA","aB","Deab","DAb","DaB","DAAb","DaBB","DAB","DAAB","DABB","DAABB")
  return(init.par)
}

solve_a_from_T_6x <- function(p, q_obs){
  # q_obs: named AAA,AAa,Aaa,aaa
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

auto6_est_LD <- function(dose){

  nts <- rep(0,49)
  nt <- table(dose)
  nts[as.numeric(names(nt))] <- as.numeric(nt)
  nn <- sum(nts)

  res4 <- EM_C6(nts=nts)

  #DA2,DB2,DA3,DB3
  pAAA <- sum(res4[1:4])
  pAAa <- sum(res4[5:8])
  pAaa <- sum(res4[9:12])
  paaa <- sum(res4[13:16])

  pBBB <- sum(res4[c(1,5,9,13)])
  pBBb <- sum(res4[c(2,6,10,14)])
  pBbb <- sum(res4[c(3,7,11,15)])
  pbbb <- sum(res4[c(4,8,12,16)])

  pA <- pAAA+pAAa*2/3+pAaa/3;pa <- paaa+pAaa*2/3+pAAa/3;
  pB <- pBBB+pBBb*2/3+pBbb/3;pb <- pbbb+pBbb*2/3+pBBb/3;

  rr1 <- solveD(pA=pA,pB=pB,pv=res4)

  estp <- c(pA,pB,rr1)
  gamA <- c(pAAA,pAAa,pAaa,paaa)
  names(gamA) <- c("AAA","AAa","Aaa","aaa")
  gamB <- c(pBBB,pBBb,pBbb,pbbb)
  names(gamB) <- c("AAA","AAa","Aaa","aaa")
  return(list(estp=estp,gamA=gamA,gamB=gamB))
}

solveD <- function(pA,pB,pv){


  HM <- H(pA=pA,pB=pB)
  HM1 <- HM[,-1]
  pv1 <- c(pv-HM[,1])

  #showEqn(HM1,pv1)
  #c(R(HM1),R(cbind(HM1,b)))

  #rr <- as.numeric(gaussianElimination(HM1, pv1, verbose = F, fractions = F)[1:13,14])
  rr <- MASS::ginv(HM1) %*% pv1
  # --- 修改结束 ---

  return(as.numeric(rr))
  #as.numeric(solve(HM1, pv1))
  #return(rr)
}

estimate_a_from_gamete6x_freq <- function(p, q_obs){
  # q_obs: named vector with names c("AAA","AAa","Aaa","aaa"), sums ~ 1

  q0 <- gamete_margin_DR_6x(p, a = 0)
  q1 <- gamete_margin_DR_6x(p, a = 1)

  # align order
  ord <- c("AAA","AAa","Aaa","aaa")
  q_obs <- q_obs[ord]; q0 <- q0[ord]; q1 <- q1[ord]

  d <- q1 - q0
  if(all(abs(d) < 1e-12)){
    return(list(a = NA_real_, reason = "Not identifiable: q(a=0) == q(a=1) for this p."))
  }

  a_hat <- sum((q_obs - q0) * d) / sum(d^2)
  a_hat <- max(0, min(1, a_hat))

  list(a = a_hat, q0 = q0, q1 = q1, d = d)
}


EM_C6 <- function(nts){

  p4 <- rep(1/16,16)
  nn <- sum(nts)
  iter <- 0
  nts11 <- matrix(rep(nts,16),nrow=16,byrow=T)
  H <- exM6()
  while(1){

    p41 <- p4
    phi <- ZIf(p=p4)
    emm <- cm6(HH=H,PHI=phi)
    p4 <- rowSums(emm*nts11)/(2*nn)

    iter <- iter + 1
    if(max(abs(p4-p41))<1e-8)
      break
  }
  return(p4)
}



exM6 <- function(){

  index <- c("PAAABBB","PAAABBb","PAAABbb","PAAAbbb","PAAaBBB","PAAaBBb","PAAaBbb","PAAabbb",
             "PAaaBBB","PAaaBBb","PAaaBbb","PAaabbb","PaaaBBB","PaaaBBb","PaaaBbb","Paaabbb")


  MM <- matrix(c("PAAABBB-PAAABBB",rep(0,7),"PAAABBB-PAAABBb",rep(0,7),"PAAABBB-PAAABbb","PAAABBb-PAAABBb",rep(0,6),
                 "PAAABBB-PAAAbbb","PAAABBb-PAAABbb",rep(0,6),"PAAABBb-PAAAbbb","PAAABbb-PAAABbb",rep(0,6),
                 "PAAABbb-PAAAbbb",rep(0,7),"PAAAbbb-PAAAbbb",rep(0,7),

                 "PAAABBB-PAAaBBB",rep(0,7),"PAAABBB-PAAaBBb","PAAaBBB-PAAABBb",rep(0,6),
                 "PAAABBB-PAAaBbb","PAAaBBB-PAAABbb","PAAABBb-PAAaBBb",rep(0,5),
                 "PAAABBB-PAAabbb","PAAaBBB-PAAAbbb","PAAABBb-PAAaBbb","PAAaBBb-PAAABbb",rep(0,4),
                 "PAAABBb-PAAabbb","PAAaBBb-PAAAbbb","PAAABbb-PAAaBbb",rep(0,5),
                 "PAAABbb-PAAabbb","PAAaBbb-PAAAbbb",rep(0,6),"PAAAbbb-PAAabbb",rep(0,7),

                 "PAAABBB-PAaaBBB","PAAaBBB-PAAaBBB",rep(0,6),
                 "PAAABBB-PAaaBBb","PAAaBBB-PAAaBBb","PAaaBBB-PAAABBb",rep(0,5),
                 "PAAaBBb-PAAaBBb","PAAABBB-PAaaBbb","PAAaBBB-PAAaBbb","PAaaBBB-PAAABbb","PAAABBb-PAaaBBb",rep(0,3),
                 "PAAABBB-PAaabbb","PAAaBBB-PAAabbb","PAaaBBB-PAAAbbb","PAAABBb-PAaaBbb","PAAaBBb-PAAaBbb","PAaaBBb-PAAABbb",rep(0,2),
                 "PAAaBbb-PAAaBbb","PAAABBb-PAaabbb","PAAaBBb-PAAabbb","PAaaBBb-PAAAbbb","PAAABbb-PAaaBbb",rep(0,3),
                 "PAAABbb-PAaabbb","PAaaBbb-PAAAbbb","PAAaBbb-PAAabbb",rep(0,5),"PAAabbb-PAAabbb","PAAAbbb-PAaabbb",rep(0,6),

                 "PAAABBB-PaaaBBB","PAAaBBB-PAaaBBB",rep(0,6),
                 "PAAABBB-PaaaBBb","PAAaBBB-PAaaBBb","PAaaBBB-PAAaBBb","PaaaBBB-PAAABBb",rep(0,4),
                 "PAAABBB-PaaaBbb","PAAaBBB-PAaaBbb","PAaaBBB-PAAaBbb","PaaaBBB-PAAABbb","PAAABBb-PaaaBBb","PAAaBBb-PAaaBBb",rep(0,2),
                 "PAAABBB-Paaabbb","PAAaBBB-PAaabbb","PAaaBBB-PAAabbb","PaaaBBB-PAAAbbb","PAAABBb-PaaaBbb","PAAaBBb-PAaaBbb","PAaaBBb-PAAaBbb","PaaaBBb-PAAABbb",
                 "PAAABBb-Paaabbb","PAAaBBb-PAaabbb","PAaaBBb-PAAabbb","PaaaBBb-PAAAbbb","PAAABbb-PaaaBbb","PAAaBbb-PAaaBbb",rep(0,2),
                 "PAAABbb-Paaabbb","PAAaBbb-PAaabbb","PAaaBbb-PAAabbb","PaaaBbb-PAAAbbb",rep(0,4),"PAAAbbb-Paaabbb","PAAabbb-PAaabbb",rep(0,6),

                 "PAaaBBB-PAaaBBB","PAAaBBB-PaaaBBB",rep(0,6),"PAAaBBB-PaaaBBb","PAaaBBB-PAaaBBb","PaaaBBB-PAAaBBb",rep(0,5),
                 "PAaaBBb-PAaaBBb","PAAaBBB-PaaaBbb","PAaaBBB-PAaaBbb","PaaaBBB-PAAaBbb","PAAaBBb-PaaaBBb",rep(0,3),
                 "PAAaBBB-Paaabbb","PAaaBBB-PAaabbb","PaaaBBB-PAAabbb","PAAaBBb-PaaaBbb","PAaaBBb-PAaaBbb","PaaaBBb-PAAaBbb",rep(0,2),
                 "PAaaBbb-PAaaBbb","PAAaBBb-Paaabbb","PAaaBBb-PAaabbb","PaaaBBb-PAAabbb","PAAaBbb-PaaaBbb",rep(0,3),
                 "PAAaBbb-Paaabbb","PAaaBbb-PAaabbb","PaaaBbb-PAAabbb",rep(0,5),
                 "PAaabbb-PAaabbb","PAAabbb-Paaabbb",rep(0,6),

                 "PAaaBBB-PaaaBBB",rep(0,7),"PAaaBBB-PaaaBBb","PaaaBBB-PAaaBBb",rep(0,6),
                 "PAaaBBB-PaaaBbb","PaaaBBB-PAaaBbb","PAaaBBb-PaaaBBb",rep(0,5),
                 "PAaaBBB-Paaabbb","PaaaBBB-PAaabbb","PAaaBBb-PaaaBbb","PaaaBBb-PAaaBbb",rep(0,4),
                 "PAaaBBb-Paaabbb","PaaaBBb-PAaabbb","PAaaBbb-PaaaBbb",rep(0,5),"PAaaBbb-Paaabbb","PaaaBbb-PAaabbb",rep(0,6),
                 "PAaabbb-Paaabbb",rep(0,7),

                 "PaaaBBB-PaaaBBB",rep(0,7),"PaaaBBB-PaaaBBb",rep(0,7),"PaaaBBb-PaaaBBb","PaaaBBB-PaaaBbb",rep(0,6),
                 "PaaaBBB-Paaabbb","PaaaBBb-PaaaBbb",rep(0,6),"PaaaBbb-PaaaBbb","PaaaBBb-Paaabbb",rep(0,6),"PaaaBbb-Paaabbb",rep(0,7),
                 "Paaabbb-Paaabbb",rep(0,7)),ncol=8,byrow=T)

  n1 <- dim(MM)[1];n2 <- dim(MM)[2];


  MIL <- list()
  ii <- 1
  for(i in index){
    MI <- matrix(0,nrow=n1,ncol=n2)
    for(j in 1:n1){
      for(k in 1:n2){
        MI[j,k] <- length(which(strsplit(MM[j,k],split="-")[[1]]==i))
      }
    }
    MIL[[ii]] <- MI
    #print(sum(MI))
    ii <- ii + 1
  }
  return(MIL)
}



ZIf <- function(p){
  PAAABBB <- p[1]
  PAAABBb <- p[2]
  PAAABbb <- p[3]
  PAAAbbb <- p[4]
  PAAaBBB <- p[5]
  PAAaBBb <- p[6]
  PAAaBbb <- p[7]
  PAAabbb <- p[8]
  PAaaBBB <- p[9]
  PAaaBBb <- p[10]
  PAaaBbb <- p[11]
  PAaabbb <- p[12]
  PaaaBBB <- p[13]
  PaaaBBb <- p[14]
  PaaaBbb <- p[15]
  Paaabbb <- p[16]

  zp <- matrix(c(PAAABBB^2,rep(0,7),2*PAAABBB*PAAABBb,rep(0,7),2*PAAABBB*PAAABbb,(PAAABBb)^2,rep(0,6),
                 2*PAAABBB*PAAAbbb,2*PAAABBb*PAAABbb,rep(0,6),2*PAAABBb*PAAAbbb,(PAAABbb)^2,rep(0,6),
                 2*PAAABbb*PAAAbbb,rep(0,7),(PAAAbbb)^2,rep(0,7),
                 2*PAAABBB*PAAaBBB,rep(0,7),2*PAAABBB*PAAaBBb,2*PAAaBBB*PAAABBb,rep(0,6),
                 2*PAAABBB*PAAaBbb,2*PAAaBBB*PAAABbb,2*PAAABBb*PAAaBBb,rep(0,5),
                 2*PAAABBB*PAAabbb,2*PAAaBBB*PAAAbbb,2*PAAABBb*PAAaBbb,2*PAAaBBb*PAAABbb,rep(0,4),
                 2*PAAABBb*PAAabbb,2*PAAaBBb*PAAAbbb,2*PAAABbb*PAAaBbb,rep(0,5),
                 2*PAAABbb*PAAabbb,2*PAAaBbb*PAAAbbb,rep(0,6),2*PAAAbbb*PAAabbb,rep(0,7),
                 2*PAAABBB*PAaaBBB,(PAAaBBB)^2,rep(0,6),
                 2*PAAABBB*PAaaBBb,2*PAAaBBB*PAAaBBb,2*PAaaBBB*PAAABBb,rep(0,5),
                 (PAAaBBb)^2,2*PAAABBB*PAaaBbb,2*PAAaBBB*PAAaBbb,2*PAaaBBB*PAAABbb,2*PAAABBb*PAaaBBb,rep(0,3),
                 2*PAAABBB*PAaabbb,2*PAAaBBB*PAAabbb,2*PAaaBBB*PAAAbbb,2*PAAABBb*PAaaBbb,2*PAAaBBb*PAAaBbb,2*PAaaBBb*PAAABbb,rep(0,2),
                 (PAAaBbb)^2,2*PAAABBb*PAaabbb,2*PAAaBBb*PAAabbb,2*PAaaBBb*PAAAbbb,2*PAAABbb*PAaaBbb,rep(0,3),
                 2*PAAABbb*PAaabbb,2*PAaaBbb*PAAAbbb,2*PAAaBbb*PAAabbb,rep(0,5),(PAAabbb)^2,2*PAAAbbb*PAaabbb,rep(0,6),
                 2*PAAABBB*PaaaBBB,2*PAAaBBB*PAaaBBB,rep(0,6),
                 2*PAAABBB*PaaaBBb,2*PAAaBBB*PAaaBBb,2*PAaaBBB*PAAaBBb,2*PaaaBBB*PAAABBb,rep(0,4),
                 2*PAAABBB*PaaaBbb,2*PAAaBBB*PAaaBbb,2*PAaaBBB*PAAaBbb,2*PaaaBBB*PAAABbb,2*PAAABBb*PaaaBBb,2*PAAaBBb*PAaaBBb,rep(0,2),
                 2*PAAABBB*Paaabbb,2*PAAaBBB*PAaabbb,2*PAaaBBB*PAAabbb,2*PaaaBBB*PAAAbbb,2*PAAABBb*PaaaBbb,2*PAAaBBb*PAaaBbb,2*PAaaBBb*PAAaBbb,2*PaaaBBb*PAAABbb,
                 2*PAAABBb*Paaabbb,2*PAAaBBb*PAaabbb,2*PAaaBBb*PAAabbb,2*PaaaBBb*PAAAbbb,2*PAAABbb*PaaaBbb,2*PAAaBbb*PAaaBbb,rep(0,2),
                 2*PAAABbb*Paaabbb,2*PAAaBbb*PAaabbb,2*PAaaBbb*PAAabbb,2*PaaaBbb*PAAAbbb,rep(0,4),2*PAAAbbb*Paaabbb,2*PAAabbb*PAaabbb,rep(0,6),
                 (PAaaBBB)^2,2*PAAaBBB*PaaaBBB,rep(0,6),2*PAAaBBB*PaaaBBb,2*PAaaBBB*PAaaBBb,2*PaaaBBB*PAAaBBb,rep(0,5),
                 (PAaaBBb)^2,2*PAAaBBB*PaaaBbb,2*PAaaBBB*PAaaBbb,2*PaaaBBB*PAAaBbb,2*PAAaBBb*PaaaBBb,rep(0,3),
                 2*PAAaBBB*Paaabbb,2*PAaaBBB*PAaabbb,2*PaaaBBB*PAAabbb,2*PAAaBBb*PaaaBbb,2*PAaaBBb*PAaaBbb,2*PaaaBBb*PAAaBbb,rep(0,2),
                 (PAaaBbb)^2,2*PAAaBBb*Paaabbb,2*PAaaBBb*PAaabbb,2*PaaaBBb*PAAabbb,2*PAAaBbb*PaaaBbb,rep(0,3),
                 2*PAAaBbb*Paaabbb,2*PAaaBbb*PAaabbb,2*PaaaBbb*PAAabbb,rep(0,5),
                 (PAaabbb)^2,2*PAAabbb*Paaabbb,rep(0,6),
                 2*PAaaBBB*PaaaBBB,rep(0,7),2*PAaaBBB*PaaaBBb,2*PaaaBBB*PAaaBBb,rep(0,6),
                 2*PAaaBBB*PaaaBbb,2*PaaaBBB*PAaaBbb,2*PAaaBBb*PaaaBBb,rep(0,5),
                 2*PAaaBBB*Paaabbb,2*PaaaBBB*PAaabbb,2*PAaaBBb*PaaaBbb,2*PaaaBBb*PAaaBbb,rep(0,4),
                 2*PAaaBBb*Paaabbb,2*PaaaBBb*PAaabbb,2*PAaaBbb*PaaaBbb,rep(0,5),2*PAaaBbb*Paaabbb,2*PaaaBbb*PAaabbb,rep(0,6),
                 2*PAaabbb*Paaabbb,rep(0,7),(PaaaBBB)^2,rep(0,7),2*PaaaBBB*PaaaBBb,rep(0,7),(PaaaBBb)^2,2*PaaaBBB*PaaaBbb,rep(0,6),
                 2*PaaaBBB*Paaabbb,2*PaaaBBb*PaaaBbb,rep(0,6),(PaaaBbb)^2,2*PaaaBBb*Paaabbb,rep(0,6),2*PaaaBbb*Paaabbb,rep(0,7),
                 (Paaabbb)^2,rep(0,7)),ncol=8,byrow=T)

  zs <-  cbind(zp,rowSums(zp))

  conp <- t(apply(zs,1,function(x){
    if(x[9]==0){
      x[9] <- 1e-200
    }
    x[1:8]/x[9]

  }))

  conp
}

cm6 <- function(HH,PHI){

  ll <- length(HH)
  pc <- c()
  for(i in 1:ll){
    pc <- rbind(pc,rowSums(HH[[i]]*PHI))
  }
  pc
}


sim_6x_unified_hap <- function(pall, n=1000, mode="partial") {

  gametes <- std_6x_gametes_hap(pall)

  prob_49 <- std_6x_zygotes_49_hap(gametes)

  if(mode == "full") {
    return(sample(1:49, n, replace=TRUE, prob=prob_49))

  } else if(mode == "partial") {
    prob_9 <- std_6x_collapse_to_9_hap(prob_49)
    return(sample(1:9, n, replace=TRUE, prob=prob_9))

  } else {
    stop("Invalid mode. Use 'full' or 'partial'.")
  }
}

est_6x_unified_hap <- function(geno_data, mode="partial") {

  nn <- length(geno_data)
  p4 <- rep(0.25, 4)

  if(mode == "full") {
    nts <- rep(0, 49)
    nt_tab <- table(geno_data)
    nts[as.numeric(names(nt_tab))] <- as.numeric(nt_tab)

    res <- run_EM_engine_6x_hap(nts, mode="full")

  } else if(mode == "partial") {
    nts <- rep(0, 9)
    nt_tab <- table(geno_data)
    nts[as.numeric(names(nt_tab))] <- as.numeric(nt_tab)

    res <- run_EM_engine_6x_hap(nts, mode="partial")

  } else {
    stop("Invalid mode. Use 'full' or 'partial'.")
  }

  pA <- sum(res[1:2])
  pB <- sum(res[c(1,3)])
  D  <- res[1] - pA*pB

  return(c(pA=pA, pB=pB, D=D))
}


run_EM_engine_6x_hap <- function(nts, mode) {
  p4 <- rep(0.25, 4)
  nn <- sum(nts)
  iter <- 0

  while(iter < 2000) {
    p4_prev <- p4

    if(mode == "full") {
      phi <- EM_Em6_hap(p4)
      emm <- EM_Mm6_hap(phi)
    } else {
      emm <- EM_Emp6_hap(p4)
    }

    # 更新频率 (M-step)
    p1 <- sum(emm[1,] * nts) / (6 * nn)
    p2 <- sum(emm[2,] * nts) / (6 * nn)
    p3 <- sum(emm[3,] * nts) / (6 * nn)
    p4_new <- sum(emm[4,] * nts) / (6 * nn)

    p4 <- c(p1, p2, p3, p4_new)

    iter <- iter + 1
    if(max(abs(p4 - p4_prev)) < 1e-8) break
  }
  return(p4)
}

std_6x_gametes_hap <- function(pall){
  pA <- pall[1]; pa <- 1-pA
  pB <- pall[2]; pb <- 1-pB
  D <- pall[3]

  hAB <- pA*pB+D; hAb <- pA*pb-D
  haB <- pa*pB-D; hab <- pa*pb+D

  g <- numeric(16)
  g[1] <- hAB^3
  g[2] <- 3*hAB^2*hAb
  g[3] <- 3*hAB*hAb^2
  g[4] <- hAb^3
  g[5] <- 3*hAB^2*haB
  g[6] <- 3*hAB^2*hab + 6*hAB*hAb*haB
  g[7] <- 3*hAb^2*haB + 6*hAB*hAb*hab
  g[8] <- 3*hAb^2*hab
  g[9] <- 3*hAB*haB^2
  g[10] <- 3*hAb*haB^2 + 6*hAB*haB*hab
  g[11] <- 3*hAB*hab^2 + 6*hAb*haB*hab
  g[12] <- 3*hAb*hab^2
  g[13] <- haB^3
  g[14] <- 3*haB^2*hab
  g[15] <- 3*haB*hab^2
  g[16] <- hab^3
  return(g)
}

std_6x_zygotes_49_hap <- function(gp){

  PAAABBB=gp[1]; PAAABBb=gp[2]; PAAABbb=gp[3]; PAAAbbb=gp[4]
  PAAaBBB=gp[5]; PAAaBBb=gp[6]; PAAaBbb=gp[7]; PAAabbb=gp[8]
  PAaaBBB=gp[9]; PAaaBBb=gp[10]; PAaaBbb=gp[11]; PAaabbb=gp[12]
  PaaaBBB=gp[13]; PaaaBBb=gp[14]; PaaaBbb=gp[15]; Paaabbb=gp[16]

  zy <- numeric(49)
  zy[1] <- (PAAABBB)^2
  zy[2] <- 2*PAAABBB*PAAABBb
  zy[3] <- 2*PAAABBB*PAAABbb+(PAAABBb)^2
  zy[4] <- 2*PAAABBB*PAAAbbb+2*PAAABBb*PAAABbb
  zy[5] <- 2*PAAABBb*PAAAbbb+(PAAABbb)^2
  zy[6] <- 2*PAAABbb*PAAAbbb
  zy[7] <- (PAAAbbb)^2
  zy[8] <- 2*PAAABBB*PAAaBBB
  zy[9] <- 2*PAAABBB*PAAaBBb+2*PAAaBBB*PAAABBb
  zy[10] <- 2*PAAABBB*PAAaBbb+2*PAAaBBB*PAAABbb+2*PAAABBb*PAAaBBb
  zy[11] <- 2*PAAABBB*PAAabbb+2*PAAaBBB*PAAAbbb+2*PAAABBb*PAAaBbb+2*PAAaBBb*PAAABbb
  zy[12] <- 2*PAAABBb*PAAabbb+2*PAAaBBb*PAAAbbb+2*PAAABbb*PAAaBbb
  zy[13] <- 2*PAAABbb*PAAabbb+2*PAAaBbb*PAAAbbb
  zy[14] <- 2*PAAAbbb*PAAabbb
  zy[15] <- 2*PAAABBB*PAaaBBB+(PAAaBBB)^2
  zy[16] <- 2*PAAABBB*PAaaBBb+2*PAAaBBB*PAAaBBb+2*PAaaBBB*PAAABBb
  zy[17] <- (PAAaBBb)^2+2*PAAABBB*PAaaBbb+2*PAAaBBB*PAAaBbb+2*PAaaBBB*PAAABbb+2*PAAABBb*PAaaBBb
  zy[18] <- 2*PAAABBB*PAaabbb+2*PAAaBBB*PAAabbb+2*PAaaBBB*PAAAbbb+2*PAAABBb*PAaaBbb+2*PAAaBBb*PAAaBbb+2*PAaaBBb*PAAABbb
  zy[19] <- (PAAaBbb)^2+2*PAAABBb*PAaabbb+2*PAAaBBb*PAAabbb+2*PAaaBBb*PAAAbbb+2*PAAABbb*PAaaBbb
  zy[20] <- 2*PAAABbb*PAaabbb+2*PAaaBbb*PAAAbbb+2*PAAaBbb*PAAabbb
  zy[21] <- (PAAabbb)^2+2*PAAAbbb*PAaabbb
  zy[22] <- 2*PAAABBB*PaaaBBB+2*PAAaBBB*PAaaBBB
  zy[23] <- 2*PAAABBB*PaaaBBb+2*PAAaBBB*PAaaBBb+2*PAaaBBB*PAAaBBb+2*PaaaBBB*PAAABBb
  zy[24] <- 2*PAAABBB*PaaaBbb+2*PAAaBBB*PAaaBbb+2*PAaaBBB*PAAaBbb+2*PaaaBBB*PAAABbb+2*PAAABBb*PaaaBBb+2*PAAaBBb*PAaaBBb
  zy[25] <- 2*PAAABBB*Paaabbb+2*PAAaBBB*PAaabbb+2*PAaaBBB*PAAabbb+2*PaaaBBB*PAAAbbb+2*PAAABBb*PaaaBbb+2*PAAaBBb*PAaaBbb+2*PAaaBBb*PAAaBbb+2*PaaaBBb*PAAABbb
  zy[26] <- 2*PAAABBb*Paaabbb+2*PAAaBBb*PAaabbb+2*PAaaBBb*PAAabbb+2*PaaaBBb*PAAAbbb+2*PAAABbb*PaaaBbb+2*PAAaBbb*PAaaBbb
  zy[27] <- 2*PAAABbb*Paaabbb+2*PAAaBbb*PAaabbb+2*PAaaBbb*PAAabbb+2*PaaaBbb*PAAAbbb
  zy[28] <- 2*PAAAbbb*Paaabbb+2*PAAabbb*PAaabbb
  zy[29] <- (PAaaBBB)^2+2*PAAaBBB*PaaaBBB
  zy[30] <- 2*PAAaBBB*PaaaBBb+2*PAaaBBB*PAaaBBb+2*PaaaBBB*PAAaBBb
  zy[31] <- (PAaaBBb)^2+2*PAAaBBB*PaaaBbb+2*PAaaBBB*PAaaBbb+2*PaaaBBB*PAAaBbb+2*PAAaBBb*PaaaBBb
  zy[32] <- 2*PAAaBBB*Paaabbb+2*PAaaBBB*PAaabbb+2*PaaaBBB*PAAabbb+2*PAAaBBb*PaaaBbb+2*PAaaBBb*PAaaBbb+2*PaaaBBb*PAAaBbb
  zy[33] <- (PAaaBbb)^2+2*PAAaBBb*Paaabbb+2*PAaaBBb*PAaabbb+2*PaaaBBb*PAAabbb+2*PAAaBbb*PaaaBbb
  zy[34] <- 2*PAAaBbb*Paaabbb+2*PAaaBbb*PAaabbb+2*PaaaBbb*PAAabbb
  zy[35] <- (PAaabbb)^2+2*PAAabbb*Paaabbb
  zy[36] <- 2*PAaaBBB*PaaaBBB
  zy[37] <- 2*PAaaBBB*PaaaBBb+2*PaaaBBB*PAaaBBb
  zy[38] <- 2*PAaaBBB*PaaaBbb+2*PaaaBBB*PAaaBbb+2*PAaaBBb*PaaaBBb
  zy[39] <- 2*PAaaBBB*Paaabbb+2*PaaaBBB*PAaabbb+2*PAaaBBb*PaaaBbb+2*PaaaBBb*PAaaBbb
  zy[40] <- 2*PAaaBBb*Paaabbb+2*PaaaBBb*PAaabbb+2*PAaaBbb*PaaaBbb
  zy[41] <- 2*PAaaBbb*Paaabbb+2*PaaaBbb*PAaabbb
  zy[42] <- 2*PAaabbb*Paaabbb
  zy[43] <- (PaaaBBB)^2
  zy[44] <- 2*PaaaBBB*PaaaBBb
  zy[45] <- (PaaaBBb)^2+2*PaaaBBB*PaaaBbb
  zy[46] <- 2*PaaaBBB*Paaabbb+2*PaaaBBb*PaaaBbb
  zy[47] <- (PaaaBbb)^2+2*PaaaBBb*Paaabbb
  zy[48] <- 2*PaaaBbb*Paaabbb
  zy[49] <- (Paaabbb)^2

  return(zy)
}

std_6x_collapse_to_9_hap <- function(prob_49){


  prob_mat <- matrix(prob_49, nrow=7, ncol=7)

  r1 <- prob_mat[1,]
  r2 <- colSums(prob_mat[2:6,])
  r3 <- prob_mat[7,]
  mat_3x7 <- rbind(r1, r2, r3)

  c1 <- mat_3x7[,1]
  c2 <- rowSums(mat_3x7[,2:6])
  c3 <- mat_3x7[,7]
  mat_3x3 <- cbind(c1, c2, c3)

  # 展平为 9 个元素的向量 (默认按列)
  return(as.numeric(mat_3x3))
}


EM_Em6_hap <- function(p){
  PAB <- p[1]; PAb <- p[2]; PaB <- p[3]; Pab <- p[4]

  phi9_AB <- (6*PAB^5*Pab)/(6*PAB^5*Pab+30*PAB^4*PAb*PaB)
  phi10_AB <- (30*PAB^4*PAb*Pab)/(30*PAB^4*PAb*Pab+60*PAB^3*PAb^2*PaB)
  phi11_AB <- (60*PAB^3*PAb^2*Pab)/(60*PAB^3*PAb^2*Pab+60*PAB^2*PAb^3*PaB)
  phi12_AB <- (60*PAB^2*PAb^3*Pab)/(60*PAB^2*PAb^3*Pab+30*PAB*PAb^4*PaB)
  phi13_AB <- (30*PAB*PAb^4*Pab)/(30*PAB*PAb^4*Pab+6*PAb^5*PaB)
  phi16_AB <- (30*PAB^4*PaB*Pab)/(30*PAB^4*PaB*Pab+60*PAB^3*PAb*PaB^2)
  phi17_AB1 <- (15*PAB^4*Pab^2)/(15*PAB^4*Pab^2+120*PAB^3*PAb*PaB*Pab+90*PAB^2*PAb^2*PaB^2)
  phi17_AB2<- (120*PAB^3*PAb*PaB*Pab)/(15*PAB^4*Pab^2+120*PAB^3*PAb*PaB*Pab+90*PAB^2*PAb^2*PaB^2)
  phi18_AB1 <- (60*PAB^3*PAb*Pab^2)/(60*PAB^3*PAb*Pab^2+180*PAB^2*PAb^2*PaB*Pab+60*PAB*PAb^3*PaB^2)
  phi18_AB2<- (180*PAB^2*PAb^2*PaB*Pab)/(60*PAB^3*PAb*Pab^2+180*PAB^2*PAb^2*PaB*Pab+60*PAB*PAb^3*PaB^2)
  phi19_AB1 <- (90*PAB^2*PAb^2*Pab^2)/(90*PAB^2*PAb^2*Pab^2+120*PAB*PAb^3*PaB*Pab+15*PAb^4*PaB^2)
  phi19_AB2<- (120*PAB*PAb^3*PaB*Pab)/(90*PAB^2*PAb^2*Pab^2+120*PAB*PAb^3*PaB*Pab+15*PAb^4*PaB^2)
  phi20_AB <- (60*PAB*PAb^3*Pab^2)/(60*PAB*PAb^3*Pab^2+30*PAb^4*PaB*Pab)
  phi23_AB <- (60*PAB^3*PaB^2*Pab)/(60*PAB^3*PaB^2*Pab+60*PAB^2*PAb*PaB^3)
  phi24_AB1 <- (60*PAB^3*PaB*Pab^2)/(60*PAB^3*PaB*Pab^2+180*PAB^2*PAb*PaB^2*Pab+60*PAB*PAb^2*PaB^3)
  phi24_AB2 <- (180*PAB^2*PAb*PaB^2*Pab)/(60*PAB^3*PaB*Pab^2+180*PAB^2*PAb*PaB^2*Pab+60*PAB*PAb^2*PaB^3)
  phi25_AB1 <- (20*PAB^3*Pab^3)/(20*PAB^3*Pab^3+180*PAB^2*PAb*PaB*Pab^2+180*PAB*PAb^2*PaB^2*Pab+20*PAb^3*PaB^3)
  phi25_AB2 <- (180*PAB^2*PAb*PaB*Pab^2)/(20*PAB^3*Pab^3+180*PAB^2*PAb*PaB*Pab^2+180*PAB*PAb^2*PaB^2*Pab+20*PAb^3*PaB^3)
  phi25_AB3 <- (180*PAB*PAb^2*PaB^2*Pab)/(20*PAB^3*Pab^3+180*PAB^2*PAb*PaB*Pab^2+180*PAB*PAb^2*PaB^2*Pab+20*PAb^3*PaB^3)
  phi26_AB1 <- (60*PAB^2*PAb*Pab^3)/(60*PAB^2*PAb*Pab^3+180*PAB*PAb^2*PaB*Pab^2+60*PAb^3*PaB^2*Pab)
  phi26_AB2 <- (180*PAB*PAb^2*PaB*Pab^2)/(60*PAB^2*PAb*Pab^3+180*PAB*PAb^2*PaB*Pab^2+60*PAb^3*PaB^2*Pab)
  phi27_AB <- (60*PAB*PAb^2*Pab^3)/(60*PAB*PAb^2*Pab^3+60*PAb^3*PaB*Pab^2)
  phi30_AB <- (60*PAB^2*PaB^3*Pab)/(60*PAB^2*PaB^3*Pab+30*PAB*PAb*PaB^4)
  phi31_AB1 <- (90*PAB^2*PaB^2*Pab^2)/(90*PAB^2*PaB^2*Pab^2+120*PAB*PAb*PaB^3*Pab+15*PAb^2*PaB^4)
  phi31_AB2 <- (120*PAB*PAb*PaB^3*Pab)/(90*PAB^2*PaB^2*Pab^2+120*PAB*PAb*PaB^3*Pab+15*PAb^2*PaB^4)
  phi32_AB1 <- (60*PAB^2*PaB*Pab^3)/(60*PAB^2*PaB*Pab^3+180*PAB*PAb*PaB^2*Pab^2+60*PAb^2*PaB^3*Pab)
  phi32_AB2 <- (180*PAB*PAb*PaB^2*Pab^2)/(60*PAB^2*PaB*Pab^3+180*PAB*PAb*PaB^2*Pab^2+60*PAb^2*PaB^3*Pab)
  phi33_AB1 <- (15*PAB^2*Pab^4)/(15*PAB^2*Pab^4+120*PAB*PAb*PaB*Pab^3+90*PAb^2*PaB^2*Pab^2)
  phi33_AB2 <- (120*PAB*PAb*PaB*Pab^3)/(15*PAB^2*Pab^4+120*PAB*PAb*PaB*Pab^3+90*PAb^2*PaB^2*Pab^2)
  phi34_AB <- (30*PAB*PAb*Pab^4)/(30*PAB*PAb*Pab^4+60*PAb^2*PaB*Pab^3)
  phi37_AB <- (30*PAB*PaB^4*Pab)/(30*PAB*PaB^4*Pab+6*PAb*PaB^5)
  phi38_AB <- (60*PAB*PaB^3*Pab^2)/(60*PAB*PaB^3*Pab^2+30*PAb*PaB^4*Pab)
  phi39_AB <- (60*PAB*PaB^2*Pab^3)/(60*PAB*PaB^2*Pab^3+60*PAb*PaB^3*Pab^2)
  phi40_AB <- (30*PAB*PaB*Pab^4)/(30*PAB*PaB*Pab^4+60*PAb*PaB^2*Pab^3)
  phi41_AB <- (6*PAB*Pab^5)/(6*PAB*Pab^5+30*PAb*PaB*Pab^4)

  vars <- ls(pattern="phi")
  for(v in vars) {
    val <- get(v)
    if(is.nan(val)) assign(v, 0)
  }

  return(c(phi9_AB,phi10_AB,phi11_AB,phi12_AB,phi13_AB,phi16_AB,phi17_AB1,phi17_AB2,phi18_AB1,phi18_AB2,phi19_AB1,phi19_AB2,phi20_AB,
           phi23_AB,phi24_AB1,phi24_AB2,phi25_AB1,phi25_AB2,phi25_AB3,phi26_AB1,phi26_AB2,phi27_AB,phi30_AB,phi31_AB1,phi31_AB2,
           phi32_AB1,phi32_AB2,phi33_AB1,phi33_AB2,phi34_AB,phi37_AB,phi38_AB,phi39_AB,phi40_AB,phi41_AB))
}


EM_Mm6_hap <- function(phi){

  phi9_AB <- phi[1];phi10_AB <- phi[2];phi11_AB <- phi[3];phi12_AB <- phi[4];phi13_AB <- phi[5]
  phi16_AB <- phi[6];phi17_AB1 <- phi[7];phi17_AB2<- phi[8];phi18_AB1 <- phi[9];phi18_AB2<- phi[10]; phi19_AB1 <- phi[11];phi19_AB2<- phi[12];phi20_AB <- phi[13]
  phi23_AB <- phi[14];phi24_AB1 <- phi[15];phi24_AB2 <- phi[16];phi25_AB1 <- phi[17];phi25_AB2 <- phi[18];phi25_AB3 <- phi[19];
  phi26_AB1 <- phi[20];phi26_AB2 <- phi[21];phi27_AB <- phi[22]
  phi30_AB <- phi[23];phi31_AB1 <- phi[24];phi31_AB2 <- phi[25];phi32_AB1 <- phi[26];phi32_AB2 <- phi[27];phi33_AB1 <- phi[28];phi33_AB2 <- phi[29];phi34_AB <- phi[30]
  phi37_AB <- phi[31];phi38_AB <- phi[32];phi39_AB <- phi[33];phi40_AB <- phi[34];phi41_AB <- phi[35]

  epAB <- c(6,5,4,3,2,1,0,
            5,4+phi9_AB,3+phi10_AB,2+phi11_AB,1+phi12_AB,phi13_AB,0,
            4,3+phi16_AB,2+2*phi17_AB1+phi17_AB2,1+2*phi18_AB1+phi18_AB2,2*phi19_AB1+phi19_AB2,phi20_AB,0,
            3,2+phi23_AB,1+2*phi24_AB1+phi24_AB2,3*phi25_AB1+2*phi25_AB2+phi25_AB3,2*phi26_AB1+phi26_AB2,phi27_AB,0,
            2,1+phi30_AB,2*phi31_AB1+phi31_AB2,2*phi32_AB1+phi32_AB2,2*phi33_AB1+phi33_AB2,phi34_AB,0,
            1,phi37_AB,phi38_AB,phi39_AB,phi40_AB,phi41_AB,0,
            0,0,0,0,0,0,0)

  epAb <- c(0,1,2,3,4,5,6,
            0,1-phi9_AB,2-phi10_AB,3-phi11_AB,4-phi12_AB,5-phi13_AB,5,
            0,1-phi16_AB,2-2*phi17_AB1-phi17_AB2,3-2*phi18_AB1-phi18_AB2,4-2*phi19_AB1-phi19_AB2,4-phi20_AB,4,
            0,1-phi23_AB,2-2*phi24_AB1-phi24_AB2,3-3*phi25_AB1-2*phi25_AB2-phi25_AB3,3-2*phi26_AB1-phi26_AB2,3-phi27_AB,3,
            0,1-phi30_AB,2-2*phi31_AB1-phi31_AB2,2-2*phi32_AB1-phi32_AB2,2-2*phi33_AB1-phi33_AB2,2-phi34_AB,2,
            0,1-phi37_AB,1-phi38_AB,1-phi39_AB,1-phi40_AB,1-phi41_AB,1,
            0,0,0,0,0,0,0)

  epaB<- c(0,0,0,0,0,0,0,
           1,1-phi9_AB,1-phi10_AB,1-phi11_AB,1-phi12_AB,1-phi13_AB,0,
           2,2-phi16_AB,2-2*phi17_AB1-phi17_AB2,2-2*phi18_AB1-phi18_AB2,2-2*phi19_AB1-phi19_AB2,1-phi20_AB,0,
           3,3-phi23_AB,3-2*phi24_AB1-phi24_AB2,3-3*phi25_AB1-2*phi25_AB2-phi25_AB3,2-2*phi26_AB1-phi26_AB2,1-phi27_AB,0,
           4,4-phi30_AB,4-2*phi31_AB1-phi31_AB2,3-2*phi32_AB1-phi32_AB2,2-2*phi33_AB1-phi33_AB2,1-phi34_AB,0,
           5,5-phi37_AB,4-phi38_AB,3-phi39_AB,2-phi40_AB,1-phi41_AB,0,
           6,5,4,3,2,1,0)

  epab<- c(0,0,0,0,0,0,0,
           0,phi9_AB,phi10_AB,phi11_AB,phi12_AB,phi13_AB,1,
           0,phi16_AB,2*phi17_AB1+phi17_AB2,2*phi18_AB1+phi18_AB2,2*phi19_AB1+phi19_AB2,1+phi20_AB,2,
           0,phi23_AB,2*phi24_AB1+phi24_AB2,3*phi25_AB1+2*phi25_AB2+phi25_AB3,1+2*phi26_AB1+phi26_AB2,2+phi27_AB,3,
           0,phi30_AB,2*phi31_AB1+phi31_AB2,1+2*phi32_AB1+phi32_AB2,2+2*phi33_AB1+phi33_AB2,3+phi34_AB,4,
           0,phi37_AB,1+phi38_AB,2+phi39_AB,3+phi40_AB,4+phi41_AB,5,
           0,1,2,3,4,5,6)

  return(rbind(epAB,epAb,epaB,epab))
}

EM_Emp6_hap <- function(p){
  PAB <- p[1]; PAb <- p[2]; PaB <- p[3]; Pab <- p[4]

  ai2 <- c(6*PAB^5*PAb,15*PAB^4*PAb^2,20*PAB^3*PAb^3,15*PAB^2*PAb^4,6*PAB^1*PAb^5)
  pr2 <- ai2/sum(ai2); if(sum(ai2)==0) pr2<-rep(0,5)
  phi2_AB <- sum(c(5,4,3,2,1)*pr2)
  phi2_Ab <- sum(c(1,2,3,4,5)*pr2)

  ai4 <- c(6*PAB^5*PaB,15*PAB^4*PaB^2,20*PAB^3*PaB^3,15*PAB^2*PaB^4,6*PAB*PaB^5)
  pr4 <- ai4/sum(ai4); if(sum(ai4)==0) pr4<-rep(0,5)
  phi4_AB <- sum(c(5,4,3,2,1)*pr4)
  phi4_aB <- sum(c(1,2,3,4,5)*pr4)

  ai5 <- c(6*PAB^5*Pab,30*PAB^4*PAb*PaB,30*PAB^4*PAb*Pab,60*PAB^3*PAb^2*PaB,
           60*PAB^3*PAb^2*Pab,60*PAB^2*PAb^3*PaB,60*PAB^2*PAb^3*Pab,30*PAB*PAb^4*PaB,
           30*PAB*PAb^4*Pab,6*PAb^5*PaB,30*PAB^4*PaB*Pab,60*PAB^3*PAb*PaB^2,
           15*PAB^4*Pab^2,120*PAB^3*PAb*PaB*Pab,90*PAB^2*PAb^2*PaB^2,
           60*PAB^3*PAb*Pab^2,180*PAB^2*PAb^2*PaB*Pab,60*PAB*PAb^3*PaB^2,
           90*PAB^2*PAb^2*Pab^2,120*PAB*PAb^3*PaB*Pab,15*PAb^4*PaB^2,60*PAB*PAb^3*Pab^2,30*PAb^4*PaB*Pab,
           60*PAB^3*PaB^2*Pab,60*PAB^2*PAb*PaB^3,60*PAB^3*PaB*Pab^2,180*PAB^2*PAb*PaB^2*Pab,60*PAB*PAb^2*PaB^3,
           20*PAB^3*Pab^3,180*PAB^2*PAb*PaB*Pab^2,180*PAB*PAb^2*PaB^2*Pab,20*PAb^3*PaB^3,
           60*PAB^2*PAb*Pab^3,180*PAB*PAb^2*PaB*Pab^2,60*PAb^3*PaB^2*Pab,60*PAB*PAb^2*Pab^3,60*PAb^3*PaB*Pab^2,
           60*PAB^2*PaB^3*Pab,30*PAB*PAb*PaB^4,90*PAB^2*PaB^2*Pab^2,120*PAB*PAb*PaB^3*Pab,15*PAb^2*PaB^4,
           60*PAB^2*PaB*Pab^3,180*PAB*PAb*PaB^2*Pab^2,60*PAb^2*PaB^3*Pab,
           15*PAB^2*Pab^4,120*PAB*PAb*PaB*Pab^3,90*PAb^2*PaB^2*Pab^2,30*PAB*PAb*Pab^4,60*PAb^2*PaB*Pab^3,
           30*PAB*PaB^4*Pab,6*PAb*PaB^5,60*PAB*PaB^3*Pab^2,30*PAb*PaB^4*Pab,
           60*PAB*PaB^2*Pab^3,60*PAb*PaB^3*Pab^2,30*PAB*PaB*Pab^4,60*PAb*PaB^2*Pab^3,
           6*PAB*Pab^5,30*PAb*PaB*Pab^4)
  pr5 <- ai5/sum(ai5); if(sum(ai5)==0) pr5<-rep(0,63) # Safety check

  iAB <- c(5,4,4,3,3,2,2,1,1,0,4,3,4,3,2,3,2,1,2,1,0,1,0,3,2,3,2,1,3,2,1,0,2,1,0,1,
           0,2,1,2,1,0,2,1,0,2,1,0,1,0,1,0,1,0,1,0,1,0,1,0)
  iAb <- c(0,1,1,2,2,3,3,4,4,5,0,1,0,1,2,1,2,3,2,3,4,3,4,0,1,0,1,2,0,1,2,3,1,2,3,2,
           3,0,1,0,1,2,0,1,2,0,1,2,1,2,0,1,0,1,0,1,0,1,0,1)
  iaB <- c(0,1,0,1,0,1,0,1,0,1,1,2,0,1,2,0,1,2,0,1,2,0,1,2,3,1,2,3,0,1,2,3,0,1,2,0,
           1,3,4,2,3,4,1,2,3,0,1,2,0,1,4,5,3,4,2,3,1,2,0,1)
  iab <- c(1,0,1,0,1,0,1,0,1,0,1,0,2,1,0,2,1,0,2,1,0,2,1,1,0,2,1,0,3,2,1,0,3,2,1,3,
           2,1,0,2,1,0,3,2,1,4,3,2,4,3,1,0,2,1,3,2,4,3,5,4)

  phi5_AB <- sum(iAB*pr5)
  phi5_Ab <- sum(iAb*pr5)
  phi5_aB <- sum(iaB*pr5)
  phi5_ab <- sum(iab*pr5)

  ai6 <- c(6*PAb^5*Pab,15*PAb^4*Pab^2,20*PAb^3*Pab^3,15*PAb^2*Pab^4,6*PAb*Pab^5)
  pr6 <- ai6/sum(ai6); if(sum(ai6)==0) pr6<-rep(0,5)
  phi6_Ab <- sum(c(5,4,3,2,1)*pr6)
  phi6_ab <- sum(c(1,2,3,4,5)*pr6)

  ai8 <- c(6*PaB^5*Pab^1,15*PaB^4*Pab^2,20*PaB^3*Pab^3,15*PaB^2*Pab^4,6*PaB*Pab^5)
  pr8 <- ai8/sum(ai8); if(sum(ai8)==0) pr8<-rep(0,5)
  phi8_aB <- sum(c(5,4,3,2,1)*pr8)
  phi8_ab <- sum(c(1,2,3,4,5)*pr8)

  epAB <- c(6,phi2_AB,0,phi4_AB,phi5_AB,0,0,0,0)
  epAb <- c(0,phi2_Ab,6,0,phi5_Ab,phi6_Ab,0,0,0)
  epaB <- c(0,0,0,phi4_aB,phi5_aB,0,6,phi8_aB,0)
  epab <- c(0,0,0,0,phi5_ab,phi6_ab,0,phi8_ab,6)

  return(rbind(epAB,epAb,epaB,epab))
}

auto6_est_LDp <- function(dose){
  nts <- rep(0,9)
  nt <- table(dose)
  nts[as.numeric(names(nt))] <- as.numeric(nt)
  nn <- sum(nts)

  res4 <- EM_C6p(nts=nts)

  #DA2,DB2,DA3,DB3
  pAAA <- sum(res4[1:4])
  pAAa <- sum(res4[5:8])
  pAaa <- sum(res4[9:12])
  paaa <- sum(res4[13:16])

  pBBB <- sum(res4[c(1,5,9,13)])
  pBBb <- sum(res4[c(2,6,10,14)])
  pBbb <- sum(res4[c(3,7,11,15)])
  pbbb <- sum(res4[c(4,8,12,16)])

  pA <- pAAA+pAAa*2/3+pAaa/3;pa <- paaa+pAaa*2/3+pAAa/3;
  pB <- pBBB+pBBb*2/3+pBbb/3;pb <- pbbb+pBbb*2/3+pBBb/3;

  rr1 <- solveD(pA=pA,pB=pB,pv=res4)

  estp <- c(pA,pB,rr1)
  gamA <- c(pAAA,pAAa,pAaa,paaa)
  names(gamA) <- c("AAA","AAa","Aaa","aaa")
  gamB <- c(pBBB,pBBb,pBbb,pbbb)
  names(gamB) <- c("AAA","AAa","Aaa","aaa")
  return(list(estp=estp,gamA=gamA,gamB=gamB))
}


EM_C6p <- function(nts){



  p4 <- rep(1/16,16)
  nn <- sum(nts)
  iter <- 0
  nts11 <- matrix(rep(nts,16),nrow=16,byrow=T)
  H1 <- exM6p()
  while(1){

    p41 <- p4
    phi <- ZIfp(p=p4)
    emm <- cm6(HH=H1,PHI=phi)
    p4 <- rowSums(emm*nts11)/(2*nn)

    iter <- iter + 1
    if(max(abs(p4-p41))<1e-8)
      break
  }
  return(p4)
}


exM6p <- function(){

  index <- c("PAAABBB","PAAABBb","PAAABbb","PAAAbbb","PAAaBBB","PAAaBBb","PAAaBbb","PAAabbb",
             "PAaaBBB","PAaaBBb","PAaaBbb","PAaabbb","PaaaBBB","PaaaBBb","PaaaBbb","Paaabbb")

  MM <- matrix(c("PAAABBB-PAAABBB",rep(0,99),"PAAABBB-PAAABBb","PAAABBB-PAAABbb","PAAABBb-PAAABBb",
                 "PAAABBB-PAAAbbb","PAAABBb-PAAABbb","PAAABBb-PAAAbbb","PAAABbb-PAAABbb",
                 "PAAABbb-PAAAbbb",rep(0,92),"PAAAbbb-PAAAbbb",rep(0,99),
                 "PAAABBB-PAAaBBB","PAAABBB-PAaaBBB","PAAaBBB-PAAaBBB","PAAABBB-PaaaBBB","PAAaBBB-PAaaBBB",
                 "PAaaBBB-PAaaBBB","PAAaBBB-PaaaBBB","PAaaBBB-PaaaBBB",rep(0,92),
                 "PAAABBB-PAAaBBb","PAAaBBB-PAAABBb","PAAABBB-PAAaBbb","PAAaBBB-PAAABbb","PAAABBb-PAAaBBb",
                 "PAAABBB-PAAabbb","PAAaBBB-PAAAbbb","PAAABBb-PAAaBbb","PAAaBBb-PAAABbb","PAAABBb-PAAabbb",
                 "PAAaBBb-PAAAbbb","PAAABbb-PAAaBbb","PAAABbb-PAAabbb","PAAaBbb-PAAAbbb","PAAABBB-PAaaBBb",
                 "PAAaBBB-PAAaBBb","PAaaBBB-PAAABBb","PAAaBBb-PAAaBBb","PAAABBB-PAaaBbb","PAAaBBB-PAAaBbb",
                 "PAaaBBB-PAAABbb","PAAABBb-PAaaBBb","PAAABBB-PAaabbb","PAAaBBB-PAAabbb","PAaaBBB-PAAAbbb",
                 "PAAABBb-PAaaBbb","PAAaBBb-PAAaBbb","PAaaBBb-PAAABbb","PAAaBbb-PAAaBbb","PAAABBb-PAaabbb","PAAaBBb-PAAabbb",
                 "PAaaBBb-PAAAbbb","PAAABbb-PAaaBbb","PAAABbb-PAaabbb","PAaaBbb-PAAAbbb","PAAaBbb-PAAabbb",
                 "PAAABBB-PaaaBBb","PAAaBBB-PAaaBBb","PAaaBBB-PAAaBBb","PaaaBBB-PAAABBb","PAAABBB-PaaaBbb","PAAaBBB-PAaaBbb",
                 "PAaaBBB-PAAaBbb","PaaaBBB-PAAABbb","PAAABBb-PaaaBBb","PAAaBBb-PAaaBBb","PAAABBB-Paaabbb","PAAaBBB-PAaabbb",
                 "PAaaBBB-PAAabbb","PaaaBBB-PAAAbbb","PAAABBb-PaaaBbb","PAAaBBb-PAaaBbb","PAaaBBb-PAAaBbb","PaaaBBb-PAAABbb",
                 "PAAABBb-Paaabbb","PAAaBBb-PAaabbb","PAaaBBb-PAAabbb","PaaaBBb-PAAAbbb","PAAABbb-PaaaBbb","PAAaBbb-PAaaBbb",
                 "PAAABbb-Paaabbb","PAAaBbb-PAaabbb","PAaaBbb-PAAabbb","PaaaBbb-PAAAbbb","PAAaBBB-PaaaBBb","PAaaBBB-PAaaBBb",
                 "PaaaBBB-PAAaBBb","PAaaBBb-PAaaBBb","PAAaBBB-PaaaBbb","PAaaBBB-PAaaBbb","PaaaBBB-PAAaBbb","PAAaBBb-PaaaBBb",
                 "PAAaBBB-Paaabbb","PAaaBBB-PAaabbb","PaaaBBB-PAAabbb","PAAaBBb-PaaaBbb","PAaaBBb-PAaaBbb","PaaaBBb-PAAaBbb",
                 "PAaaBbb-PAaaBbb","PAAaBBb-Paaabbb","PAaaBBb-PAaabbb","PaaaBBb-PAAabbb","PAAaBbb-PaaaBbb",
                 "PAAaBbb-Paaabbb","PAaaBbb-PAaabbb","PaaaBbb-PAAabbb","PAaaBBB-PaaaBBb","PaaaBBB-PAaaBBb",
                 "PAaaBBB-PaaaBbb","PaaaBBB-PAaaBbb","PAaaBBb-PaaaBBb",
                 "PAaaBBB-Paaabbb","PaaaBBB-PAaabbb","PAaaBBb-PaaaBbb","PaaaBBb-PAaaBbb",
                 "PAaaBBb-Paaabbb","PaaaBBb-PAaabbb","PAaaBbb-PaaaBbb","PAaaBbb-Paaabbb","PaaaBbb-PAaabbb",
                 "PAAAbbb-PAAabbb","PAAabbb-PAAabbb","PAAAbbb-PAaabbb","PAAAbbb-Paaabbb","PAAabbb-PAaabbb",
                 "PAaabbb-PAaabbb","PAAabbb-Paaabbb", "PAaabbb-Paaabbb",rep(0,92),
                 "PaaaBBB-PaaaBBB",rep(0,99),"PaaaBBB-PaaaBBb","PaaaBBb-PaaaBBb","PaaaBBB-PaaaBbb",
                 "PaaaBBB-Paaabbb","PaaaBBb-PaaaBbb","PaaaBbb-PaaaBbb","PaaaBBb-Paaabbb","PaaaBbb-Paaabbb",
                 rep(0,92),"Paaabbb-Paaabbb",rep(0,99)),ncol=100,byrow=T)

  n1 <- dim(MM)[1];n2 <- dim(MM)[2];


  MIL <- list()
  ii <- 1
  for(i in index){
    MI <- matrix(0,nrow=n1,ncol=n2)
    for(j in 1:n1){
      for(k in 1:n2){
        MI[j,k] <- length(which(strsplit(MM[j,k],split="-")[[1]]==i))
      }
    }
    MIL[[ii]] <- MI
    #print(sum(MI))
    ii <- ii + 1
  }
  return(MIL)
}

ZIfp <- function(p){
  PAAABBB <- p[1]
  PAAABBb <- p[2]
  PAAABbb <- p[3]
  PAAAbbb <- p[4]
  PAAaBBB <- p[5]
  PAAaBBb <- p[6]
  PAAaBbb <- p[7]
  PAAabbb <- p[8]
  PAaaBBB <- p[9]
  PAaaBBb <- p[10]
  PAaaBbb <- p[11]
  PAaabbb <- p[12]
  PaaaBBB <- p[13]
  PaaaBBb <- p[14]
  PaaaBbb <- p[15]
  Paaabbb <- p[16]

  t1 <- c(PAAABBB^2)
  t2 <- c(2*PAAABBB*PAAABBb,2*PAAABBB*PAAABbb,(PAAABBb)^2,
          2*PAAABBB*PAAAbbb,2*PAAABBb*PAAABbb,2*PAAABBb*PAAAbbb,(PAAABbb)^2,
          2*PAAABbb*PAAAbbb)
  t3 <- c(PAAAbbb^2)
  t4 <- c(2*PAAABBB*PAAaBBB,2*PAAABBB*PAaaBBB,(PAAaBBB)^2,2*PAAABBB*PaaaBBB,2*PAAaBBB*PAaaBBB,
          (PAaaBBB)^2,2*PAAaBBB*PaaaBBB,2*PAaaBBB*PaaaBBB)
  t5 <- c(2*PAAABBB*PAAaBBb,2*PAAaBBB*PAAABBb,2*PAAABBB*PAAaBbb,2*PAAaBBB*PAAABbb,2*PAAABBb*PAAaBBb,
          2*PAAABBB*PAAabbb,2*PAAaBBB*PAAAbbb,2*PAAABBb*PAAaBbb,2*PAAaBBb*PAAABbb,2*PAAABBb*PAAabbb,
          2*PAAaBBb*PAAAbbb,2*PAAABbb*PAAaBbb,2*PAAABbb*PAAabbb,2*PAAaBbb*PAAAbbb,2*PAAABBB*PAaaBBb,
          2*PAAaBBB*PAAaBBb,2*PAaaBBB*PAAABBb,(PAAaBBb)^2,2*PAAABBB*PAaaBbb,2*PAAaBBB*PAAaBbb,
          2*PAaaBBB*PAAABbb,2*PAAABBb*PAaaBBb,2*PAAABBB*PAaabbb,2*PAAaBBB*PAAabbb,2*PAaaBBB*PAAAbbb,
          2*PAAABBb*PAaaBbb,2*PAAaBBb*PAAaBbb,2*PAaaBBb*PAAABbb,(PAAaBbb)^2,2*PAAABBb*PAaabbb,2*PAAaBBb*PAAabbb,
          2*PAaaBBb*PAAAbbb,2*PAAABbb*PAaaBbb,2*PAAABbb*PAaabbb,2*PAaaBbb*PAAAbbb,2*PAAaBbb*PAAabbb,
          2*PAAABBB*PaaaBBb,2*PAAaBBB*PAaaBBb,2*PAaaBBB*PAAaBBb,2*PaaaBBB*PAAABBb,2*PAAABBB*PaaaBbb,2*PAAaBBB*PAaaBbb,
          2*PAaaBBB*PAAaBbb,2*PaaaBBB*PAAABbb,2*PAAABBb*PaaaBBb,2*PAAaBBb*PAaaBBb,2*PAAABBB*Paaabbb,2*PAAaBBB*PAaabbb,
          2*PAaaBBB*PAAabbb,2*PaaaBBB*PAAAbbb,2*PAAABBb*PaaaBbb,2*PAAaBBb*PAaaBbb,2*PAaaBBb*PAAaBbb,2*PaaaBBb*PAAABbb,
          2*PAAABBb*Paaabbb,2*PAAaBBb*PAaabbb,2*PAaaBBb*PAAabbb,2*PaaaBBb*PAAAbbb,2*PAAABbb*PaaaBbb,2*PAAaBbb*PAaaBbb,
          2*PAAABbb*Paaabbb,2*PAAaBbb*PAaabbb,2*PAaaBbb*PAAabbb,2*PaaaBbb*PAAAbbb,2*PAAaBBB*PaaaBBb,2*PAaaBBB*PAaaBBb,
          2*PaaaBBB*PAAaBBb,(PAaaBBb)^2,2*PAAaBBB*PaaaBbb,2*PAaaBBB*PAaaBbb,2*PaaaBBB*PAAaBbb,2*PAAaBBb*PaaaBBb,
          2*PAAaBBB*Paaabbb,2*PAaaBBB*PAaabbb,2*PaaaBBB*PAAabbb,2*PAAaBBb*PaaaBbb,2*PAaaBBb*PAaaBbb,2*PaaaBBb*PAAaBbb,
          (PAaaBbb)^2,2*PAAaBBb*Paaabbb,2*PAaaBBb*PAaabbb,2*PaaaBBb*PAAabbb,2*PAAaBbb*PaaaBbb,
          2*PAAaBbb*Paaabbb,2*PAaaBbb*PAaabbb,2*PaaaBbb*PAAabbb,2*PAaaBBB*PaaaBBb,2*PaaaBBB*PAaaBBb,
          2*PAaaBBB*PaaaBbb,2*PaaaBBB*PAaaBbb,2*PAaaBBb*PaaaBBb,
          2*PAaaBBB*Paaabbb,2*PaaaBBB*PAaabbb,2*PAaaBBb*PaaaBbb,2*PaaaBBb*PAaaBbb,
          2*PAaaBBb*Paaabbb,2*PaaaBBb*PAaabbb,2*PAaaBbb*PaaaBbb,2*PAaaBbb*Paaabbb,2*PaaaBbb*PAaabbb)
  t6 <- c(2*PAAAbbb*PAAabbb,(PAAabbb)^2,2*PAAAbbb*PAaabbb,2*PAAAbbb*Paaabbb,2*PAAabbb*PAaabbb,
          (PAaabbb)^2,2*PAAabbb*Paaabbb, 2*PAaabbb*Paaabbb)
  t7 <- c((PaaaBBB)^2)
  t8 <-  c(2*PaaaBBB*PaaaBBb,(PaaaBBb)^2,2*PaaaBBB*PaaaBbb,
           2*PaaaBBB*Paaabbb,2*PaaaBBb*PaaaBbb,(PaaaBbb)^2,2*PaaaBBb*Paaabbb,2*PaaaBbb*Paaabbb)
  t9 <-  c((Paaabbb)^2)

  zp <- matrix(c(t1,rep(0,100-length(t1)),t2,rep(0,100-length(t2)),t3,rep(0,100-length(t3)),
                 t4,rep(0,100-length(t4)),t5,t6,rep(0,100-length(t6)),
                 t7,rep(0,100-length(t7)),t8,rep(0,100-length(t8)),t9,rep(0,100-length(t9))),ncol=100,byrow=T)

  zs <-  cbind(zp,rowSums(zp))

  conp <- t(apply(zs,1,function(x){
    if(x[101]==0){
      x[101] <- 1e-200
    }
    x[1:100]/x[101]

  }))

  conp
}
