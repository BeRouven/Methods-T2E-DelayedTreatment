################################################################################
####        Costume functions needed for analysis of generated data         ####
################################################################################

library(tidyverse)

#------------- Function to calculate weighted logrank statistics ---------------
#---------------------------- and their covariance------------------------------
weightedLR <- function(data, 
                       vars = list(time= "time", event = "event", group = "group"),
                       params = list(FH = list(rho = c(1,1,0,0), gamma=c(0,1,0,1)),
                                     Cheng = 0.5,
                                     V0 = threshold,
                                     MERT = list(t.low, t.up),
                                     gen.lin = list(t.low, t.up),
                                     threshold = threshold,
                                     logit = list(w.low = 0.1, t.low, t.up),
                                     par.group = threshold)){
  ##Necessary packages:
  # survival, YPmodel
  #
  ##Parameters:
  # data:     input dataset of survival data
  # vars:     names of the time, event and group variables in the input dataset
  # params:   parameters for the analysis of the methods
  #
  ##Output:
  # res.all:  List with the following entries
  #   FH:         (std.) test statistic G (Z) and covariance of the specified FH tests
  #   mMaxCombo:  (std.) test statistic G (Z) and covariance of the LR, FH(0,1), FH(1,0) and Cheng test
  #   mScore:     (std.) test statistic G (Z) and covariance of the LR and NA weighted LR test
  #   mZm3:       (std.) test statistic G (Z) and covariance of the LR, G(0,1) and the modified G(1,0)
  #   AWLRT:      (std.) test statistic G (Z) and covariance of the LR, YP1 and YP2
  #   V0:         (std.) test statistic G (Z) and covariance of the LR and threshold test
  #   par.group:  (std.) test statistic G (Z) and variance of the ParGroup test
  #   KS.XXX:     (std.) test statistic G (Z) and variance of the XXX test (LR, FH, GB and Cheng) at each event time
  #   MERT:       (std.) test statistic G (Z) and variance of the MERT test
  #   GB:         (std.) test statistic G (Z) and variance of the GB test
  #   TW:         (std.) test statistic G (Z) and variance of the TW test
  #   gen.lin:    (std.) test statistic G (Z) and variance of the GenLin test
  #   threshold:  (std.) test statistic G (Z) and variance of the Thres test
  #   PP:         (std.) test statistic G (Z) and variance of the PP test
  #   mPP:        (std.) test statistic G (Z) and variance of the mPP test
  #   asymLR:     (std.) test statistic G (Z) and variance of the asymLR test
  #   logit:      (std.) test statistic G (Z) and variance of the Logit test
  #   mlogit:     (std.) test statistic G (Z) and variance of the mLogit test
  
  data.frame(X = data[,vars$time], delta = data[,vars$event], Z1 = as.numeric(as.character(data[,vars$group])),
             Z0 = 1 - as.numeric(as.character(data[,vars$group]))) ->
    data.simu
  
  sample.size <- nrow(data.simu)
  
  fit_tmp<-survival::survfit(survival::Surv(X,delta)~1,data=data.simu, conf.type = "none")
  fit_tmp2<-survival::survfit(survival::Surv(X,delta)~Z1,data=data.simu, conf.type = "none")
  
  n.times0 <- fit_tmp2$strata[1]
  n.times1 <- fit_tmp2$strata[2]
  
  df.tmp.all <- data.frame(time = fit_tmp$time, n.risk = fit_tmp$n.risk, 
                           n.event = fit_tmp$n.event, surv = fit_tmp$surv,
                           Nelson.Aalen = fit_tmp$cumhaz)
  df.tmp0 <- data.frame(time = fit_tmp2$time, n.risk0 = fit_tmp2$n.risk, 
                        n.event0 = fit_tmp2$n.event)[1:n.times0,]
  df.tmp1 <- data.frame(time = fit_tmp2$time, n.risk1 = fit_tmp2$n.risk, 
                        n.event1 = fit_tmp2$n.event)[(n.times0+1):(n.times0+n.times1),]
  
  #Merge together number at risk and number of events overall and by group
  df.tmp.all[, c("n.risk0", "n.event0")] <- df.tmp0[match(df.tmp.all$time, df.tmp0$time), c("n.risk0", "n.event0")]
  df.tmp.all[, c("n.risk1", "n.event1")] <- df.tmp1[match(df.tmp.all$time, df.tmp1$time), c("n.risk1", "n.event1")]
  
  df.tmp <- df.tmp.all[df.tmp.all$n.event!=0,]
  df.tmp$n.event0[is.na(df.tmp$n.event0)] <- 0
  df.tmp$n.event1[is.na(df.tmp$n.event1)] <- 0
  df.tmp$n.risk0[is.na(df.tmp$n.risk0)] <- df.tmp$n.risk[is.na(df.tmp$n.risk0)] - df.tmp$n.risk1[is.na(df.tmp$n.risk0)]
  df.tmp$n.risk1[is.na(df.tmp$n.risk1)] <- df.tmp$n.risk[is.na(df.tmp$n.risk1)] - df.tmp$n.risk0[is.na(df.tmp$n.risk1)]
  df.tmp$rank <- rank(df.tmp$time)
  
  n.times <- length(df.tmp$time)
  
  #IF the overall number at risk at the last timepoint is 1, this time-point
  #does NOT contribute to the estimation of variance since the number at risk
  #in one of the groups must then be 0.
  stat <- (df.tmp$n.event0-df.tmp$n.event*df.tmp$n.risk0/df.tmp$n.risk)
  var <- ifelse(df.tmp$n.risk==1, 0,
                (df.tmp$n.risk0*df.tmp$n.risk1*df.tmp$n.event*(df.tmp$n.risk-df.tmp$n.event)/(df.tmp$n.risk^2*(df.tmp$n.risk-1)))
  )
  res.all <- list()
  
  # Fleming-Harrington ------------------------------------------------------
  rho <- params$FH$rho
  gamma <- params$FH$gamma
  if(length(rho) != length(gamma)){
    print("Error: Weight vectors must be of the same length")
  } else {
    #Make it left continuous
    surv_t <- c(1,df.tmp$surv)[1:n.times]
    
    G <- c()
    Cov <- matrix(data=NA, nrow=length(rho), ncol = length(gamma))
    for (i in 1:length(rho)){
      weight1<-(surv_t)^rho[i]*(1-surv_t)^gamma[i]
      G[i] <- sum(weight1*stat)
      for (j in 1:length(gamma)){
        weight2 <- (surv_t)^rho[j]*(1-surv_t)^gamma[j]
        Cov[i,j] <- sum(weight1*weight2*var)
      }
    }
    Z <- G/sqrt(diag(Cov))
    res.all$FH <- cbind(G, Z, Cov)
    rownames(res.all$FH) <- paste0(rho, gamma)
    colnames(res.all$FH) <- c("G", "Z", paste0("Cov", rho, gamma))
    
  }
  # mMaxCombo -------------------------------------------------------------------
  #Make it left continuous
  surv_t <- c(1,df.tmp$surv)[1:n.times]
  
  weight1 <- 1
  weight2 <- surv_t
  weight3 <- 1-surv_t
  #Define weight proposed by Cheng
  theta <- 0.5
  weight4 <- sapply(surv_t, function(x) 1*(x <= theta)*(x-theta)/theta + 1*(x>theta)*(x-theta)/(1-theta))
  
  G <- c(sum(weight1*stat),
         sum(weight2*stat),
         sum(weight3*stat),
         sum(weight4*stat))
  
  Cov <- matrix(data=NA, nrow=4, ncol = 4)
  Cov[1,1] <- sum(weight1*weight1*var)
  Cov[2,2] <- sum(weight2*weight2*var)
  Cov[3,3] <- sum(weight3*weight3*var)
  Cov[4,4] <- sum(weight4*weight4*var)
  Cov[1,2] <- Cov[2,1] <- sum(weight1*weight2*var)
  Cov[1,3] <- Cov[3,1] <- sum(weight1*weight3*var)
  Cov[2,3] <- Cov[3,2] <- sum(weight2*weight3*var)
  Cov[1,4] <- Cov[4,1] <- sum(weight1*weight4*var)
  Cov[2,4] <- Cov[4,2] <- sum(weight2*weight4*var)
  Cov[3,4] <- Cov[4,3] <- sum(weight3*weight4*var)
  
  Z <- G/sqrt(diag(Cov))
  res.all$mMaxCombo <- cbind(G, Z, Cov)
  rownames(res.all$mMaxCombo) <- c("LR", "FH(1,0)", "FH(0,1)", "Cheng")
  colnames(res.all$mMaxCombo) <- c("G", "Z", paste0("Cov.", c("LR", "FH(1,0)", "FH(0,1)", "Cheng")))
  
  # Modified Score ----------------------------------------------------------
  #Make it left continuous
  NA_t <- c(0,df.tmp$Nelson.Aalen)[1:n.times]
  
  weight1 <- 1
  weight2 <- log(1+NA_t)
  
  G <- c(sum(weight1*stat),
         sum(weight2*stat))
  
  Cov <- matrix(data=NA, nrow=2, ncol = 2)
  Cov[1,1] <- sum(weight1*weight1*var)
  Cov[2,2] <- sum(weight2*weight2*var)
  Cov[1,2] <- Cov[2,1] <- sum(weight1*weight2*var)
  
  Z <- G/sqrt(diag(Cov))
  res.all$mScore <- cbind(G, Z, Cov)
  rownames(res.all$mScore) <- c("LR", "NA")
  colnames(res.all$mScore) <- c("G", "Z", paste0("Cov.", c("LR", "NA")))
  
  # modified Zm3 -------------------------------------------------------
  #Make it left continuous
  surv_t <- c(1,df.tmp$surv)[1:n.times]
  
  weight1 <- 1
  weight2 <- 1-surv_t
  weight3 <- pmax(0.001, (surv_t - min(surv_t))/(1-min(surv_t)))
  
  G <- c(sum(weight1*stat),
         sum(weight2*stat),
         sum(weight3*stat))
  
  Cov <- matrix(data=NA, nrow=3, ncol = 3)
  Cov[1,1] <- sum(weight1*weight1*var)
  Cov[2,2] <- sum(weight2*weight2*var)
  Cov[3,3] <- sum(weight3*weight3*var)
  Cov[1,2] <- Cov[2,1] <- sum(weight1*weight2*var)
  Cov[1,3] <- Cov[3,1] <- sum(weight1*weight3*var)
  Cov[2,3] <- Cov[3,2] <- sum(weight2*weight3*var)
  
  Z <- G/sqrt(diag(Cov))
  res.all$mZm3 <- cbind(G, Z, Cov)
  rownames(res.all$mZm3) <- c("LR", "FH(0,1)", "mFH(1,0)")
  colnames(res.all$mZm3) <- c("G", "Z", paste0("Cov.", c("LR", "FH(0,1)", "mFH(1,0)")))
  
  # Yang-Prentice model -----------------------------------------------------
  data.simu.order <- data.simu[order(data.simu$X),]
  data.frame(V1 = data.simu.order$X, V2 = data.simu.order$delta,
             V3 = data.simu.order$Z1, Z0 = data.simu.order$Z0, Y=sample.size:1,
             Y1=cumsum(data.simu.order$Z1[sample.size:1])[sample.size:1],
             Y0=cumsum(data.simu.order$Z0[sample.size:1])[sample.size:1]) ->
    data.rev
  
  YP.stat <- data.rev$V2*(data.rev$V3*data.rev$Y0-data.rev$Z0*data.rev$Y1)/data.rev$Y
  YP.var <- data.rev$V2*data.rev$Y0*data.rev$Y1/data.rev$Y^2
  YPestim <- YPmodel::YPmodel.estimate(data.rev)
  
  beta <- YPestim$beta
  r <- YPestim$r
  res.all$AWLRT <- matrix(data = NA, nrow = 3, ncol = 5)
  rownames(res.all$AWLRT) <- c("LR", "YP1", "YP2")
  colnames(res.all$AWLRT) <- c("G", "Z", "Cov.LR", "Cov.YP1", "Cov.YP2")
  
  weight0 <- 1
  res.all$AWLRT[1,1] <- sum(weight0*YP.stat)
  res.all$AWLRT[1,3] <- sum(weight0^2*YP.var)
  res.all$AWLRT[1,2] <- res.all$AWLRT[1,1]/sqrt(res.all$AWLRT[1,3])
  
  weight1 <- (1 + r)/(exp(-beta[1]) + exp(-beta[2]) * r)
  res.all$AWLRT[2,1] <- sum(weight1*YP.stat)
  res.all$AWLRT[2,4] <- sum(weight1^2*YP.var)
  res.all$AWLRT[2,2] <-res.all$AWLRT[2,1]/sqrt(res.all$AWLRT[2,4])
  res.all$AWLRT[2,3] <- res.all$AWLRT[1,4] <- sum(weight1*YP.var)
  
  weight2 <- 1/weight1
  res.all$AWLRT[3,1] <- sum(weight2*YP.stat)
  res.all$AWLRT[3,5] <- sum(weight2^2*YP.var)
  res.all$AWLRT[3,2] <-res.all$AWLRT[3,1]/sqrt(res.all$AWLRT[3,5])
  res.all$AWLRT[3,3] <- res.all$AWLRT[1,5] <- sum(weight2*YP.var)
  res.all$AWLRT[3,4] <- res.all$AWLRT[2,5] <- sum(YP.var)
  # V0 ----------------------------------------------------------------------
  res.all$V0 <- list()
  
  V0 <- params$V0
  
  for (i in 1: length(V0)){
    #V0 is the sum of the logrank test and the threshold weighted LR test
    res <- matrix(data = NA, nrow = 2, ncol = 5)
    rownames(res) <- c("LR", paste0("threshold=", V0[i]))
    colnames(res) <- c("G", "Z", "Cov.LR", "Cov.threshold", "Z.V0")
    
    weight1 <- 1
    res[1,1] <- sum(weight1*stat)
    res[1,3] <- sum(weight1^2*var)
    res[1,2] <-res[1,1]/sqrt(res[1,3])
    
    weight2 <- 1*(df.tmp$time > V0[i])
    res[2,1] <- sum(weight2*stat)
    #Since the threshold weight squared is the same, the variance equals the covariance with the LR test
    res[2,4] <- res[1,4] <- res[2,3] <- sum(weight2^2*var)
    res[2,2] <- res[2,1]/sqrt(res[2,4])
    
    cor.hat <- sqrt(res[2,4])/sqrt(res[1,3])
    res[1,5] <- (res[1,2] + res[2,2])/sqrt(2*(1+cor.hat))
    
    res.all$V0[[i]] <- res
  } 
  names(res.all$V0) <- paste0("threshold=", V0)
  
  # Partially grouped LR test -----------------------------------------------
  res.all$par.group <- list()
  
  par.group <- params$par.group
  for (i in 1: length(par.group)){
    res <- matrix(NA, nrow=1, ncol=3)
    colnames(res) <- c("G", "Var", "Z")
    
    if (par.group[i] < min(fit_tmp2$time[fit_tmp2$strata[1]], fit_tmp2$time[sum(fit_tmp2$strata)])){
      n.risk <- summary(fit_tmp2, times = 0)$n.risk
      surv_threshold <- summary(fit_tmp2, times = par.group[i])
      
      x <- n.risk*(1-surv_threshold$surv)
      G0 <- x[1]-sum(x)*n.risk[1]/(sum(n.risk))
      V0 <- prod(n.risk)*summary(fit_tmp, times = par.group[i])$std.err^2
      
      #Part of the threshold LR test for the partially grouped LR test
      weight <- 1*(df.tmp$time > par.group[i])
      G1 <- sum(weight*stat)
      V1 <- sum(weight^2*var)
      
      res[1,1] <- G0 + G1
      res[1,2] <- V0 + V1
      res[1,3] <- res[1,1] / sqrt(res[1,2])
    }
    res.all$par.group[[i]] <- res
  } 
  names(res.all$par.group) <- paste0("threshold=", par.group)
  
  # Kolmogorov-Smirnov test -----------------------------------------------
  #LR
  res.all$KS.LR <- matrix(NA, nrow=n.times, ncol=3)
  colnames(res.all$KS.LR) <- c("G", "Var", "Z")
  for (ind in 1:n.times){
    weight.tmp <- rep(1, n.times)*c(rep(1,ind), rep(0, n.times-ind))
    res.all$KS.LR[ind, "G"] <- sum(weight.tmp*stat)
    res.all$KS.LR[ind, "Var"] <- sum(weight.tmp^2*var)
    res.all$KS.LR[ind, "Z"] <- res.all$KS.LR[ind, "G"]/sqrt(res.all$KS.LR[ind, "Var"])
  }
  #GB
  res.all$KS.GB <- matrix(NA, nrow=n.times, ncol=3)
  colnames(res.all$KS.GB) <- c("G", "Var", "Z")
  for (ind in 1:n.times){
    weight.tmp <- df.tmp$n.risk*c(rep(1,ind), rep(0, n.times-ind))
    res.all$KS.GB[ind, "G"] <- sum(weight.tmp*stat)
    res.all$KS.GB[ind, "Var"] <- sum(weight.tmp^2*var)
    res.all$KS.GB[ind, "Z"] <- res.all$KS.GB[ind, "G"]/sqrt(res.all$KS.GB[ind, "Var"])
  }
  
  #Cheng
  res.all$KS.Cheng <- matrix(NA, nrow=n.times, ncol=3)
  colnames(res.all$KS.Cheng) <- c("G", "Var", "Z")
  for (ind in 1:n.times){
    weight.tmp <- (2*c(1,df.tmp$surv)[1:n.times]-1)*c(rep(1,ind), rep(0, n.times-ind))
    res.all$KS.Cheng[ind, "G"] <- sum(weight.tmp*stat)
    res.all$KS.Cheng[ind, "Var"] <- sum(weight.tmp^2*var)
    res.all$KS.Cheng[ind, "Z"] <- res.all$KS.Cheng[ind, "G"]/sqrt(res.all$KS.Cheng[ind, "Var"])
  }
  
  #FH
  res.all$KS.FH <- matrix(NA, nrow=n.times, ncol=3)
  colnames(res.all$KS.FH) <- c("G", "Var", "Z")
  for (ind in 1:n.times){
    weight.tmp <- c(1,df.tmp$surv)[1:n.times]*c(rep(1,ind), rep(0, n.times-ind))
    res.all$KS.FH[ind, "G"] <- sum(weight.tmp*stat)
    res.all$KS.FH[ind, "Var"] <- sum(weight.tmp^2*var)
    res.all$KS.FH[ind, "Z"] <- res.all$KS.FH[ind, "G"]/sqrt(res.all$KS.FH[ind, "Var"])
  }
  # MERT --------------------------------------------------------------------
  res.all$MERT <- list()
  
  t1.tilde <- params$MERT$t.low
  t2.tilde <- params$MERT$t.up
  
  for (i in 1:length(t1.tilde)){
    Psi <- 1/sample.size*cumsum(df.tmp$n.event*df.tmp$n.risk0*df.tmp$n.risk1/(df.tmp$n.risk)^2)
    Psi.h.tau <- max(Psi)
    Psi.h.t1.tilde <- ifelse(df.tmp$time[1]>t1.tilde[i],
                             0,
                             Psi[max(which(df.tmp$time<=t1.tilde[i]))])
    Psi.h.t2.tilde <- Psi[max(which(df.tmp$time<=t2.tilde[i]))]
    Psi.rev <- Psi*(df.tmp$time<=t2.tilde[i]) #only use the Psi before t2**
    W.star <- ((Psi.h.tau-Psi.rev)/(Psi.h.tau-Psi.h.t1.tilde))^(-1/2)*(df.tmp$time<=t2.tilde[i]&df.tmp$time>=t1.tilde[i])
    + 2*((Psi.h.tau-Psi.h.t2.tilde)/(Psi.h.tau-Psi.h.t1.tilde))^(-1/2)*(df.tmp$timeX>t2.tilde[i])
    
    res <- matrix(NA, nrow=1, ncol=3)
    colnames(res) <- c("G", "Var", "Z")
    res[1,1] <- sum(W.star*stat)
    res[1,2] <- sum(W.star^2*var)
    res[1,3] <- res[1,1]/sqrt(res[1,2])
    
    res.all$MERT[[i]] <- res
  } 
  names(res.all$MERT) <- paste0("t.low=", t1.tilde, ".t.up=", t2.tilde)
  # Gehan-Breslow -----------------------------------------------------------
  #Calculation of weighted logrank statistic
  res.all$GB <- matrix(NA, nrow=1, ncol=3)
  colnames(res.all$GB) <- c("G", "Var", "Z")
  res.all$GB[1,1] <- sum(df.tmp$n.risk*stat)
  res.all$GB[1,2] <- sum(df.tmp$n.risk^2*var)
  res.all$GB[1,3] <- res.all$GB[1,1]/sqrt(res.all$GB[1,2])
  # Tarone-Ware -------------------------------------------------------------
  res.all$TW <- matrix(NA, nrow=1, ncol=3)
  colnames(res.all$TW) <- c("G", "Var", "Z")
  res.all$TW[1,1] <- sum(sqrt(df.tmp$n.risk)*stat)
  res.all$TW[1,2] <- sum(sqrt(df.tmp$n.risk)^2*var)
  res.all$TW[1,3] <- res.all$TW[1,1]/sqrt(res.all$TW[1,2])
  # Generalized linear lag model --------------------------------------------
  res.all$gen.lin <- list()
  
  t1.tilde <- params$gen.lin$t.low
  t2.tilde <- params$gen.lin$t.up
  for (i in 1:length(t1.tilde)){
    weight <- sapply(df.tmp$time,
                     function(t) (t-t1.tilde[i])/(t2.tilde[i]-t1.tilde[i])*(t<=t2.tilde[i]&t>t1.tilde[i])+(t>t2.tilde[i]))
    
    res <- matrix(NA, nrow=1, ncol=3)
    colnames(res) <- c("G", "Var", "Z")
    res[1,1] <- sum(weight*stat)
    res[1,2] <- sum(weight^2*var)
    res[1,3] <- res[1,1]/sqrt(res[1,2])
    
    res.all$gen.lin[[i]] <- res
  }
  
  names(res.all$gen.lin) <- paste0("t.low=", t1.tilde, ".t.up=", t2.tilde)
  # Linear threshold model --------------------------------------------------
  res.all$threshold <- list()
  
  threshold <- params$threshold
  for (i in 1:length(threshold)){
    weight <- 1*(df.tmp$time > threshold[i])
    
    res <- matrix(NA, nrow=1, ncol=3)
    colnames(res) <- c("G", "Var", "Z")
    res[1,1] <- sum(weight*stat)
    res[1,2] <- sum(weight^2*var)
    res[1,3] <- res[1,1]/sqrt(res[1,2])
    
    res.all$threshold[[i]] <- res
  }
  
  names(res.all$threshold) <- paste0("threshold=", threshold)
  # (Modified) Peto-Peto ----------------------------------------------------
  mKM.PP <- c()
  for (i in 1:n.times){
    if (i == 1){
      mKM.PP[1] <- 1-df.tmp$n.event[1]/(df.tmp$n.risk[1]+1)
    } else {
      mKM.PP[i] <- mKM.PP[(i-1)] * (1-df.tmp$n.event[i]/(df.tmp$n.risk[i]+1))
    }
  }
  
  weight.PP <- sort(mKM.PP,decreasing = TRUE)
  res.all$PP <- matrix(NA, nrow=1, ncol=3)
  colnames(res.all$PP) <- c("G", "Var", "Z")
  res.all$PP[1,1] <- sum(weight.PP*stat)
  res.all$PP[1,2] <- sum(weight.PP^2*var)
  res.all$PP[1,3] <- res.all$PP[1,1]/sqrt(res.all$PP[1,2])
  
  weight.mPP <- sort(mKM.PP * df.tmp$n.risk/(df.tmp$n.risk+1),decreasing = TRUE)
  res.all$mPP <- matrix(NA, nrow=1, ncol=3)
  colnames(res.all$mPP) <- c("G", "Var", "Z")
  res.all$mPP[1,1] <- sum(weight.mPP*stat)
  res.all$mPP[1,2] <- sum(weight.mPP^2*var)
  res.all$mPP[1,3] <- res.all$mPP[1,1]/sqrt(res.all$mPP[1,2])    
  
  # Asymptotic Logrank ------------------------------------------------------
  Eu <- c()
  for (i in 1:n.times){
    if (i == 1){
      Eu[1] <- df.tmp$n.risk[1]/(df.tmp$n.risk[1] + df.tmp$n.event[1])
    } else {
      Eu[i] <- Eu[(i-1)] * df.tmp$n.risk[i]/(df.tmp$n.risk[i] + df.tmp$n.event[i])
    }
  }
  Eu <- sort(Eu,decreasing = TRUE)
  
  weight.asymLR <- 1 + log(-log(Eu))
  res.all$asymLR <- matrix(NA, nrow=1, ncol=3)
  colnames(res.all$asymLR) <- c("G", "Var", "Z")
  res.all$asymLR[1,1] <- sum(weight.asymLR*stat)
  res.all$asymLR[1,2] <- sum(weight.asymLR^2*var)
  res.all$asymLR[1,3] <- res.all$asymLR[1,1]/sqrt(res.all$asymLR[1,2])   
  
  # (Modified) logit --------------------------------------------------------
  #Define logit weight function
  Yu.logit <- Vectorize(function(x, a, tau) {exp(a*(x-tau))/(1+exp(a*(x-tau)))}, "x")
  
  res.all$logit <- list()
  res.all$mlogit <- list()
  
  t.low <- params$logit$t.low
  t.up <- params$logit$t.up
  w.low <- params$logit$w.low
  for (i in 1:length(t.low)){
    #Define parameters of logit function based on the transition period and the starting weight      
    tau <- (t.low[i] + t.up[i])/2
    a <- log(w.low/(1-w.low))/(t.low[i]-tau)
    
    #Logit model
    weight <- Yu.logit(x=df.tmp$time, a = a, tau = tau)
    res <- matrix(NA, nrow=1, ncol=3)
    colnames(res) <- c("G", "Var", "Z")
    res[1,1] <- sum(weight*stat)
    res[1,2] <- sum(weight^2*var)
    res[1,3] <- res[1,1]/sqrt(res[1,2])
    
    res.all$logit[[i]] <- res
    
    #Modified logit model
    weight <- (Yu.logit(x=df.tmp$time, a = a, tau = tau) - Yu.logit(x=t.low[i], a = a, tau = tau))/
      (Yu.logit(x=t.up[i], a = a, tau = tau) - Yu.logit(x=t.low[i], a = a, tau = tau))
    weight[weight<0] <- 0
    weight[weight>1] <- 1
    res <- matrix(NA, nrow=1, ncol=3)
    colnames(res) <- c("G", "Var", "Z")
    res[1,1] <- sum(weight*stat)
    res[1,2] <- sum(weight^2*var)
    res[1,3] <- res[1,1]/sqrt(res[1,2])
    
    res.all$mlogit[[i]] <- res
  }
  
  names(res.all$logit) <- paste0("t.low=", t.low, ".t.up=", t.up)
  names(res.all$mlogit) <- paste0("t.low=", t.low, ".t.up=", t.up)
  
  return(res.all)
}


#------------- Function to calculate ABC statistic -----------------------------
#function is taken from RBT4TCSC and the error due to missing group sizes fixed
LinStatABC.cor <- function (sample1, sample2) 
{
  sur.fit1 = survfit(Surv(sample1[, 1], sample1[, 2]) ~ 1)
  sur.fit2 = survfit(Surv(sample2[, 1], sample2[, 2]) ~ 1)
  sur.est1 = stepfun(sur.fit1$time, c(1, sur.fit1$surv))
  sur.est2 = stepfun(sur.fit2$time, c(1, sur.fit2$surv))
  pooledSample = rbind(sample1, sample2)
  event.time = sort(pooledSample[pooledSample[, 2] == 1, 1])
  sta1 = subset(sample1, sample1[, 1] == max(sample1[, 1]))[2]
  sta2 = subset(sample2, sample2[, 1] == max(sample2[, 1]))[2]
  if (sta1 == 0 && sta2 == 0) {
    tau = min(max(sample1[, 1]), max(sample2[, 1]))
  }
  else {
    if (sta1 == 1 && sta2 == 1) {
      tau = max(max(sample1[, 1]), max(sample2[, 1]))
    }
    else {
      tau = max(max(sample1[, 1]) * (1 - sta1), max(sample2[, 
                                                            1]) * (1 - sta2))
    }
  }
  gapTime = diff(c(event.time[event.time < tau], tau))
  Lin_ABC = sum(abs(sur.est1(event.time[event.time < tau]) - 
                      sur.est2(event.time[event.time < tau])) * gapTime)
  GW.std.error = function(samsize, surv, risk, event) {
    sigma = NULL
    for (i in 1:samsize) {
      S = 0
      for (j in 1:i) {
        S = S + event[j]/(risk[j] * (risk[j] - event[j]))
      }
      sigma[i] = surv[i] * sqrt(S)
    }
    return(sigma)
  }
  Sigma1 = GW.std.error(sur.fit1$n, sur.fit1$surv, sur.fit1$n.risk, 
                        sur.fit1$n.event)
  Sigma2 = GW.std.error(sur.fit2$n, sur.fit2$surv, sur.fit2$n.risk, 
                        sur.fit2$n.event)
  sigma.est1 = stepfun(sort(sample1[, 1]), c(0, Sigma1))
  sigma.est2 = stepfun(sort(sample2[, 1]), c(0, Sigma2))
  Sigma1 = sigma.est1(event.time[event.time < tau])
  Sigma2 = sigma.est2(event.time[event.time < tau])
  Sigma1[is.na(Sigma1)] <- 0
  Sigma2[is.na(Sigma2)] <- 0
  E_Lin_ABC = sum(sqrt(2/pi * (Sigma1^2 + Sigma2^2)) * gapTime)
  rho = 1/2
  Cov = 0
  L = length(gapTime)
  if (L >= 2) {
    for (j in 2:L) {
      for (i in 1:(j - 1)) {
        Cov = Cov + 2 * rho * (gapTime[j]) * (gapTime[i]) * 
          (1 - 2/pi) * sqrt((Sigma1[j]^2 + Sigma2[j]^2) * 
                              (Sigma1[i]^2 + Sigma2[i]^2))
      }
    }
  }
  else {
    Cov = 0
  }
  V_Lin_ABC = (1 - 2/pi) * sum((Sigma1^2 + Sigma2^2) * (gapTime^2)) + 
    Cov
  stat_Lin = (Lin_ABC - E_Lin_ABC)/sqrt(V_Lin_ABC)
  return(stat_Lin)
}

#---------------- Function to calculate p-values of all tests ------------------
MethodsSurvival <- function(data, 
                            vars = list(time="obs.time", 
                                        event = "event", 
                                        group = "group"),
                            params = list(AFTdist = c("weibull",
                                                      "exponential", 
                                                      "gaussian", 
                                                      "logistic",
                                                      "lognormal", 
                                                      "loglogistic"),
                                          FH = list(rho = c(0,0,1,1,-1,0,0,0.5),
                                                    gamma = c(0,1,0,1,0,0.5,2,0.5)),
                                          mWLR.Delay,
                                          RMST.cutoff, 
                                          milestone.cutoff,
                                          landmark.cutoff,
                                          gen.lin = list(t.low, t.up),
                                          PW.exp.lag = list(lag, CP),
                                          threshold, 
                                          ABC,
                                          PW.exp,
                                          V0, 
                                          MERT = list(t.low, t.up),
                                          Cheng = 0.5,
                                          par.group,
                                          logit = list(w.low = 0.1, t.low, t.up))){
  ##Necessary packages:
  # tidyverse, survival, modestWLRT, timereg, survRM2, mnormt, nphsim, coxphw,
  # CauchyCP, flexsurv, RBT4TCSC
  #
  ##Parameters:
  # data:     input dataset of survival data
  # vars:     names of the time, event and group variables in the input dataset
  # params:   parameters for the analysis of the methods
  #
  ##Output:
  # res:      dataframe with name and p-value of each test
  if(length(params$gen.lin$t.low) != length(params$gen.lin$t.up)){
    print(paste("ERROR: The length of the parameters for general linear lag model
                must be the same. First parameter has length", 
                length(params$gen.lin$t.low),
                "and second parameter has length", 
                length(params$gen.lin$t.up)))
  } else if(any(params$gen.lin$t.low>=params$gen.lin$t.up)){
    print("ERROR: For the general linear lag model the time when the hazards
                start to diverge must be lower than the time when the full effect is reached.")
    print(paste("The first parameter specified is", 
                params$gen.lin$t.low,
                "and the second parameter is", 
                params$gen.lin$t.up))
  } else if(length(params$MERT$t.low) != length(params$MERT$t.up)){
    print(paste("ERROR: The length of the parameters for MERT must be the same.
                First parameter has length", 
                length(params$MERT$t.low),
                "and second parameter has length", 
                length(params$MERT$t.up)))
  } else if(any(params$MERT$t.low>=params$MERT$t.up)){
    print("ERROR: For MERT the time when the hazards start to diverge must
    be lower than the time when the full effect is reached.")
    print(paste("The first parameter specified is", 
                params$MERT$t.low,
                "and the second parameter is", 
                params$MERT$t.up))
  } else{
    df.work <- data.frame(time = data[,vars$time],
                          event = data[,vars$event],
                          group = data[,vars$group])
    #If an event at time 0 is observed, this can cause problems with the survSplit
    #function within the CauchyCP function. Hence the observed times are set to
    #1% of the lowest non-zero event time
    if (any(df.work$time==0)){
      df.work[df.work$time==0,"time"] <- min(df.work[df.work$time != 0, "time"])/100
    }
    
    nWLR <-length(params$mWLR.Delay)
    nAFT <- length(params$AFTdist)
    nRMST <- length(params$RMST.cutoff)
    nmile <- length(params$milestone.cutoff)
    ngen.lin <- length(params$gen.lin$t.low)
    nthresh <- length(params$threshold)
    nABC <- length(params$ABC)
    nMERT <- length(params$MERT$t.low)
    nV0 <- length(params$V0.thresh)
    nCheng <- 1
    npar.group <- length(params$par.group)
    nlogit <- length(params$logit$t.low)
    nlandmark <- length(params$landmark.cutoff)
    npwexp <- length(params$PW.exp)
    npwexplag <- length(params$PW.exp.lag$lag)
    
    #Create matrix for results with 1 row for the statistic and 1 for the p-value
    res <- matrix(NA, ncol = 2,
                  nrow = 42 + nWLR + nRMST + nAFT + 3*nmile + 3*nlandmark + ngen.lin +
                    npwexplag + nthresh + npwexp + nMERT + nV0 + nCheng + npar.group + 
                    2*nlogit)
    colnames(res) <- c("test", "pval")
    res <- as.data.frame(res)
    res$test <- c("G(0,0)", "G(0,1)", "G(1,0)", "G(1,1)", "G(-1,0)", "G(0,0.5)",
                  "G(0,2)", "G(0.5,0.5)", "Zm3", "MaxCombo", "Lee1", "Lee2", "Lee3",
                  "mLee2", "mLee3", "ProjTest", "mScore", "YP", "Cox", 
                  "CoxTD", "CheckPH", "RP.PH", "RP.TD", "AHR", "Aalen", "ABC",
                  "jointTest", "combTest", "CauchyCP", "KS LR", "KS GB", "KS FH", 
                  "KS Cheng", "mZm3", "GB", "TW", "PP", "mPP", "asymLR", "WKM",
                  paste0("GenLin_", 1:length(params$gen.lin$t.low)),
                  paste0("PWExpLag_", 1:length(params$PW.exp.lag$CP)),
                  paste0("Thres_", 1:length(params$threshold)),
                  paste0("PWExp_", 1:length(params$PW.exp)),
                  paste0("MERT_", 1:length(params$MERT$t.low)),
                  paste0("V0_", 1:length(params$V0.thresh)),
                  paste0("MWLRT_", 1:length(params$mWLR.Delay)),
                  paste0("RMST_", 1:length(params$RMST.cutoff)),
                  paste0("AFT_", params$AFTdist),
                  paste0("Mile_", 1:length(params$milestone.cutoff)),
                  paste0("MileCLL_", 1:length(params$milestone.cutoff)),
                  paste0("MileNA_", 1:length(params$milestone.cutoff)),
                  paste0("LLRNA_", 1:length(params$landmark.cutoff)),
                  paste0("QLRNA_", 1:length(params$landmark.cutoff)),
                  paste0("Landmark_", 1:length(params$landmark.cutoff)),
                  paste0("mMaxCombo_", 0.5),
                  paste0("ParGroup_", 1:length(params$par.group)),
                  paste0("Logit_", 1:length(params$logit$t.low)),
                  paste0("mLogit_", 1:length(params$logit$t.low)))
    
    # Define function to calculate p-value of Lee 2 statistics
    Lee2.int <- function(x,w,rho){
      (pnorm((2*x-abs(w)-rho*w)/sqrt(1-rho**2)) - pnorm((-2*x+abs(w)-rho*w)/sqrt(1-rho**2))) * dnorm(w)
    }
    
    #Apply function to calculate self-implemented weighted logrank statistics
    WLR.all <- weightedLR(data = df.work,
                          vars = list(time = 'time', 
                                      event = 'event', 
                                      group = 'group'),
                          params = params)
    
    # FH and combinations-------------------------------------------------
    # MaxCombo and other combinations of the Fleming-Harrington class based on the
    # estimated correlation of the statistics
    
    FlemHar <- WLR.all$FH
    # Calculate 2-sided p-values
    #G(0,0)
    LR.stat <- FlemHar["00","Z"]
    res[res$test=="G(0,0)", "pval"] <- p.LR <- 2*pnorm(-abs(FlemHar["00","Z"]), 
                                                       mean = 0, sd = 1, lower.tail = TRUE)
    #G(0,1)
    res[res$test=="G(0,1)", "pval"] <- 2*pnorm(-abs(FlemHar["01","Z"]), 
                                               mean = 0, sd = 1, lower.tail = TRUE)
    #G(1,0)
    res[res$test=="G(1,0)", "pval"] <- 2*pnorm(-abs(FlemHar["10","Z"]), 
                                               mean = 0, sd = 1, lower.tail = TRUE)
    #G(1,1)
    res[res$test=="G(1,1)", "pval"] <- 2*pnorm(-abs(FlemHar["11","Z"]), 
                                               mean = 0, sd = 1, lower.tail = TRUE)
    
    #Norm covariance matrix to obtain correlation matrix
    cor.mat <- cov.mat <- FlemHar[c("00", "01", "10", "11"),
                                  c("Cov00", "Cov01", "Cov10", "Cov11")]
    for (i in 1:4){
      for (j in 1:4){
        cor.mat[i,j] <- cov.mat[i,j]/sqrt(cov.mat[i,i]*cov.mat[j,j])
      }
    }
    #Zm3
    res[res$test=="Zm3", "pval"] <- (1-mnormt::sadmvn(
      rep(-max(abs(FlemHar[c("00", "01", "10"),"Z"])),3),
      rep(max(abs(FlemHar[c("00", "01", "10"),"Z"])),3),
      mean = rep(0,3), 
      varcov = as.matrix(cor.mat[c("00", "01", "10"),
                                 c("Cov00", "Cov01", "Cov10")])))
    #MaxCombo
    res[res$test == "MaxCombo", "pval"] <- (1-mnormt::sadmvn(
      rep(-max(abs(FlemHar[c("00", "01", "10", "11"),"Z"])),4),
      rep(max(abs(FlemHar[c("00", "01", "10", "11"),"Z"])),4),
      mean = rep(0,4), 
      varcov = as.matrix(cor.mat)))
    #Lee1
    res[res$test=="Lee1", "pval"] <- 2*pnorm(
      -abs(sum(FlemHar[c("01","10"),"Z"]))/sqrt(2+2*cor.mat["01","Cov10"]))
    #Lee2
    res[res$test=="Lee2", "pval"] <- 1-integrate(Lee2.int, 
                                                 lower=-sum(abs(FlemHar[c("01","10"),"Z"])),
                                                 upper = sum(abs(FlemHar[c("01","10"),"Z"])),
                                                 x=(sum(abs(FlemHar[c("01","10"),"Z"])))/2,
                                                 rho=cor.mat["01","Cov10"])$value
    #mLee2
    res[res$test=="mLee2", "pval"] <- 1-integrate(Lee2.int, 
                                                  lower=-sum(abs(FlemHar[c("00","10"),"Z"])),
                                                  upper = sum(abs(FlemHar[c("00","10"),"Z"])),
                                                  x=(sum(abs(FlemHar[c("00","10"),"Z"])))/2,
                                                  rho=cor.mat["00","Cov10"])$value
    #Lee3
    res[res$test == "Lee3", "pval"] <- (1-mnormt::sadmvn(
      rep(-max(abs(FlemHar[c("01","10"),"Z"])),2),
      rep(max(abs(FlemHar[c("01","10"),"Z"])),2),
      mean = rep(0,2), 
      varcov = as.matrix(cor.mat[c("01","10"),c("Cov01", "Cov10")])))
    
    #mLee3
    res[res$test=="mLee3", "pval"] <- (1-mnormt::sadmvn(
      rep(-max(abs(FlemHar[c("00","01"),"Z"])),2),
      rep(max(abs(FlemHar[c("00","01"),"Z"])),2),
      mean = rep(0,2), 
      varcov = cor.mat[c("00","01"),c("Cov00", "Cov01")]))
    
    #G(-1,0)
    res[res$test=="G(-1,0)", "pval"] <- 2*pnorm(-abs(FlemHar["-10","Z"]), 
                                                mean = 0, sd = 1, lower.tail = TRUE)
    #G(0,0.5)
    res[res$test=="G(0,0.5)", "pval"] <- 2*pnorm(-abs(FlemHar["00.5","Z"]), 
                                                 mean = 0, sd = 1, lower.tail = TRUE)
    #G(0,2)
    res[res$test=="G(0,2)", "pval"] <- 2*pnorm(-abs(FlemHar["02","Z"]), 
                                               mean = 0, sd = 1, lower.tail = TRUE)
    #G(0.5,0.5)
    res[res$test=="G(0.5,0.5)", "pval"] <- 2*pnorm(-abs(FlemHar["0.50.5","Z"]), 
                                                   mean = 0, sd = 1, lower.tail = TRUE)
    
    # YP  -----------------------------------------------------
    YP <- WLR.all$AWLRT
    #Calculate correlation matrix
    cor.mat <- cov.mat <- YP[,c("Cov.LR", "Cov.YP1", "Cov.YP2")]
    for (i in 1:3){
      for (j in 1:3){
        cor.mat[i,j] <- cov.mat[i,j]/sqrt(cov.mat[i,i]*cov.mat[j,j])
      }
    }
    # Calculate 2-sided p-values
    res[res$test == "YP", "pval"] <- (1-mnormt::sadmvn(
      rep(-max(abs(YP[c("YP1", "YP2"),"Z"])),2),
      rep(max(abs(YP[c("YP1", "YP2"),"Z"])),2),
      mean = rep(0,2), 
      varcov = cor.mat[c("YP1", "YP2"),c("Cov.YP1", "Cov.YP2")]))
    
    # Remaining weighted LR tests -------------------------------------------------------------
    #Gehan-Breslow
    res[res$test == "GB", "pval"] <- 2*pnorm(-abs(WLR.all$GB[, "Z"]))
    
    #Tarone-Ware
    res[res$test == "TW", "pval"] <- 2*pnorm(-abs(WLR.all$TW[, "Z"]))
    
    #modified Zm3 according to Royston
    mZm3 <- WLR.all$mZm3
    
    cor.mat <- cov.mat <- mZm3[,c("Cov.LR", "Cov.FH(0,1)", "Cov.mFH(1,0)")]
    for (i in 1:3){
      for (j in 1:3){
        cor.mat[i,j] <- cov.mat[i,j]/sqrt(cov.mat[i,i]*cov.mat[j,j])
      }
    }
    res[res$test == "mZm3", "pval"] <- (1-mnormt::sadmvn(
      rep(-max(abs(mZm3[,"Z"])),3),
      rep(max(abs(mZm3[,"Z"])),3),
      mean = rep(0,3), cor.mat))
    
    #modified MaxCombo according to Cheng
    Cheng <- WLR.all$mMaxCombo
    
    cor.mat <- cov.mat <- Cheng[,c("Cov.LR", "Cov.FH(1,0)", "Cov.FH(0,1)", "Cov.Cheng")]
    for (i in 1:4){
      for (j in 1:4){
        cor.mat[i,j] <- cov.mat[i,j]/sqrt(cov.mat[i,i]*cov.mat[j,j])
      }
    }
    res[res$test==paste0("mMaxCombo_", 0.5), "pval"] <- (1-mnormt::sadmvn(
      rep(-max(abs(Cheng[,"Z"])),4),
      rep(max(abs(Cheng[,"Z"])),4),
      mean = rep(0,4), cor.mat))
    #Projection test according to Cheng
    ProjCheng <- Cheng[c("LR", "FH(1,0)", "Cheng"),
                       c("G", "Z", "Cov.LR", "Cov.FH(1,0)", "Cov.Cheng")]
    
    cov.mat <- ProjCheng[,c("Cov.LR", "Cov.FH(1,0)", "Cov.Cheng")]
    
    res[res$test=="ProjTest", "pval"] <- pchisq(
      ProjCheng[,"G"] %*% MASS::ginv(cov.mat) %*% ProjCheng[,"G"],
      df = pracma::Rank(cov.mat),
      lower.tail = FALSE)
    #Modified Score
    MS <- WLR.all$mScore
    MS.stat <- t(MS[,"G"]) %*% solve(MS[,c("Cov.LR", "Cov.NA")]) %*% MS[,"G"]
    res[res$test=="mScore", "pval"] <- pchisq(MS.stat, df=2, lower.tail=FALSE)
    
    #Peto-Peto
    res[res$test == "PP", "pval"] <- 2*pnorm(-abs(WLR.all$PP[, "Z"]))
    
    #modified Peto-Peto
    res[res$test == "mPP", "pval"] <- 2*pnorm(-abs(WLR.all$mPP[, "Z"]))
    
    #asymptotisch LR
    res[res$test == "asymLR", "pval"] <- 2*pnorm(-abs(WLR.all$asymLR[, "Z"]))
    
    # Generalized linear lag model --------------------------------------------
    gen.lin1 <- params$gen.lin$t.low
    gen.lin2 <- params$gen.lin$t.up
    
    res[res$test %in% paste0("GenLin_", 1:length(gen.lin1)), "pval"] <- 
      sapply(WLR.all$gen.lin,
             function(x) 2*pnorm(-abs(x[,"Z"])))
    
    # Threshold lag model -----------------------------------------------------
    res[res$test %in% paste0("Thres_", 1:length(params$threshold)), "pval"] <- 
      sapply(WLR.all$threshold,
             function(x) 2*pnorm(-abs(x[,3])))
    
    # Area between curves (ABC)  -----------------------------------------------------
    #Resampling calculation by Liu
    samples <- split(df.work[, c("time", "event")], f=df.work$group)
    # Naive calculation by Lin and Xu
    Lin.ABC <- LinStatABC.cor(samples$`0`, samples$`1`)
    res[res$test == "ABC", "pval"] <- 2*pnorm(-abs(Lin.ABC), mean = 0, sd = 1)

    # MERT --------------------------------------------------------------------
    #Maximin efficiency robust test
    MERT1 <- params$MERT$t.low
    MERT2 <- params$MERT$t.up
    
    res[res$test %in% paste0("MERT_", 1:length(MERT1)), "pval"] <- 
      sapply(WLR.all$MERT,
             function(x) 2*pnorm(-abs(x[,3])))
    
    # V0 ----------------------------------------------------------------------
    #V0 test/Zucker-Lakatos
    res[res$test %in% paste0("V0_", 1:length(params$V0.thresh)), "pval"] <- 
      sapply(WLR.all$V0,
             function(x) 2*pnorm(-abs(x["LR","Z.V0"])))
    
    
    # Partially grouped LR ----------------------------------------------------
    res[res$test %in% paste0("ParGroup_", 1:length(params$par.group)), "pval"] <- 
      sapply(WLR.all$par.group,
             function(x) 2*pnorm(-abs(x[,"Z"])))
    
    # Logit and modified logit model ------------------------------------------
    Logit1 <- params$logit$t.low
    Logit2 <- params$logit$t.up
    
    res[res$test %in% paste0("Logit_", 1:length(Logit1)), "pval"] <- 
      sapply(WLR.all$logit,
             function(x) 2*pnorm(-abs(x[,"Z"])))
    
    res[res$test %in% paste0("mLogit_", 1:length(Logit1)), "pval"] <- 
      sapply(WLR.all$mlogit,
             function(x) 2*pnorm(-abs(x[,"Z"])))
    
    # Kolmogorov-Smirnov ------------------------------------------------------
    Sun.stat.LR <- max(abs(WLR.all$KS.LR[,"G"]))
    KS.stat.LR <- Sun.stat.LR/sqrt(max(WLR.all$KS.LR[,"Var"]))
    
    Sun.stat.GB <- max(abs(WLR.all$KS.GB[,"G"]))
    KS.stat.GB <- Sun.stat.GB/sqrt(max(WLR.all$KS.GB[,"Var"]))
    
    Sun.stat.FH <- max(abs(WLR.all$KS.FH[,"G"]))
    KS.stat.FH <- Sun.stat.FH/sqrt(max(WLR.all$KS.FH[,"Var"]))
    
    Sun.stat.Cheng <- max(abs(WLR.all$KS.Cheng[,"G"]))
    KS.stat.Cheng <- Sun.stat.Cheng/sqrt(max(WLR.all$KS.Cheng[,"Var"]))
    #Determine p-value by asymptotic distribution, i.e. sup|B(t)| where B is a 
    #standard Brownian motion on [0,1] (Fleming et al. Supremum Versions of the 
    #Log-Rank and Generalized Wilcoxon Statistics)
    BM.pval <- function(x){
      k <- 0:100
      res <- 1 - 4/pi * sum((-1)^k/(2*k+1) * exp(-pi^2*(2*k+1)^2/(8*x^2)))
      return(res)
    }
    
    res[res$test=="KS LR", "pval"] <- BM.pval(KS.stat.LR)
    res[res$test=="KS GB", "pval"] <- BM.pval(KS.stat.GB)
    res[res$test=="KS FH", "pval"] <- BM.pval(KS.stat.FH)
    res[res$test=="KS Cheng", "pval"] <- BM.pval(KS.stat.Cheng)
    
    # AHR ---------------------------------------------------------------------
    #Weighted Cox regression based on Average HR
    Cox.AHR <- tryCatch({coxphw::coxphw(Surv(time, event) ~ group, 
                                        data = df.work, template = "AHR")},
                        error = function(e) NA)
    res[res$test=="AHR", "pval"] <- ifelse(any(class(Cox.AHR)=="coxphw"), 
                                           Cox.AHR$prob, NA)
    
    # Weighted Kaplan-Meier ---------------------------------------------------
    WKM <- nphsim::wkm.Stat(surv=df.work$time, 
                            cnsr=1-df.work$event, 
                            trt=factor(df.work$group, 
                                       levels = 0:1, 
                                       labels=(c("control", "experimental"))))
    res[res$test=="WKM", "pval"] <- 2*pnorm(-abs(WKM$z))
    
    # CauchyCP ----------------------------------------------------------------
    CCP <- CauchyCP::CauchyCP(time = df.work$time, status=df.work$event, x=df.work$group, 
                              covar = rep(1, length(df.work$time)),
                              cutpoints = c(0, quantile(df.work$time[df.work$event== 1])[2:4]))
    res[res$test == "CauchyCP", "pval"] <- CCP$pval
    
    # Modestly weighted LR test -----------------------------------------------
    library(tidyverse)
    df.work %>%
      mutate(event = as.logical(event),
             group = ifelse(group==0, "control", "experimental")) %>%
      select(time, event, group) %>%
      modestWLRT::get_risk_table() ->
      mWLR.risk.table
    modWLRT <- lapply(params$mWLR.Delay,
                      function(x) modestWLRT::add_weights(mWLR.risk.table,  method = "fixed_c",
                                                          delay = x, plot_weights = FALSE))
    res[res$test %in% paste0("MWLRT_", 1:length(params$mWLR.Delay)), "pval"] <- 
      sapply(modWLRT,
             function(x) 2 * pnorm(-abs(modestWLRT::get_zs(x)), mean = 0, sd = 1))
    
    # RMST --------------------------------------------------------------------
    RMST <- lapply(params$RMST.cutoff,
                   function(x) survRM2::rmst2(time = df.work$time,
                                              status = df.work$event,
                                              arm = df.work$group,
                                              tau = x))
    
    res[res$test %in% paste0("RMST_", 1:length(params$RMST.cutoff)), "pval"] <- 
      sapply(RMST, function(x) x$unadjusted.result[1, "p"])
    
    # Cox ---------------------------------------------------------------------
    # Cox model
    Cox <- survival::coxph(Surv(time, event) ~ group, data = df.work)
    res[res$test == "Cox", "pval"] <- p.Cox <- 
      pchisq(Cox$wald.test, df = 1, lower.tail = FALSE) 
    
    #Cox model with time interaction
    CoxTD <- survival::coxph(Surv(time, event) ~ group + tt(group), data = df.work,
                             tt = function(x, t, ...) {
                               mtrx <- model.matrix(~x)[,-1]
                               mtrx*log(t+1)})
    res[res$test=="CoxTD", "pval"] <- p.CoxTD <- 
      pchisq(-2*(CoxTD$loglik[1]-CoxTD$loglik[2]), df = 2, lower.tail = FALSE) 
    
    #Checking PH method (naive)
    alpha_GT <- 0.05
    if (cox.zph(Cox)$table[1, "p"] < alpha_GT){
      res[res$test=="CheckPH", "pval"] <- p.naive <- p.CoxTD
    } else {
      res[res$test=="CheckPH", "pval"] <- p.naive <-  p.Cox
    }
    
    # Royston-Parmar ----------------------------------------------------------
    RP <- tryCatch({flexsurv::flexsurvspline(Surv(time, event) ~ 1, scale = "hazard",
                                             data = df.work, k = 4)},
                   error = function(e) NA)
    # PH model
    RP.PH <- tryCatch({flexsurv::flexsurvspline(Surv(time, event) ~ group, scale = "hazard",
                                                data = df.work, k = 4)},
                      error = function(e) NA)
    res[res$test=="RP.PH", "pval"] <- 
      ifelse(class(RP) == "flexsurvreg" & class(RP.PH) == "flexsurvreg",
             pchisq(-2*(RP$loglik-RP.PH$loglik), df = 1, lower.tail = FALSE), 
             NA)
    
    # TD model
    RP.TD <- tryCatch({flexsurv::flexsurvspline(Surv(time, event) ~ group + gamma1(group), 
                                                scale = "hazard", data = df.work, k = 4)},
                      error = function(e) NA)
    res[res$test=="RP.TD", "pval"] <- 
      ifelse(class(RP) == "flexsurvreg" & class(RP.TD) == "flexsurvreg",
             pchisq(-2*(RP$loglik-RP.TD$loglik), df = 2, lower.tail = FALSE),
             NA)
    
    # Joint and combined test -------------------------------------------------
    #Joint test
    #res[res$test == "joint.test", "Stat"] <- 
    JT <- Cox$score + survival::cox.zph(Cox)$table[1]
    res[res$test == "jointTest", "pval"] <- 1-pchisq(JT, df = 2)
    
    # Combined test and max RMST permutation test
    # Define search grid for maximal RMST statistic
    t.lower <- as.numeric(quantile(df.work$time, 0.3))
    t.upper <- min(max(df.work$time[df.work$group==1]),
                   max(df.work$time[df.work$group==0]))
    
    Cmax <- 0
    for (i in seq(t.lower, t.upper, by = (t.upper - t.lower)/10)){
      hh <- survRM2::rmst2(time = df.work$time, status = df.work$event,
                           arm = df.work$group, tau = i)
      hh.stat <- (hh$RMST.arm1$rmst["Est."] - hh$RMST.arm0$rmst["Est."])^2 / 
        (hh$RMST.arm1$rmst.var + hh$RMST.arm0$rmst.var)
      if (hh.stat > Cmax){
        Cmax <- hh.stat
      }
    }
    
    p.max <- 1 - pchisq(Cmax, df = 1)
    if (p.max > 0.85) {
      p.perm <- 0.9963
    } else {
      p.perm <- 1.762 * p.max^(0.885) - 0.802 * p.max^(2.547)
    }
    p.min <- min(p.perm, p.Cox)
    res[res$test == "combTest", "pval"] <- pbeta(p.min, shape1 = 1, shape2 = 1.5)
    
    # AFT ---------------------------------------------------------------------
    AFT <- lapply(params$AFTdist,
                  function(x) survival::survreg(Surv(time, event) ~ group,
                                                dist=x,
                                                data = df.work))
    
    res[res$test %in% paste0("AFT_", params$AFTdist), "pval"] <- 
      sapply(AFT, function(x) 
        2*pnorm(-abs(x$coefficients["group"] / sqrt(diag(x$var))["group"]),
                mean = 0, sd = 1))
    
    # Piecewise-exponential ---------------------------------------------------
    # PW exp model
    PW <- lapply(params$PW.exp, function(x) 
      if (x>min(df.work$time)){
        summary(eha::pchreg(Surv(time, event) ~ group, 
                            data = df.work, cuts = c(0,x)))
      } else {NA})
    res[res$test %in% paste0("PWExp_", 1:length(params$PW.exp)), "pval"] <- 
      sapply(PW, function(x) 
        if (class(x) == "logical"){NA} else{x$coefficients[,"Wald p"]})
    
    #Data for lagged PW exp model
    lag <- params$PW.exp.lag$lag
    CP <- params$PW.exp.lag$CP - params$PW.exp.lag$lag
    PW.lag.p <- c()
    for (i in 1:length(CP)){
      df.temp <- df.work[df.work$time >= lag[i],]
      df.temp$time <- df.temp$time - lag[i]
      
      PW.lag.p[i] <- tryCatch({
        summary(eha::pchreg(Surv(time, event) ~ group, 
                            data = df.temp, cuts = c(0,CP[i])))$coefficients[,"Wald p"]},
        error = function(e) NA)
    }
    
    res[res$test %in% paste0("PWExpLag_", 1:length(params$PW.exp.lag$CP)), "pval"] <- 
      PW.lag.p
    
    # Aalen -------------------------------------------------------------------
    Aalen <- timereg::aalen(Surv(time, event) ~ const(group), data = df.work)
    res[res$test== "Aalen", "pval"] <- timereg::coef.aalen(Aalen)[,"P-val"]
    
    # Milestone ---------------------------------------------------------------
    #Based on survival function
    fit <- survival::survfit(Surv(time, event) ~ group, data = df.work, conf.type = "none")
    info <- summary(fit, time = params$milestone.cutoff, extend = TRUE)
    df.info0 <- data.frame(t = info$time[info$strata== "group=0"], 
                           p0 = info$surv[info$strata== "group=0"], 
                           se0 = info$std.err[info$strata== "group=0"],
                           std0 = info$std.err[info$strata== "group=0"]^2/
                             info$surv[info$strata== "group=0"]^2)
    df.info1 <- data.frame(t = info$time[info$strata== "group=1"], 
                           p1 = info$surv[info$strata== "group=1"], 
                           se1 = info$std.err[info$strata== "group=1"],
                           std1 = info$std.err[info$strata== "group=1"]^2/
                             info$surv[info$strata== "group=1"]^2)
    
    dplyr::left_join(
      dplyr::left_join(
        data.frame(t = params$milestone.cutoff), df.info0, by = "t"),
      df.info1, by = "t") ->
      df.info
    
    z_stats <- data.frame(
      t = df.info$t,
      diff = df.info$p1 - df.info$p0,
      se = sqrt(df.info$se0^2+df.info$se1^2),
      diff.cloglog = log(-log(df.info$p1)) - log(-log(df.info$p0)),
      se.cloglog = sqrt((df.info$std0/(log(df.info$p0)^2)) + 
                          (df.info$std1/(log(df.info$p1)^2)))
    )
    z_stats$z <- z_stats$diff/z_stats$se
    z_stats$z.cloglog <- z_stats$diff.cloglog/z_stats$se.cloglog
    
    #Naive
    res[res$test %in% paste0("Mile_", 1:length(params$milestone.cutoff)), "pval"] <- 
      2*pnorm(-abs(z_stats$z), mean = 0, sd = 1)
    
    #Cloglog
    res[res$test %in% paste0("MileCLL_", 1:length(params$milestone.cutoff)), "pval"] <- 
      2*pnorm(-abs(z_stats$z.cloglog), mean = 0, sd = 1)
    
    #Based on cumulative hazard
    cumhaz <- data.frame(time = fit$time, cumhaz = fit$cumhaz, std.err = fit$std.err)
    cumhaz0 <- cumhaz[1:fit$strata[1],]
    ind0 <- sapply(params$milestone.cutoff, function(x) max(which(cumhaz0$time<=x)))
    cumhaz1 <- cumhaz[(1+fit$strata[1]):(fit$strata[1]+fit$strata[2]),]
    ind1 <- sapply(params$milestone.cutoff, function(x) max(which(cumhaz1$time<=x)))
    
    NA.stat <- cumhaz0[ind0, "cumhaz"] - cumhaz1[ind1, "cumhaz"]
    NA.var <- sqrt(cumhaz0[ind0, "std.err"]^2 + cumhaz1[ind1, "std.err"]^2)
    
    NA.Z <- NA.stat/NA.var
    res[res$test %in% paste0("MileNA_", 1:length(params$milestone.cutoff)), "pval"] <- 
      2*pnorm(-abs(NA.Z), mean = 0, sd = 1)
    
    #Linear and quadratic combination of threshold LR and cumulative hazard test
    #at the landmark cutoff points
    for (LM.i in 1:length(params$landmark.cutoff)){
      #If CO smaller than first event time in one of the arms NA cannot be estimated
      #at this timepoint and a LR test is performed
      if (params$landmark.cutoff[LM.i]< max(min(cumhaz0$time), min(cumhaz1$time))){
        #Linear combination of cumulative hazard test at fixed timepoint and LR
        res[res$test %in% paste0("LLRNA_", LM.i), "pval"] <- p.LR
        
        #Quadratic combination of cumulative hazard test at fixed timepoint and LR
        res[res$test %in% paste0("QLRNA_", LM.i), "pval"] <- p.LR
      } else {
        #Recalculate for landmark timepoints
        cumhaz <- data.frame(time = fit$time, cumhaz = fit$cumhaz, std.err = fit$std.err)
        cumhaz0 <- cumhaz[1:fit$strata[1],]
        ind0 <- max(which(cumhaz0$time<=params$landmark.cutoff[LM.i]))
        cumhaz1 <- cumhaz[(1+fit$strata[1]):(fit$strata[1]+fit$strata[2]),]
        ind1 <- max(which(cumhaz1$time<=params$landmark.cutoff[LM.i]))
        
        LNA.stat <- cumhaz0[ind0, "cumhaz"] - cumhaz1[ind1, "cumhaz"]
        LNA.var <- sqrt(cumhaz0[ind0, "std.err"]^2 + cumhaz1[ind1, "std.err"]^2)
        LNA.Z <- LNA.stat/LNA.var
        
        #Combine with Z statistic for threshold LR
        LR.thres <- WLR.all$threshold[[LM.i]][,"Z"]
        
        #Linear combination of cumulative hazard test at fixed timepoint and LR
        res[res$test %in% paste0("LLRNA_", LM.i), "pval"] <- 
          2*pnorm(-abs((LNA.Z + LR.thres)/sqrt(2)), mean = 0, sd = 1)
        
        #Quadratic combination of cumulative hazard test at fixed timepoint and LR
        res[res$test %in% paste0("QLRNA_", LM.i), "pval"] <- 
          pchisq(LNA.Z^2 + LR.thres^2, df=2, lower.tail = FALSE)
      }   
    }
    # Landmark ----------------------------------------------------------------
    LM.df <- lapply(params$landmark.cutoff, function(x){
      df.work %>% filter(time >= x) %>% mutate(time = time - x)
    })
    res[res$test %in% paste0("Landmark_", 1:length(params$landmark.cutoff)), "pval"] <- 
      sapply(LM.df, function(x){
        tryCatch(summary(survival::coxph(Surv(time, event)~group,data = x))$coef[,"Pr(>|z|)"],
                 error = function(e) NA)
      })
    return(res)
  }
}


# Exemplary analysis ------------------------------------------------------
df <- data.frame(
  time = c(rexp(150, log(2)/10), rexp(150, log(2)/20)),
  event = rbinom(300, 1, 0.8),
  group = rep(0:1, each=150)
)

MethodsSurvival(data = df, 
                vars = list(time="time", 
                            event = "event", 
                            group = "group"),
                params = list(AFTdist = c("weibull",
                                          "exponential", 
                                          "gaussian", 
                                          "logistic",
                                          "lognormal", 
                                          "loglogistic"),
                              FH = list(rho = c(0,0,1,1,-1,0,0,0.5),
                                        gamma = c(0,1,0,1,0,0.5,2,0.5)),
                              mWLR.Delay = 5,
                              RMST.cutoff = 5, 
                              milestone.cutoff = 5,
                              landmark.cutoff = 5,
                              gen.lin = list(t.low = 5, t.up = 10),
                              PW.exp.lag = list(lag = 5, CP = 10),
                              threshold=5, 
                              ABC=5,
                              PW.exp=5,
                              V0=5, 
                              MERT = list(t.low=5, t.up=10),
                              Cheng = 0.5,
                              par.group=5,
                              logit = list(w.low = 0.1, t.low=5, t.up=10)))
