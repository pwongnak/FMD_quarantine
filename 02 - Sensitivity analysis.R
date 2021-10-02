#############################################################
# Modelling FMD risk through quarantine center
# SEIRV model: Sensitivity analysis
# P. WONGNAK; Last updated; 11 November 2020
# Loading packages ##########################################
library(mc2d)
library(deSolve)
library(ggpubr)
library(ggplot2)
library(epiR)
# Customized function #######################################
round_preserve_sum <- function(x, digits = 0) {
  up <- 10^digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x - y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}
# Parameters ################################################
Duration <- 3
REP <- 5000
Imm <- runif(REP, min = 0, max = 1)
BETA <- runif(REP, min = 0.005, max = 0.015)
R <- runif(REP, min = 0, max = 1)
DELTA <- runif(REP, min = 0, max = 1)
PI <- runif(REP, min = 0, max = 1)
PREV <- runif(REP, min = 0.01, max = 1)
LATENT <- runif(REP, min = 1, max = 5)
SUBCLIN <- runif(REP, min = 1, max = 5)
RECOV <- runif(REP, min = 1, max = 5)
TIMM <- runif(REP, min = 7, max = 28)
#############################################################
ALL_list<- list()
for(I in 1:length(Imm)){
  duration <- Duration
  imm <- Imm[I]
  
  beta <- matrix(NA, ncol = duration, nrow = duration)
  beta[,] <- BETA[I]
  #To make different beta between batches
  for(i in 1:duration){
    for(j in 1:duration){
      if (i != j){
        beta[i,j] <- R[I]*beta[i,j]
      }
    }
  }
  delta <- DELTA[I]
  
  latent <- LATENT[I]
  subclin <- SUBCLIN[I]
  recov <- RECOV[I]
  
  incubation <- latent + subclin
  infectiousP <- subclin + recov
  
  alpha <- 1/TIMM[I]
  sigma <- 1/latent
  pi <- PI[I]
  gammaAR <- 1/infectiousP
  gammaAS <- 1/subclin
  gammaSR <- 1/recov
  
  prev <- PREV[I]
  
  
  #model---------------------------------------------------
  ja.multistage.model <- function (t, x, ...) { 
    s <- x[sindex]
    e <- x[eindex]
    ia <- x[iaindex]
    is <- x[isindex]
    r <- x[rindex]
    v <- x[vindex]
    
    N <- sum(s+e+ia+is+r+v)
    
    dsdt <- ((-s)*(beta%*%(delta * ia + is)))  - (alpha*s)
    dedt <- ((s)*(beta%*%(delta * ia + is))) - (sigma*e)
    diadt <- (sigma*e) - (gammaAR*ia*(1-pi)) - (gammaAS*ia*pi)
    disdt <- (gammaAS*ia*pi) - (gammaSR*is)
    drdt <- (gammaAR*ia*(1-pi)) + (gammaSR*is)
    dvdt <- (alpha*s)
    list(c(dsdt, dedt, diadt, disdt, drdt, dvdt)) 
  }
  
  sindex <- 1:duration
  eindex <- (duration+1):(2*duration)
  iaindex <- (2*duration+1):(3*duration)
  isindex <- (3*duration+1):(4*duration)
  rindex <-  (4*duration+1):(5*duration)
  vindex <- (5*duration+1):(6*duration)
  ########################################################################
  rep <- 1
  marketList <- NULL
  outputALL <- NULL
  for(r in 1:rep){
    nbatch <- round(rtriang(1, min = 1, mode = 3, max = 9))
    ncattle <- NULL
    for(b in 1:nbatch){
      ncattle[b] <- round(rpois(1, lambda = 30))
    }
    cattleinit <- sum(ncattle)
    cattleinf <- rbinom(1, size = cattleinit, prob = prev)
    cattleE <- rbinom(1, size = cattleinf, prob = 1/3)
    cattleIS <- 0
    cattleIA <- cattleinf- cattleE - cattleIS
    cattleR <- rbinom(1, size = cattleinit, prob = imm)
    cattleS <- cattleinit-cattleinf-cattleR
    if(cattleS < 0) {cattleS <- 0}
    
    yinit <- c(
      S=c(cattleS, rep(0,(duration-1))),
      E=c(cattleE, rep(0, duration-1)),
      IA=c(cattleIA, rep(0,(duration-1))),
      IS=c(cattleIS, rep(0,(duration-1))),
      R=c(cattleR, rep(0,(duration-1))),
      V=c(rep(0,(duration)))
    )
    
    output <- matrix(NA, nrow = 1, ncol = length(yinit)+3)
    output <- data.frame(output)
    names(output) <- c("rep", "time", "week", names(yinit))
    output[1,] <- c(r, 0, 0, yinit)
    
    runtime <- 52*1
    
    market <- matrix(NA, nrow = runtime, ncol = 8)
    for(i in 1:runtime){
      sol <- ode(y=yinit, times=seq(0,7,by=0.1), func=ja.multistage.model, method = "euler", parms = NULL)
      weekoutput <- sol[-1,-1]
      week <- i
      time <- (((i-1)*7)+1):(((i-1)*7)+7)
      sol <- data.frame(sol)
      sol <- sol[sol$time%%1==0,-1]
      weekoutput <- data.frame(cbind("rep" = r, "time" = time, "week" = rep(week, 7), sol[-1,]))
      output <- rbind(output, weekoutput)
      
      weekstatus <- as.numeric(sol[8,])
      market[i,1:7] <- c(i, round_preserve_sum(weekstatus[seq(duration,6*duration,duration)]))
      market[i,8] <- r
      send <- weekstatus[sindex]
      eend <- weekstatus[eindex]
      iaend <- weekstatus[iaindex]
      isend <- weekstatus[isindex]
      rend <- weekstatus[rindex]
      vend <- weekstatus[vindex]
      
      
      nbatch <- round(rtriang(1, min = 1, mode = 3, max = 9))
      ncattle <- NULL
      for(b in 1:nbatch){
        ncattle[b] <- round(rpois(1, lambda = 30))
      }
      
      cattleinit <- sum(ncattle)
      cattleinf <- rbinom(1, size = cattleinit, prob = prev)
      cattleE <- rbinom(1, size = cattleinf, prob = 1/3)
      cattleIS <- 0
      cattleIA <- cattleinf- cattleE - cattleIS
      cattleR <- rbinom(1, size = cattleinit, prob = imm)
      cattleS <- cattleinit-cattleinf-cattleR
      if(cattleS < 0) {cattleS <- 0}
      
      simp <- cattleS
      eimp <- cattleE
      iaimp <- cattleIA
      isimp <- 0
      rimp <- cattleR
      vimp <- 0
      
      if(duration == 1){
        yinit <-  c(
          S=c(simp),
          E=c(eimp),
          IA=c(iaimp),
          IS=c(isimp),
          R=c(rimp),
          V=c(vimp)
        ) }else {
          yinit <- c(
            S=c(simp, send[1:(length(sindex)-1)]),
            E=c(eimp, eend[1:(length(eindex)-1)]),
            IA=c(iaimp, iaend[1:(length(iaindex)-1)]),
            IS=c(isimp, isend[1:(length(iaindex)-1)]),
            R=c(rimp, rend[1:(length(rindex)-1)]),
            V=c(vimp, vend[1:(length(vindex)-1)])
          )
        }
    }
    market <- data.frame(market)
    outputALL <- rbind(outputALL, output)
    colnames(market) <- c("week", "S", "E", "IA", "IS", "R", "V", "Rep")
    marketList <- rbind(marketList, market)
    conditionList <- c(Imm[I], BETA[I], R[I], DELTA[I], PI[I], PREV[I], gammaAR, gammaAS, gammaSR, alpha)
    ALL <- list(conditionList, marketList)
    ALL_list[[I]] <- ALL
  }
}
########################################################################
All_summary <- NULL
for(i in 1:length(ALL_list)){
  
  weeksummary <- ALL_list[[i]][[2]][-1:-5,]
  risk <- sum(weeksummary$IA)/sum((weeksummary$S + weeksummary$E + weeksummary$IA + weeksummary$R + weeksummary$V))
  prot <- sum(weeksummary$R + weeksummary$V)/sum((weeksummary$S + weeksummary$E + weeksummary$IA + weeksummary$R + weeksummary$V))
  erad  <- sum(weeksummary$IS)/sum((weeksummary$S + weeksummary$E + weeksummary$IA + weeksummary$R + weeksummary$V + weeksummary$IS))
  
  summary <- data.frame("Imm" = ALL_list[[i]][[1]][1],
                        "Beta" = ALL_list[[i]][[1]][2],
                        "R" = ALL_list[[i]][[1]][3],
                        "Delta" = ALL_list[[i]][[1]][4],
                        "Pi" = ALL_list[[i]][[1]][5],
                        "Prev" = ALL_list[[i]][[1]][6],
                        "Gamma1" = ALL_list[[i]][[1]][7],
                        "Gamma2" = ALL_list[[i]][[1]][8],
                        "Gamma3" = ALL_list[[i]][[1]][9],
                        "Alpha" = ALL_list[[i]][[1]][10],
                        "risk" = risk,
                        "prot" = prot,
                        "erad" = erad)
  All_summary <- rbind(All_summary, summary)
  
}

All_summary$Imm_r <- rank(All_summary$Imm)
All_summary$Beta_r <- rank(All_summary$Beta)
All_summary$R_r <- rank(All_summary$R)
All_summary$Delta_r <- rank(All_summary$Delta)
All_summary$Pi_r <- rank(All_summary$Pi)
All_summary$Prev_r <- rank(All_summary$Prev)
All_summary$Gamma1_r <- rank(All_summary$Gamma1)
All_summary$Gamma2_r <- rank(All_summary$Gamma2)
All_summary$Gamma3_r <- rank(All_summary$Gamma3)
All_summary$Alpha_r <- rank(All_summary$Alpha)

All_summary$risk_r <- rank(All_summary$risk)
All_summary$prot_r <- rank(All_summary$prot)
All_summary$erad_r <- rank(All_summary$erad)

par(mfrow = c(3,3))

plot(All_summary$Imm_r, All_summary$risk_r)
plot(All_summary$Beta_r, All_summary$risk_r)
plot(All_summary$R_r, All_summary$risk_r)
plot(All_summary$Delta_r, All_summary$risk_r)
plot(All_summary$Pi_r, All_summary$risk_r)
plot(All_summary$Prev_r, All_summary$risk_r)
plot(All_summary$Gamma1_r, All_summary$risk_r)
plot(All_summary$Gamma2_r, All_summary$risk_r)
plot(All_summary$Gamma3_r, All_summary$risk_r)

plot(All_summary$Imm_r, All_summary$prot_r)
plot(All_summary$Beta_r, All_summary$prot_r)
plot(All_summary$R_r, All_summary$prot_r)
plot(All_summary$Delta_r, All_summary$prot_r)
plot(All_summary$Pi_r, All_summary$prot_r)
plot(All_summary$Prev_r, All_summary$prot_r)
plot(All_summary$Gamma1_r, All_summary$prot_r)
plot(All_summary$Gamma2_r, All_summary$prot_r)
plot(All_summary$Gamma3_r, All_summary$prot_r)

plot(All_summary$Imm_r, All_summary$erad_r)
plot(All_summary$Beta_r, All_summary$erad_r)
plot(All_summary$R_r, All_summary$erad_r)
plot(All_summary$Delta_r, All_summary$erad_r)
plot(All_summary$Pi_r, All_summary$erad_r)
plot(All_summary$Prev_r, All_summary$erad_r)
plot(All_summary$Gamma1_r, All_summary$erad_r)
plot(All_summary$Gamma2_r, All_summary$erad_r)
plot(All_summary$Gamma3_r, All_summary$erad_r)
################################################
col1 <- c("Imm_r"   , "Beta_r"  , "R_r"   ,   "Delta_r" , "Pi_r"  ,   "Prev_r",
          "Gamma1_r" ,"Gamma2_r" ,"Gamma3_r", "Alpha", "risk_r")
col2 <- c("Imm_r"   , "Beta_r"  , "R_r"   ,   "Delta_r" , "Pi_r"  ,   "Prev_r",
          "Gamma1_r" ,"Gamma2_r" ,"Gamma3_r", "Alpha", "prot_r")
col3 <- c("Imm_r"   , "Beta_r"  , "R_r"   ,   "Delta_r" , "Pi_r"  ,   "Prev_r",
          "Gamma1_r" ,"Gamma2_r" ,"Gamma3_r", "Alpha", "erad_r")

testData1 <- All_summary[,col1]
testData2 <- All_summary[,col2]
testData3 <- All_summary[,col3]

#PRCC#############################################
P1<-cbind(epi.prcc(testData1), "param" = col1[-length(col1)])
P2<-cbind(epi.prcc(testData2), "param" = col1[-length(col2)])
P3<-cbind(epi.prcc(testData3), "param" = col1[-length(col3)])

P1$col <- P1$gamma > 0
P1$sig <- ifelse(P1$p.value < 0.05, "*", "ns")
P1$hjust <- ifelse(P1$gamma > 0, 1.4, -0.4)

P2$col <- P2$gamma > 0
P2$sig <- ifelse(P2$p.value < 0.05, "*", "ns")
P2$hjust <- ifelse(P2$gamma > 0, 1.4, -0.4)

P3$col <- P3$gamma > 0
P3$sig <- ifelse(P3$p.value < 0.05, "*", "ns")
P3$hjust <- ifelse(P3$gamma > 0, 1.4, -0.4)

G1 <- ggplot(data = P1, aes(x = param, y = gamma, fill = col, label = sig, hjust = hjust)) + geom_bar(stat = "identity") +
  theme_bw() +
  coord_flip() +
  ylim(-0.65,0.65) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = c(expression(alpha), expression(beta), expression(delta), expression(gamma[1]), expression(gamma[2]),
                              expression(gamma[3]), expression("P"["imm"]), expression(pi),
                              expression("P"["inf"]), expression("r"))) +
  xlab("Parameter") +
  ylab("Partial rank correlation coefficients") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#ea2c62", "#7579e7")) +
  ggtitle(expression("Sensitivity for P"["exp"])) +
  geom_text(aes(y = 0))
G1

G2 <- ggplot(data = P2, aes(x = param, y = gamma, fill = col, label = sig, hjust = hjust)) + geom_bar(stat = "identity") +
  theme_bw() +
  coord_flip() +
  ylim(-0.65,0.65) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = c(expression(alpha), expression(beta), expression(delta), expression(gamma[1]), expression(gamma[2]),
                              expression(gamma[3]), expression("P"["imm"]), expression(pi),
                              expression("P"["inf"]), expression("r"))) +
  xlab("Parameter") +
  ylab("Partial rank correlation coefficients") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#ea2c62", "#7579e7")) +
  ggtitle(expression("Sensitivity for P"["pro"])) +
  geom_text(aes(y = 0))
G2

G3 <- ggplot(data = P3, aes(x = param, y = gamma, fill = col, label = sig, hjust = hjust)) + geom_bar(stat = "identity") +
  theme_bw() +
  coord_flip() +
  ylim(-0.65,0.65) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = c(expression(alpha), expression(beta), expression(delta), expression(gamma[1]), expression(gamma[2]),
                              expression(gamma[3]), expression("P"["imm"]), expression(pi),
                              expression("P"["inf"]), expression("r"))) +
  xlab("Parameter") +
  ylab("Partial rank correlation coefficients") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#ea2c62", "#7579e7")) +
  ggtitle(expression("Sensitivity for P"["rem"])) +
  geom_text(aes(y = 0))
G3

ggarrange(G1, G2, G3, ncol = 3, nrow = 1, labels = "AUTO")









