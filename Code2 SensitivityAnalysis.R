########################################################################################################
# P. Wongnak
# Scripts for the sensivity analysis for the article "A stochastic modeling study of quarantine strategies against 
# foot-and-mouth disease risks through cattle trades across the Thailand-Myanmar border"
########################################################################################################
mainDir <- "C:/" ### User define
setwd(mainDir)
########################################################################################################
library(odin)
library(dde)
library(mc2d)
library(stringr)
library(epiR)
########################################################################################################
#Model
seir <- odin({
  ## model
  update(S[]) <- S[i] - immunized_s[i] - new_infection[i] 
  update(E[]) <- E[i] + new_infection[i] - new_asymptomatic1[i] - new_asymptomatic2[i]
  update(IA1[]) <- IA1[i] + new_asymptomatic1[i] - recover_from_IA1[i] 
  update(IA2[]) <- IA2[i] + new_asymptomatic2[i] - new_symptomatic[i] 
  update(Is[]) <- Is[i] + new_symptomatic[i] - recover_from_Is[i] 
  update(R[]) <- R[i] + recover_from_Is[i] + recover_from_IA1[i] + immunized_s[i]
  ### stochastic process
  ## immunization process
  nWeek <- length(initS)
  p_immunized[] <- interpolate(time, p_immuned, "linear") #Time to immunized
  immunized_s[] <- rbinom(S[i], vaccine_effectiveness * p_immunized[i])
  ## infection process
  betaI[,] <- beta[i,j] * (Is[j] +  delta_asymptomatic * (IA1[j] + IA2[j])) #transmission between batch i and j
  p_new_infection[] <- 1 - exp(-sum(betaI[i,])) #Prob S >> E
  new_infection[] <- rbinom(S[i]-immunized_s[i], p_new_infection[i]) # Number S >> E
  ## infectious process
  p_new_asymptomatic <- 1 - exp(-alpha) # Prob E >>> IA 
  new_asymptomatic[] <- rbinom(E[i], p_new_asymptomatic) # Number E >>> IA = IA1 + IA2
  new_asymptomatic1[] <- rbinom(new_asymptomatic[i], (1-pi)) # Number of IA1 from IA
  new_asymptomatic2[] <- new_asymptomatic[i] - new_asymptomatic1[i] # Number of IA2 from IA
  ## recovery process from asymptomatic
  p_recover_from_IA1 <- 1 - exp(-gamma_IA1_R) #Prob IA1 >>> R 
  recover_from_IA1[] <- rbinom(IA1[i], p_recover_from_IA1) # Number IA1 >>> R
  ## process of symptomatic cattle
  P_new_symptomatic <- 1 - exp(-gamma_IA2_Is) # Prob IA2 >>> Is
  new_symptomatic[] <- rbinom(IA2[i], P_new_symptomatic) # Number IA2 >>> Is
  p_recover_from_Is <- 1 - exp(-gamma_Is_R) #Prob Is >>> R 
  recover_from_Is[] <- rbinom(Is[i], p_recover_from_Is) # Number Is >>> R
  ## Total cattle
  N <- sum(S[]) + sum(E[]) + sum(IA1[]) + sum(IA2[]) + sum(Is[]) + sum(R[])
  ## Initial values
  initial(S[]) <- initS[i]
  initial(E[]) <- initE[i]
  initial(IA1[]) <- initIA1[i]
  initial(IA2[]) <- initIA2[i]
  initial(Is[]) <- initIs[i]
  initial(R[]) <- initR[i]
  
  initS[] <- user()
  initE[] <- user()
  initIA1[] <- user()
  initIA2[] <- user()
  initIs[] <- user()
  initR[] <- user()
  
  time[] <- user()
  p_immuned[,] <- user()
  
  ## Parameters
  delta_asymptomatic <- user(0.99) #Proportional infectious of asymptomatic cattle
  pi <- user(0.68) #Proportion of infectious cattle become symptomatic
  beta[,] <- user() #Transmission coefficient (Density dependent)
  alpha <- 1/latent_period #Rate of latent becoming infectious
  gamma_IA1_R <- 1/infectious_period #Recovery rate from IA1
  gamma_IA2_Is <- 1/(incubation_period - latent_period) #Rate of showing clinical signs
  gamma_Is_R <- 1/(infectious_period - (incubation_period - latent_period)) #Recovery rate from IA1
  
  latent_period <- user(3.6) #latent periods; time from exposed to infectious
  infectious_period <- user(4.4) # infectious period; time stay infectious
  incubation_period <- user(5.9) # incubation period; time from exposed to clinical sign
  vaccine_effectiveness <- user(0.80)
  
  ## Dimensions
  dim(S) <- nWeek
  dim(E) <- nWeek
  dim(IA1) <- nWeek
  dim(IA2) <- nWeek
  dim(Is) <- nWeek
  dim(R) <- nWeek
  dim(immunized_s) <- nWeek
  dim(p_new_infection) <- nWeek
  
  dim(new_infection) <- nWeek
  dim(new_asymptomatic) <- nWeek
  dim(new_asymptomatic1) <- nWeek
  dim(new_asymptomatic2) <- nWeek
  dim(recover_from_IA1) <- nWeek
  dim(new_symptomatic) <- nWeek
  dim(recover_from_Is) <- nWeek
  dim(p_immunized) <- nWeek
  
  
  dim(time) <- user()
  dim(p_immuned) <- c(length(time), nWeek)
  dim(initS) <- user()
  dim(initE) <- user()
  dim(initIA1) <- user()
  dim(initIA2) <- user()
  dim(initIs) <- user()
  dim(initR) <- user()
  dim(beta) <- c(nWeek, nWeek)
  dim(betaI) <- c(nWeek, nWeek)
  
  output(N) <- N
  
})

###########################################################################################################
nWeek <- 3
# running time per week
tmax <- 7
tt <- seq(0, tmax, length.out = tmax+1)
time <- 0:tmax

# Transmission rate
b <- 0.01 #transmission rate (density-dependent)
RList <- seq(0,1,0.1) #reduction factor for between batch transmission

REP <- 5000
IMM <- runif(REP, min = 0, max = 1)
BETA <- runif(REP, min = 0.005, max = 0.015)
R <- runif(REP, min = 0, max = 1)
DELTA_ASYM <- runif(REP, min = 0, max = 1)
DELTA_VAC <- runif(REP, min = 0, max = 1)
PI <- runif(REP, min = 0, max = 1)
PREV <- runif(REP, min = 0.01, max = 1)
LATENT <- runif(REP, min = 1, max = 5)
INFECTIOUS <- runif(REP, min = 1, max = 5)
INCUBATION <- rep(NA, REP)
for(i in 1:REP){INCUBATION[i] <- runif(1, min = LATENT[i], max =  LATENT[i]+INFECTIOUS[i])}
VACCINE_EFFECTIVENESS <- runif(REP, min = 0, max = 1)

# Importing cattle ################################################
Import_cattle <- function(Imm, prev, pi){
  nbatch <- round(rtriang(1, min = 1, mode = 3, max = 9))
  ncattle <- NULL
  for(batch in 1:nbatch){
    ncattle[batch] <- round(rpois(1, lambda = runif(1, 10, 30)))
  }
  cattleinit <- sum(ncattle) # total number of cattle upon importation
  cattleInf <- rbinom(1, cattleinit, prev)
  cattleE <- rbinom(1, size = cattleInf, prob =  0.45) #3.6/(3.6+4.4)
  cattleI <- cattleInf - cattleE
  cattleIA1 <- rbinom(1, size = cattleI, prob = 1-pi)
  cattleIs <- rbinom(1, size = cattleI - cattleIA1, prob = 0.5) #(5.9-3.6)/4.4
  cattleIA2 <- cattleI - cattleIA1 - cattleIs
  cattleR <- rbinom(1, size = (cattleinit - cattleInf), prob = Imm)
  cattleS <- cattleinit - cattleInf - cattleR
  cattleIs <- 0
  output <- data.frame(cattleS, cattleE, cattleIA1, cattleIA2, cattleIs, cattleR)
  return(output)
}
###########################################################################
ALL_list <- list()
for(i in 1:REP){
  print(i)
  # Transmission rate
  b <- BETA[i] #transmission rate (density-dependent)
  r <- R[i] #reduction factor for between batch transmission
  br <- b * r
  
  # reduction of transmission from asymptpmatic cattle
  delta_asymptomatic <- DELTA_ASYM[i]
  # Probability of imported cattle being FMDV infected; P_inf
  prev <- PREV[i]
  # Proportions of immune cattle upon arrival; Imm_1(0)
  Imm <- IMM[i]
  # Proportion of infectious individual becoming symptomatic; pi
  pi <- PI[i]
  #vaccine_effectiveness
  vaccine_effectiveness <- VACCINE_EFFECTIVENESS[i]
  
  latent_period <- LATENT[i]
  incubation_period <- INCUBATION[i]
  infectious_period <- INFECTIOUS[i]
  
  # Vaccine-induced immunity
  immune_v_function <- function(time, t_imm_min = 4, t_imm_max = 11){
    output <- (punif(time+1, t_imm_min, t_imm_max) - punif(time, t_imm_min, t_imm_max)) / (1 - punif(time, t_imm_min, t_imm_max))
    output[which(is.nan(output))] <- 0
    return(output)
  } # Function of vaccine-induced immunization
  p_immuned <- matrix(NA, nrow = length(time), ncol = nWeek)
  for(nw in 1:nWeek){
    p_immuned[,nw] <- immune_v_function(time + (7*(nw-1)))
  }
  
  beta <- matrix(b, nrow = nWeek, ncol = nWeek)
  for(k in 1:nWeek){
    for(j in 1:nWeek){
      if (k != j){
        beta[k,j] <- br
      }
    }
  }
  
  init_zero <- rep(0, (nWeek-1))
  ## initiate the model
  imported_cattle <- Import_cattle(Imm, prev, pi)
  initS <- c(imported_cattle$cattleS, init_zero)
  initE <- c(imported_cattle$cattleE, init_zero)
  initIA1 <- c(imported_cattle$cattleIA1, init_zero)
  initIA2 <- c(imported_cattle$cattleIA2, init_zero)
  initIs <- c(imported_cattle$cattleIs, init_zero)
  initR <- c(imported_cattle$cattleR, init_zero)
  
  maxit <- 52 # maximum iteration
  ToMarket <- NULL
  AllResults <- NULL
  
  for(w in 1:maxit){
    week <- w
    params <- list(initS = initS,
                   initE = initE,
                   initIA1 = initIA1,
                   initIA2 = initIA2,
                   initIs = initIs,
                   initR = initR,
                   p_immuned = p_immuned,
                   time = time,
                   beta = beta,
                   delta_asymptomatic = delta_asymptomatic,
                   latent_period = latent_period,
                   infectious_period = infectious_period,
                   incubation_period = incubation_period,
                   vaccine_effectiveness = vaccine_effectiveness
    )
    #model
    mod <- seir$new(user = params)
    #output
    (y <- data.frame(mod$run(tt)))
    
    if(w == 1){y = y}else{y = y[-1,]}
    y$day <- y$step + (7*(w-1))
    y <- cbind("rep" = i, y)
    AllResults <- rbind(AllResults, y)
    
    flagLastWeek <- which(str_detect(colnames(y), paste0(".", nWeek, ".")))
    ToMarket_sub <- cbind("rep" = i, "week" = week, y[nrow(y),flagLastWeek])
    #ToMarket_sub$N <- sum(ToMarket_sub[,-1])
    ToMarket <- rbind(ToMarket, ToMarket_sub)
    conditionList <- c("b" = b, "r" = r, "delta_asymptomatic" = delta_asymptomatic, 
                       "vaccine_effectiveness" = vaccine_effectiveness, "prev" = prev, "Imm" = Imm, "pi" = pi,  "gamma_IA1_R" = 1/infectious_period,
                       "gamma_IA2_IS" = 1/(incubation_period - latent_period), "gamma_IS_R" = 1/(infectious_period - (incubation_period - latent_period)))
    #changing the class
    S_end <- c(as.matrix(y[nrow(y), which(str_detect(colnames(y), paste0("S.")))]))
    E_end <- c(as.matrix(y[nrow(y), which(str_detect(colnames(y), paste0("E.")))]))
    IA1_end <- c(as.matrix(y[nrow(y), which(str_detect(colnames(y), paste0("IA1.")))]))
    IA2_end <- c(as.matrix(y[nrow(y), which(str_detect(colnames(y), paste0("IA2.")))]))
    Is_end <- c(as.matrix(y[nrow(y), which(str_detect(colnames(y), paste0("Is.")))]))
    R_end <- c(as.matrix(y[nrow(y), which(str_detect(colnames(y), paste0("R.")))]))
    
    #importing the next batch
    imported_cattle <- Import_cattle(Imm, prev, pi)
    
    initS <- c(imported_cattle$cattleS, S_end[-nWeek])
    initE <- c(imported_cattle$cattleE, E_end[-nWeek])
    initIA1 <- c(imported_cattle$cattleIA1, IA1_end[-nWeek])
    initIA2 <- c(imported_cattle$cattleIA2, IA2_end[-nWeek])
    initIs <- c(imported_cattle$cattleIs, Is_end[-nWeek])
    initR <- c(imported_cattle$cattleR, R_end[-nWeek])
    
  }
  ALL <- list(conditionList, ToMarket)
  ALL_list[[i]] <- ALL
}

str(ALL_list)
########################################################################
All_summary <- NULL
i <- 1
for(i in 1:length(ALL_list)){
  
  data <- ALL_list[[i]][[2]][-1:-5,]
  
  flagS <- which(str_detect(colnames(data), "S"))
  flagE <- which(str_detect(colnames(data), "E"))
  flagIA <- which(str_detect(colnames(data), "IA"))
  flagIS <- which(str_detect(colnames(data), "Is"))
  flagR <- which(str_detect(colnames(data), "R"))

  data$S <- data[,flagS]
  data$E <- data[,flagE]
  data$IA <- rowSums(data[,flagIA])
  data$IS <- data[,flagIS]
  data$R <- data[,flagR]
  
  risk <- sum(data$IA)/sum((data$S + data$E + data$IA + data$R))
  prot <- sum(data$R)/sum((data$S + data$E + data$IA + data$R))
  erad  <- sum(data$IS)/sum((data$S + data$E + data$IA + data$R + data$IS))
  
  weekNotSafe <- sum(!(data$E == 0 & data$IA == 0))/nrow(data)
  
  summary <- data.frame("Beta" = ALL_list[[i]][[1]][1],
                        "R" = ALL_list[[i]][[1]][2],
                        "Delta_asymp" = ALL_list[[i]][[1]][3],
                        "vaccine_effectiveness" = ALL_list[[i]][[1]][4],
                        "prev" = ALL_list[[i]][[1]][5],
                        "Imm" = ALL_list[[i]][[1]][6],
                        "pi" = ALL_list[[i]][[1]][7],
                        "gamma_IA1_R" = ALL_list[[i]][[1]][8],
                        "gamma_IA2_IS" = ALL_list[[i]][[1]][9],
                        "gamma_IS_R" = ALL_list[[i]][[1]][10],
                        "risk" = risk,
                        "prot" = prot,
                        "erad" = erad,
                        "weekNotSafe" = weekNotSafe)
  All_summary <- rbind(All_summary, summary)
}

All_summary$Beta_r <- rank(All_summary$Beta)
All_summary$R_r <- rank(All_summary$R)
All_summary$Delta_asymp_r <- rank(All_summary$Delta_asymp)
All_summary$vaccine_effectiveness_r <- rank(All_summary$vaccine_effectiveness)

All_summary$prev_r <- rank(All_summary$prev)
All_summary$Imm_r <- rank(All_summary$Imm)
All_summary$pi_r <- rank(All_summary$pi)
All_summary$gamma_IA1_R_r <- rank(All_summary$gamma_IA1_R)
All_summary$gamma_IA2_IS_r <- rank(All_summary$gamma_IA2_IS)
All_summary$gamma_IS_R_r <- rank(All_summary$gamma_IS_R)
All_summary$risk_r <- rank(All_summary$risk)
All_summary$prot_r <- rank(All_summary$prot)
All_summary$erad_r <- rank(All_summary$erad)
All_summary$weekNotSafe_r <- rank(All_summary$weekNotSafe)

write.csv(All_summary, "Sensitivity_results.csv", row.names = F)
colnames(All_summary)

################################################
All_summary <- read.csv("C:/Users/WONGNAK/Dropbox/Project/JSPS FMD/NEW2_SDE/Sensitivity_results.csv")


col1 <- c("Beta_r", "R_r", "Delta_asymp_r",  "vaccine_effectiveness_r",          
          "prev_r", "Imm_r", "pi_r", "gamma_IA2_IS_r", "gamma_IA1_R_r", "gamma_IS_R_r",  "risk_r")
col2 <- c("Beta_r", "R_r", "Delta_asymp_r",  "vaccine_effectiveness_r",           
          "prev_r", "Imm_r", "pi_r", "gamma_IA2_IS_r", "gamma_IA1_R_r", "gamma_IS_R_r",  "prot_r")
col3 <- c("Beta_r", "R_r", "Delta_asymp_r",  "vaccine_effectiveness_r",          
          "prev_r", "Imm_r", "pi_r", "gamma_IA2_IS_r", "gamma_IA1_R_r", "gamma_IS_R_r",  "erad_r")
col4 <- c("Beta_r", "R_r", "Delta_asymp_r",  "vaccine_effectiveness_r",           
          "prev_r", "Imm_r", "pi_r", "gamma_IA2_IS_r", "gamma_IA1_R_r", "gamma_IS_R_r", "weekNotSafe_r")

testData1 <- All_summary[,col1]
testData2 <- All_summary[,col2]
testData3 <- All_summary[,col3]
testData4 <- All_summary[,col4]
#PRCC#############################################
P1<-cbind(epi.prcc(testData1), "param" = col1[-length(col1)])
P2<-cbind(epi.prcc(testData2), "param" = col1[-length(col2)])
P3<-cbind(epi.prcc(testData3), "param" = col1[-length(col3)])
P4<-cbind(epi.prcc(testData4), "param" = col1[-length(col4)])

P1$col <- P1$gamma > 0
P1$sig <- ifelse(P1$p.value < 0.05, "*", "ns")
P1$hjust <- ifelse(P1$gamma > 0, 1.4, -0.4)

P2$col <- P2$gamma > 0
P2$sig <- ifelse(P2$p.value < 0.05, "*", "ns")
P2$hjust <- ifelse(P2$gamma > 0, 1.4, -0.4)

P3$col <- P3$gamma > 0
P3$sig <- ifelse(P3$p.value < 0.05, "*", "ns")
P3$hjust <- ifelse(P3$gamma > 0, 1.4, -0.4)

P4$col <- P4$gamma > 0
P4$sig <- ifelse(P4$p.value < 0.05, "*", "ns")
P4$hjust <- ifelse(P4$gamma > 0, 1.4, -0.4)


P1$param <- factor(P1$param, levels = c("Beta_r", "R_r", "Delta_asymp_r",  "vaccine_effectiveness_r",          
                                        "prev_r", "Imm_r", "pi_r", "gamma_IA2_IS_r", "gamma_IA1_R_r", "gamma_IS_R_r"))
P2$param <- factor(P2$param, levels = c("Beta_r", "R_r", "Delta_asymp_r",  "vaccine_effectiveness_r",          
                                        "prev_r", "Imm_r", "pi_r", "gamma_IA2_IS_r", "gamma_IA1_R_r", "gamma_IS_R_r"))
P3$param <- factor(P3$param, levels = c("Beta_r", "R_r", "Delta_asymp_r",  "vaccine_effectiveness_r",          
                                        "prev_r", "Imm_r", "pi_r", "gamma_IA2_IS_r", "gamma_IA1_R_r", "gamma_IS_R_r"))
P4$param <- factor(P4$param, levels = c("Beta_r", "R_r", "Delta_asymp_r",  "vaccine_effectiveness_r",          
                                        "prev_r", "Imm_r", "pi_r", "gamma_IA2_IS_r", "gamma_IA1_R_r", "gamma_IS_R_r"))


G1 <- ggplot(data = P1, aes(x = param, y = gamma, fill = col, label = sig, hjust = hjust)) + geom_bar(stat = "identity") +
  theme_bw() +
  coord_flip() +
  ylim(-1,1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = c(expression(beta), expression("r"), expression(delta),expression("e"["v"]),
                              expression("P"["inf"]), expression("P"["R(0)"]), expression(pi),
                                  expression(gamma[1]), expression(gamma[2]), expression(gamma[3]))) +
  xlab("Parameter") +
  ylab("Partial rank correlation coefficients") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#ea2c62", "#7579e7")) +
  ggtitle(expression("Risk of released cattle being infected, " *italic(P["cattle"]))) +
  geom_text(aes(y = 0))
G1

G2 <- ggplot(data = P2, aes(x = param, y = gamma, fill = col, label = sig, hjust = hjust)) + geom_bar(stat = "identity") +
  theme_bw() +
  coord_flip() +
  ylim(-1,1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = c(expression(beta), expression("r"), expression(delta),expression("e"["v"]),
                              expression("P"["inf"]), expression("P"["R(0)"]), expression(pi),
                              expression(gamma[1]), expression(gamma[2]), expression(gamma[3]))) +
  xlab("Parameter") +
  ylab("Partial rank correlation coefficients") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#ea2c62", "#7579e7")) +
  ggtitle(expression("Probability of releasing immunized cattle, " *italic(P[immunized]))) +
  geom_text(aes(y = 0))
G2

G3 <- ggplot(data = P3, aes(x = param, y = gamma, fill = col, label = sig, hjust = hjust)) + geom_bar(stat = "identity") +
  theme_bw() +
  coord_flip() +
  ylim(-1,1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = c(expression(beta), expression("r"), expression(delta),expression("e"["v"]),
                              expression("P"["inf"]), expression("P"["R(0)"]), expression(pi),
                              expression(gamma[1]), expression(gamma[2]), expression(gamma[3]))) +
  xlab("Parameter") +
  ylab("Partial rank correlation coefficients") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#ea2c62", "#7579e7")) +
  ggtitle(expression("Proportion of symptomatic cattle removed from cattle batches, " *italic(P[remove]))) +
  geom_text(aes(y = 0))
G3

G4 <- ggplot(data = P4, aes(x = param, y = gamma, fill = col, label = sig, hjust = hjust)) + geom_bar(stat = "identity") +
  theme_bw() +
  coord_flip() +
  ylim(-1,1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = c(expression(beta), expression("r"), expression(delta),expression("e"["v"]),
                              expression("P"["inf"]), expression("P"["R(0)"]), expression(pi),
                              expression(gamma[1]), expression(gamma[2]), expression(gamma[3]))) +
  xlab("Parameter") +
  ylab("Partial rank correlation coefficients") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#ea2c62", "#7579e7")) +
  ggtitle(expression("Risk of infected weeks, " *italic(P["week"]))) +
  geom_text(aes(y = 0))
G4

ggarrange(G1, G4, G3, ncol = 2, nrow = 2, labels = "AUTO")
ggarrange(G1, G4, G3, G2, ncol = 2, nrow = 2, labels = "AUTO")

