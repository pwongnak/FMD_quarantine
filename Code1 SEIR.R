########################################################################################################
# P. Wongnak
# Scripts for the stochastic SEIR model for the article "A stochastic modeling study of quarantine strategies against 
# foot-and-mouth disease risks through cattle trades across the Thailand-Myanmar border"
########################################################################################################
mainDir <- "C:/" ### User define
setwd(mainDir)
########################################################################################################
library(odin)
library(dde)
library(mc2d)
library(stringr)
########################################################################################################
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
############################################################################################################
# Quarantine period (in weeks)
nWeekList <- c(1,2,3,4)
# running time per week
tmax <- 7
tt <- seq(0, tmax, length.out = tmax+1)
time <- 0:tmax

# Transmission rate
b <- 0.01 #transmission rate (density-dependent)
RList <- seq(0,1,0.1) #reduction factor for between batch transmission


# reduction of transmission from asymptpmatic cattle
delta_asymptomatic <- 1
# vaccine_effectiveness
vaccine_effectivenessList <- seq(0,1, 0.1)
# Probability of imported cattle being FMDV infected; P_inf
prev <- 0.1
# Proportions of immune cattle upon arrival; Imm_1(0)
ImmList <- seq(0,1, 0.1) 
# Proportion of infectious individual becoming symptomatic; pi
pi <- 0.68
# Importing cattle ################################################
Import_cattle <- function(Imm, prev){
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

#################################################################################################
nrep <- 100
for(nWeek in nWeekList){
  for(Imm in ImmList){
    for(r in RList){
     for(vaccine_effectiveness in vaccine_effectivenessList){
    
    #transmission   
    br <- b * r
    beta <- matrix(b, nrow = nWeek, ncol = nWeek)
    for(i in 1:nWeek){
      for(j in 1:nWeek){
        if (i != j){
          beta[i,j] <- br
        }
      }
    }

    # Vaccine-induced immunity
    immune_v_function <- function(time, t_imm_min = 4, t_imm_max = 11){
      output <- (punif(time+1, t_imm_min, t_imm_max) - punif(time, t_imm_min, t_imm_max)) / (1 - punif(time, t_imm_min, t_imm_max))
      output[which(is.nan(output))] <- 0
      return(output)
    } # Function of vaccine-induced immunization
    p_immuned <- matrix(NA, nrow = length(time), ncol = nWeek)
    for(i in 1:nWeek){
      p_immuned[,i] <- immune_v_function(time + (7*(i-1)))
    }
    
    marketList <-NULL
    outputALL <-NULL
    for(j in 1:nrep){
      init_zero <- rep(0, (nWeek-1))
      ## initiate the model
      imported_cattle <- Import_cattle(Imm, prev)
      initS <- c(imported_cattle$cattleS, init_zero)
      initE <- c(imported_cattle$cattleE, init_zero)
      initIA1 <- c(imported_cattle$cattleIA1, init_zero)
      initIA2 <- c(imported_cattle$cattleIA2, init_zero)
      initIs <- c(imported_cattle$cattleIs, init_zero)
      initR <- c(imported_cattle$cattleR, init_zero)

      maxit <- 52 # maximum iteration
      ToMarket <- NULL
      AllResults <- NULL
      
      for(i in 1:maxit){
        week <- i
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
                       vaccine_effectiveness = vaccine_effectiveness
        )
        
        #model
        mod <- seir$new(user = params)
        #output
        (y <- data.frame(mod$run(tt)))
        
        if(i == 1){y = y}else{y = y[-1,]}
        y$day <- y$step + (7*(i-1))
        y <- cbind("rep" = j, y)
        AllResults <- rbind(AllResults, y)
        
        flagLastWeek <- which(str_detect(colnames(y), paste0(".", nWeek, ".")))
        ToMarket_sub <- cbind("rep" = j, "week" = week, y[nrow(y),flagLastWeek])
        ToMarket <- rbind(ToMarket, ToMarket_sub)
        
        #changing the class
        S_end <- c(as.matrix(y[nrow(y), which(str_detect(colnames(y), paste0("S.")))]))
        E_end <- c(as.matrix(y[nrow(y), which(str_detect(colnames(y), paste0("E.")))]))
        IA1_end <- c(as.matrix(y[nrow(y), which(str_detect(colnames(y), paste0("IA1.")))]))
        IA2_end <- c(as.matrix(y[nrow(y), which(str_detect(colnames(y), paste0("IA2.")))]))
        Is_end <- c(as.matrix(y[nrow(y), which(str_detect(colnames(y), paste0("Is.")))]))
        R_end <- c(as.matrix(y[nrow(y), which(str_detect(colnames(y), paste0("R.")))]))
        
        #importing the next batch
        imported_cattle <- Import_cattle(Imm, prev)
        
        initS <- c(imported_cattle$cattleS, S_end[-nWeek])
        initE <- c(imported_cattle$cattleE, E_end[-nWeek])
        initIA1 <- c(imported_cattle$cattleIA1, IA1_end[-nWeek])
        initIA2 <- c(imported_cattle$cattleIA2, IA2_end[-nWeek])
        initIs <- c(imported_cattle$cattleIs, Is_end[-nWeek])
        initR <- c(imported_cattle$cattleR, R_end[-nWeek])
      }
      marketList <-rbind(marketList, ToMarket)
      outputALL <-rbind(outputALL, AllResults)
      print(c(Imm, r, nWeek,vaccine_effectiveness, j))
    }
    #Exporting the results
    dir.create(file.path(mainDir, "mainResults_new"), showWarnings = FALSE)
    setwd(file.path(mainDir, "mainResults_new"))
    write.csv(marketList, paste0("marketList_Imm_", Imm, "_r_", r, "_vac_", vaccine_effectiveness, "_", nWeek, "w.csv"))
    write.csv(outputALL, paste0("outputALL_Imm_", Imm, "_r_",  r, "_vac_", vaccine_effectiveness, "_", nWeek, "w.csv"))
   }
  }
  }
}

