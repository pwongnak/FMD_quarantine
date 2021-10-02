#############################################################
# Modelling FMD risk through quarantine center
# SEIRV model (Base model)
# P. WONGNAK; Last updated; 11 September 2020
# Loading packages ##########################################
library(mc2d)
library(deSolve)

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
# Duration of quarantine; T (in weeks)
Duration <- c(1,2,3,4)
# Proportions of immune cattle upon arrival; Imm_1(0)
Imm <- seq(0.025,0.9, 0.05) 
# Other parameters was placed in the loops to make it easier when performing sensitivity analysis
# SEIRV model ###############################################
for(d in 1:length(Duration)){
  for(I in 1:length(Imm)){
    # Defining parameters -------------------------------------------
    duration <- Duration[d]
    imm <- Imm[I]
    # Within-herd transmission coefficient; beta
    beta <- matrix(NA, ncol = duration, nrow = duration)
    beta[,] <- 0.010 # Density-dependent 
    # Reduction factor for between-batch transmission; r
    red <- 1
    # To make different beta between batches
    for(i in 1:duration){
      for(j in 1:duration){
        if (i != j){
          beta[i,j] <- red*beta[i,j]
        }
      }
    }
    # Proportional infectiousness of asymptomatic individuals; delta
    delta <- 0.99
    # Latent period; t_lat (in days)
    latent <- 3.6
    # Infectious period; t_inf (in days)
    infectiousP <- 4.4
    # Incubation period; t_inc (in days)
    incubation <- 5.9
    # Vaccine protection rate; 1/(time to immunity after vaccination); alpha
    alpha <- 1/28
    # Transition rate from asymptomatic to symptomatic; sigma
    sigma <- 1/latent
    # Proportion of infectious individual becoming symptomatic; pi
    pi <- 0.68
    # Recovery rate of asymptomatic individuals; gamma_1 (gammaAR)
    gammaAR <- 1/infectiousP
    # Transition rate from asymptomatic to symptomatic individuals; gamma_2 (gammaAS)
    gammaAS <- 1/(incubation-latent)
    # Recovery rate of symptomatic individuals; gamma_3 (gammaSR)
    gammaSR <- 1/(infectiousP-(incubation-latent))
    # Probability of imported cattle being FMDV infected; P_inf
    prev <- 0.1
    
    # SEIRV model function ---------------------------------------------------
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
    
    # Defining index for each compartment; depending of the duration of quarantine
    sindex <- 1:duration
    eindex <- (duration+1):(2*duration)
    iaindex <- (2*duration+1):(3*duration)
    isindex <- (3*duration+1):(4*duration)
    rindex <-  (4*duration+1):(5*duration)
    vindex <- (5*duration+1):(6*duration)
    
    # Number of interations
    rep <- 100
    
    marketList <- NULL
    outputALL <- NULL
    for(r in 1:rep){
      # Intitial states-------------------------------------
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
      
      runtime <- 52 # 52 weeks
      
      market <- matrix(NA, nrow = runtime, ncol = 8)
      for(i in 1:runtime){
        # Solving SEIRV model for 7 days
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
        
        # Importing new cattle batch
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
      print(r)
    }
    #Exporting the results
    dir.create(file.path(mainDir, "mainResults"), showWarnings = FALSE)
    setwd(file.path(mainDir, "mainResults"))
    write.csv(marketList, paste0("marketList_", imm, "_", duration, "w.csv"))
    write.csv(outputALL, paste0("outputALL_", imm, "_", duration, "w.csv"))
  }
}





