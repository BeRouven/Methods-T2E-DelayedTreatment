################################################################################
####       R-Script to generate data for type 1 error assessment            ####
################################################################################

library(tidyverse)
library(rpact)

# Create dataframe that contains all parameter constellations -------------
tau <- c(12,24,48,60)
acc.prop <- 0.2 * 1:2
lag2 <- 0
lag1 <- 0
theta <- 0.5
kC <- c(0.5,1,2)
medC <- c(5, 15, 20)
params <- expand.grid(
  tau = tau, acc = acc.prop, lag2 = lag2, lag1 =lag1, theta = theta,
  kC = kC, medC = medC)
#Calculate number of observations with rpact
set.seed(1042024)
params %>%
  rowwise() %>%
  mutate(nobs = 2*ceiling(
    rpact::getSampleSizeSurvival(alpha=0.05, sided=2, beta=0.2,
                                 median2 = medC, hazardRatio = theta,
                                 kappa = kC,
                                 followUpTime = (1-acc)*tau,
                                 accrualTime = c(0,tau*acc),
                                 typeOfComputation = "Schoenfeld")$nFixed1),
    nevents = ceiling(
      rpact::getSampleSizeSurvival(alpha=0.05, sided=2, beta=0.2,
                                   median2 = medC, hazardRatio = theta,
                                   kappa = kC,
                                   followUpTime = (1-acc)*tau,
                                   accrualTime = c(0,tau*acc),
                                   typeOfComputation = "Schoenfeld")$eventsFixed)) %>%
  ungroup() %>%
  mutate(nsim=2500,
         seed = sample(10000:99999999, 72)) %>%
  as.data.frame() ->
params

# Generate data for each power scenario -----------------------------------

library(tidyverse)
library(parallel)
library(doMC)
library(foreach)

# Simulate datasets and save to disk
core.number <- 20
registerDoMC(core.number)
foreach (i = 1:nrow(params),
         .packages = c("stats", "tidyverse"),
         .combine = rbind
) %dopar% {
  #Set path for output data
  out.path <- "./T1E/SimulationData/"
  
  #Set seed from params dataframe for simulations
  set.seed(params[i, "seed"])
  #Derive parameters from params dataframe
  nobs <- params[i, "nobs"]
  kC <- params[i, "kC"] 
  medC <- params[i, "medC"]
  lambdaC <- log(2)^(1/kC)/medC
  tau <- params[i, "tau"]
  acc <- params[i, "acc"]
  #Initialize list to save simulated datasets
  ls <- list()
  #Simulate data
  for (j in 1:params[i, "nsim"]){
    df <- data.frame(
      group = rep(0:1, each = nobs/2),
      fail.time = rweibull(nobs, shape = kC, scale = 1/lambdaC),
      acc.time = runif(nobs, min = 0, max = acc * tau)
    )
    
    df %>%
      mutate(
        obs.time = pmin(fail.time, params[i, "tau"]-acc.time),
        event = ifelse(fail.time<= params[i, "tau"]-acc.time, 1, 0)
      ) ->
      ls[[j]]
  }
  saveRDS(ls, file = paste0(out.path,
                            "tau_", tau,
                            "_acc_", acc,
                            "_medC_", medC,
                            "_kC_", kC,
                            "_nobs_", nobs, ".RData"))
  rm(ls)
}
