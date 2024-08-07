################################################################################
####           R-Script to generate data for power assessment               ####
################################################################################

library(tidyverse)
library(rpact)

# Create dataframe that contains all parameter constellations -------------
tau <- c(12,24,48,60)
acc.prop <- 0.2 * 1:2
lag2 <- 0.1 * 0:4
lag1 <- c(0, 0.3, 0.7, 1)
theta <- 0.1 * 5:8
kC <- c(0.5,1,2)
medC <- c(5, 15, 20)
params <- expand.grid(
  tau = tau, acc = acc.prop, lag2 = lag2, lag1 =lag1, theta = theta,
  kC = kC, medC = medC) %>%
  filter((lag2==0 & lag1==0) | (lag2 != 0 & lag1 !=0))
#Calculate number of observations with rpact
set.seed(67057533)
params %>%
  mutate(t1star = lag1*lag2*tau,
         t2star = lag2*tau) %>%
  mutate(naive.avgHR = 1 + ((theta - 1)* (t2star-t1star))/(2*tau) + ((theta-1)*(tau-t2star))/tau) %>%
  rowwise() %>%
  mutate(naive.nobs.avg = 2*ceiling(
    rpact::getSampleSizeSurvival(alpha=0.05, sided=2, beta=0.2,
                                 median2 = medC, hazardRatio = naive.avgHR,
                                 kappa = kC,
                                 followUpTime = (1-acc)*tau,
                                 accrualTime = c(0,tau*acc),
                                 typeOfComputation = "Schoenfeld")$nFixed1),
    naive.nevents.avg = ceiling(
      rpact::getSampleSizeSurvival(alpha=0.05, sided=2, beta=0.2,
                                   median2 = medC, hazardRatio = naive.avgHR,
                                   kappa = kC,
                                   followUpTime = (1-acc)*tau,
                                   accrualTime = c(0,tau*acc),
                                   typeOfComputation = "Schoenfeld")$eventsFixed)) %>%
  ungroup() %>%
  mutate(nsim=2500,
         seed = sample(10000:99999999, 3744)) ->
params

# Generate data for each power scenario -----------------------------------

library(tidyverse)
library(parallel)
library(doMC)
library(foreach)

params %>%
  arrange(desc(naive.nobs.avg)) %>%
  as.data.frame() ->
params

# Simulate datasets and save to disk
core.number <- 20
registerDoMC(core.number)
foreach (i = 1:nrow(params),
         .packages = c("stats", "tidyverse"),
         .combine = rbind
) %dopar% {
  source("./0_CostumeFunctions_DataGeneration.R")
  #Set path for output data
  out.path <- "./Power/SimulationData/"
  
  #Set seed from params dataframe for simulations
  set.seed(params[i, "seed"])
  #Derive parameters from params dataframe
  nobs <- params[i, "naive.nobs.avg"]
  kC <- params[i, "kC"] 
  medC <- params[i, "medC"]
  lambdaC <- log(2)^(1/kC)/medC
  tau <- params[i, "tau"]
  acc <- params[i, "acc"]
  acc.time <- acc * tau
  lag2 <- params[i, "lag2"]
  lag1 <- params[i, "lag1"]
  t2star <- lag2 * tau
  t1star <- lag1 * t2star
  theta <- params[i, "theta"]
  #Initialize list to save simulated datasets
  ls <- list()
  #Simulate data
  for (j in 1:params[i, "nsim"]){
    if (t2star==0){
      df <- data.frame(
        group = rep(0:1, each = nobs/2),
        fail.time = c(rweibull(nobs/2, shape = kC, scale = 1/lambdaC),
                      rweibull(nobs/2, shape = kC, scale = 1/(theta^(1/kC)*lambdaC))),
        acc.time = runif(nobs, min = 0, max = params[i, "acc"] * params[i, "tau"])
      )
    } else {
      df <- data.frame(
        group = rep(0:1, each = nobs/2),
        fail.time = c(rweibull(nobs/2, shape = kC, scale = 1/lambdaC),
                      vinvLambdaE(-log(1-runif(nobs/2, min = 0, max = 1)),
                             theta = theta, t1star = t1star, t2star = t2star, 
                             lambdaC = lambdaC, kC = kC)),
        acc.time = runif(nobs, min = 0, max = params[i, "acc"] * params[i, "tau"])
      )
    }
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
                            "_HR_", theta,
                            "_medC_", medC,
                            "_kC_", kC,
                            "_lag2_", lag2,
                            "_lag1_", lag1, ".RData"))
  rm(ls)
}
