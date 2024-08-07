################################################################################
####       R-Script to analyse data for type 1 error assessment            ####
################################################################################

library(tidyverse)
library(rpact)
library(parallel)
library(doFuture)
library(foreach)

# Create dataframe that contains all parameter constellations -------------
tau <- c(12,24,48,60)
acc.prop <- 0.2 * 1:2
lag2 <- 0
lag1 <- 0
theta <- 0.5
kC <- c(0.5,1,2)
medC <- c(5, 15, 20)
expand.grid(
  tau = tau, acc = acc.prop, lag2 = lag2, lag1 =lag1, theta = theta,
  kC = kC, medC = medC) %>%
  rowwise() %>%
  mutate(nobs = 2*ceiling(
    rpact::getSampleSizeSurvival(alpha=0.05, sided=2, beta=0.2,
                                 median2 = medC, hazardRatio = theta,
                                 kappa = kC,
                                 followUpTime = (1-acc)*tau,
                                 accrualTime = c(0,tau*acc),
                                 typeOfComputation = "Schoenfeld")$nFixed1)) %>%
  ungroup() %>%
  as.data.frame() ->
params

#Assess the maximal CO for RMST and Milestone survival for which all simulated data are evaluable
core.number <- 10
plan(multisession, workers = core.number)
foreach (i = 1:nrow(params),
         .options.future = list(packages = c("stats", "tidyverse")),
         .combine = rbind
) %dofuture% {
  #Set path for output data
  in.path <- "./T1E/SimulationData/"
  #Derive parameters from params dataframe
  kC <- params[i, "kC"]
  medC <- params[i, "medC"]
  tau <- params[i, "tau"]
  acc <- params[i, "acc"]
  nobs <- params[i, "nobs"]
  
  #Read in simulated datasets
  ls <- readRDS(file = paste0(in.path,
                              "tau_", tau,
                              "_acc_", acc,
                              "_medC_", medC,
                              "_kC_", kC,
                              "_nobs_", nobs, ".RData"))
  #Calculate the number of datasets in which the milestone rate is estimable
  max.time <- c()
  for (j in 1:length(ls)){
    df <- ls[[j]]
    #Calculate minimum of maximum time of both groups
    max.time[j] <- min(max(df[df$group==0, "obs.time"]), max(df[df$group==1, "obs.time"]))
  }

  c(kC, medC, tau, acc, nobs, min(max.time))
} ->
max.CO.evaluable

colnames(max.CO.evaluable) <- c("kC", "medC", "tau", "acc", "nobs", "MinMaxTime")

params %>%
  left_join(., as.data.frame(max.CO.evaluable), by = c("kC", "medC", "tau", "acc", "nobs")) %>%
  mutate(CO.RMST = floor(MinMaxTime), 
         CO.mile = floor(MinMaxTime),
         CO.land = 0.2*tau,
         CO.genlin.low = 0, 
         CO.genlin.up = 0.2*tau) ->
params


# Analyze datasets and save to disk ---------------------------------------
core.number <- 16
plan(multisession, workers = core.number)
foreach (i = 1:nrow(params),
         .options.future = list(packages = c("stats")),
         .combine = append
) %dofuture% {
  source("./0_CostumeFunctions_Analysis.R", local = TRUE)
  
  #Set path to read in data
  in.path <- "./T1E/SimulationData/"
  #Derive parameters from params dataframe
  nobs <- params[i, "nobs"]
  kC <- params[i, "kC"] 
  medC <- params[i, "medC"]
  tau <- params[i, "tau"]
  acc <- params[i, "acc"]
  ls <- readRDS(paste0(in.path,
                       "tau_", tau,
                       "_acc_", acc,
                       "_medC_", medC,
                       "_kC_", kC,
                       "_nobs_", nobs, ".RData"))
  
  out.path <- "./T1E/Results/"
  out.file <- paste0(out.path,
                     "tau_", tau,
                     "_acc_", acc,
                     "_medC_", medC,
                     "_kC_", kC, 
                     "_nobs_", nobs, ".csv")
  for (df.i in 1:length(ls)){
    res <- MethodsSurvival(data = ls[[df.i]],
                           vars = list(time="obs.time", event = "event", group = "group"),
                           params = list(AFTdist = c("weibull", "exponential", "gaussian", 
                                                     "logistic", "lognormal", "loglogistic"),
                                         FH = list(rho = c(0,0,1,1,-1,0,0,0.5),
                                                   gamma = c(0,1,0,1,0,0.5,2,0.5)),
                                         mWLR.Delay = params[i, "CO.land"],
                                         RMST.cutoff = params[i,"CO.RMST"], 
                                         milestone.cutoff = unlist(params[i,"CO.mile"]),
                                         landmark.cutoff = params[i, "CO.land"],
                                         gen.lin = list(t.low = params[i, "CO.genlin.low"],
                                                        t.up = params[i, "CO.genlin.up"]),
                                         PW.exp.lag = list(lag = params[i, "CO.genlin.low"], 
                                                           CP = params[i, "CO.genlin.up"]),
                                         threshold = params[i, "CO.land"],
                                         ABC = unlist(params[i,"CO.RMST"]),
                                         PW.exp = params[i, "CO.land"],
                                         V0.thresh = params[i, "CO.land"], 
                                         MERT = list(t.low = params[i, "CO.genlin.low"], 
                                                     t.up = params[i, "CO.genlin.up"]),
                                         par.group = params[i, "CO.land"],
                                         logit = list(w.low = 0.1, 
                                                      t.low = params[i, "CO.genlin.low"], 
                                                      t.up = params[i, "CO.genlin.up"])
                                  )
    )
    write.table(res, file = out.file, sep = ";", append = TRUE, col.names = FALSE, row.names = FALSE)
    
  }
}