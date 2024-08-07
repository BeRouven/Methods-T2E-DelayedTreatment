################################################################################
####       R-Script to analyse data for type 1 error assessment            ####
################################################################################

library(tidyverse)
library(rpact)
library(doFuture)
library(foreach)
library(progressr)

# Create dataframe that contains all parameter constellations -------------
tau <- c(12,24,48,60)
acc.prop <- 0.2 * 1:2
lag2 <- 0.1 * 0:4
lag1 <- c(0, 0.3, 0.7, 1)
theta <- 0.1 * 5:8
kC <- c(0.5,1,2)
medC <- c(5, 15, 20)
expand.grid(
  tau = tau, acc = acc.prop, lag2 = lag2, lag1 =lag1, theta = theta,
  kC = kC, medC = medC) %>%
  filter((lag2==0 & lag1==0) | (lag2 != 0 & lag1 !=0)) %>%
  mutate(t1star = lag1*lag2*tau,
         t2star = lag2*tau,
         naive.avgHR = 1 + ((theta - 1)* (t2star-t1star))/(2*tau) + ((theta-1)*(tau-t2star))/tau) %>%
  rowwise() %>%
  mutate(naive.nobs.avg = 2*ceiling(
    rpact::getSampleSizeSurvival(alpha=0.05, sided=2, beta=0.2,
                                 median2 = medC, hazardRatio = naive.avgHR,
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
  in.path <- "./Power/SimulationData/"
  #Derive parameters from params dataframe
  tau <- params[i, "tau"]
  acc <- params[i, "acc"]
  theta <- params[i, "theta"]
  medC <- params[i, "medC"]
  kC <- params[i, "kC"]
  lag2 <- params[i, "lag2"]
  lag1 <- params[i, "lag1"]
  
  #Read in simulated datasets
  ls <- readRDS(ls, file = paste0(in.path,
                            "tau_", tau,
                            "_acc_", acc,
                            "_HR_", theta,
                            "_medC_", medC,
                            "_kC_", kC,
                            "_lag2_", lag2,
                            "_lag1_", lag1, ".RData"))
  #Calculate the number of datasets in which the milestone rate is estimable
  max.time <- c()
  for (j in 1:length(ls)){
    df <- ls[[j]]
    #Calculate minimum of maximum time of both groups
    max.time[j] <- min(max(df[df$group==0, "obs.time"]), max(df[df$group==1, "obs.time"]))
  }
  
  c(tau, acc, theta, medC, kC, lag2, lag1, min(max.time))
} ->
  max.CO.evaluable

colnames(max.CO.evaluable) <- c("tau", "acc", "theta", "medC", "kC", "lag2", "lag1", "MinMaxTime")

params %>%
  left_join(., as.data.frame(max.CO.evaluable), by = c("tau", "acc", "theta", "medC", "kC", "lag2", "lag1")) %>%
  mutate(CO.RMST = floor(MinMaxTime),
         CO.mile = floor(MinMaxTime),
         CO.land1 = ifelse(lag2==0, 0.1*tau, 0.9 * lag1 *lag2 * tau),
         CO.land2 = ifelse(lag2==0, 0.2*tau, lag1 *lag2 * tau),
         CO.land3 = ifelse(lag2==0, 0.3*tau, 1.1 * lag1 *lag2 * tau),
         CO.genlin.low1 = ifelse(lag2==0, 0*tau,
                                 ifelse(lag1==1, 0.9* lag1 * lag2 * tau,
                                        lag2*tau*(0.5*(1+lag1) - 1.1 * 0.5 * (1-lag1)))), 
         CO.genlin.up1 = ifelse(lag2==0, 0.3*tau,
                                ifelse(lag1==1, 1.1* lag1 * lag2 * tau,
                                       lag2*tau*(0.5*(1+lag1) + 1.1 * 0.5 * (1-lag1)))),
         CO.genlin.low2 = ifelse(lag2==0, 0*tau,
                                 ifelse(lag1==1, 0.8* lag1 * lag2 * tau,
                                        lag1*lag2*tau)), 
         CO.genlin.up2 = ifelse(lag2==0, 0.2*tau,
                                ifelse(lag1==1, 1.2* lag1 * lag2 * tau,
                                       lag2*tau)),
         CO.genlin.low3 = ifelse(lag2==0, 0.1*tau,
                                 ifelse(lag1==1, 0.7* lag1 * lag2 * tau,
                                        lag2*tau*(0.5*(1+lag1) - 0.9 * 0.5* (1-lag1)))), 
         CO.genlin.up3 = ifelse(lag2==0, 0.3*tau,
                                ifelse(lag1==1, 1.3* lag1 * lag2 * tau,
                                       lag2*tau*(0.5*(1+lag1) + 0.9 * 0.5 * (1-lag1))))
         ) ->
params.analysis

# Analyze datasets and save to disk ---------------------------------------
#Due to the big sample size in some scenarios the power calculations is not
#parallelized over the simulation scenarios but over the simulated datasets
#in each scenario
for (scen.i in 1:nrow(params)){
  #Set path to read in data
  in.path <- "./Power/SimulationData/"
  #Derive parameters from params dataframe
  kC <- params.analysis[scen.i, "kC"] 
  medC <- params.analysis[scen.i, "medC"]
  tau <- params.analysis[scen.i, "tau"]
  acc <- params.analysis[scen.i, "acc"]
  lag2 <- params.analysis[scen.i, "lag2"]
  lag1 <- params.analysis[scen.i, "lag1"]
  theta <- params.analysis[scen.i, "theta"]
  #Read in data
  ls <- readRDS(paste0(in.path,
                       "tau_", tau,
                       "_acc_", acc,
                       "_HR_", theta,
                       "_medC_", medC,
                       "_kC_", kC,
                       "_lag2_", lag2,
                       "_lag1_", lag1, ".RData"))
  #Save each simulated dataset to disk
  for (df.i in 1:length(ls)){
    saveRDS(ls[df.i], file=paste0(in.path,
                                  "BD/tau_", tau,
                                  "_acc_", acc,
                                  "_HR_", theta,
                                  "_medC_", medC,
                                  "_kC_", kC,
                                  "_lag2_", lag2,
                                  "_lag1_", lag1,
                                  "_", df.i, ".RData"))
  }
  rm(ls, df.i)
  
  #Assign parameters for analysis
  Landmark.CO <- unlist(params.analysis[scen.i, paste0("CO.land", 1:3)])
  RMST.CO <- unlist(params.analysis[scen.i,"CO.RMST"])
  Milestone.CO <- unlist(unique(params.analysis[scen.i,"CO.mile"]))
  GenLin.CO.low <- unlist(params.analysis[scen.i, paste0("CO.genlin.low", 1:3)])
  GenLin.CO.up <- unlist(params.analysis[scen.i, paste0("CO.genlin.up", 1:3)])
  
  handlers(global=TRUE)
  handlers("progress")
  plan(multisession, workers = 20)
  
  my_fcn <- function(xs){
    p <- progressor(along=xs)
    foreach (df.i = xs,
             .options.future = list(packages = c("stats", "survival", "mnormt", 
                                                 "timereg", "tidyverse")),
             .combine = rbind
    ) %dofuture% {
      source("./0_CostumeFunctions_Analysis.R", local=TRUE)
      
      df.sim <- as.data.frame(readRDS(file=paste0(in.path,
                                                  "BD/tau_", tau,
                                                  "_acc_", acc,
                                                  "_HR_", theta,
                                                  "_medC_", medC,
                                                  "_kC_", kC,
                                                  "_lag2_", lag2,
                                                  "_lag1_", lag1,
                                                  "_", df.i, ".RData")))
    
      res.temp <- MethodsSurvival(data = df.sim,
                                  vars = list(time="obs.time", event = "event", group = "group"),
                                  params = list(AFTdist = c("weibull", "exponential", "gaussian", 
                                                            "logistic", "lognormal", "loglogistic"),
                                                FH = list(rho = c(0,0,1,1,-1,0,0,0.5),
                                                          gamma = c(0,1,0,1,0,0.5,2,0.5)),
                                                mWLR.Delay = Landmark.CO,
                                                RMST.cutoff = RMST.CO, 
                                                milestone.cutoff = Milestone.CO,
                                                landmark.cutoff = Landmark.CO,
                                                gen.lin = list(t.low = GenLin.CO.low, 
                                                               t.up = GenLin.CO.up),
                                                PW.exp.lag = list(lag = GenLin.CO.low, 
                                                                  CP = GenLin.CO.up),
                                                threshold = Landmark.CO,
                                                ABC = RMST.CO,
                                                PW.exp = Landmark.CO,
                                                V0.thresh = Landmark.CO, 
                                                MERT = list(t.low = GenLin.CO.low, 
                                                            t.up = GenLin.CO.up),
                                                par.group = Landmark.CO,
                                                logit = list(w.low = 0.1, 
                                                             t.low = GenLin.CO.low, 
                                                             t.up = GenLin.CO.up)
                                         )
      )
      p()
      res.temp
    } ->
      res.sim
    return(res.sim)
  }
  
  res.sim1 <- my_fcn(1:2500)
  out.path <- "./Power/Results/"
  out.file <- paste0(out.path, "tau_", tau, "_acc_", acc, "_HR_", theta, 
                     "_medC_", medC, "_kC_", kC, "_lag2_", lag2, "_lag1_", 
                     lag1, ".csv")
  
  write.table(res.sim1, file = out.file, sep = ";", append = TRUE, col.names = FALSE, row.names = FALSE)
  
  #Delete single datasets of this scenario
  for (df.i in 1:2500){
    file.remove(paste0(in.path,
                       "BD/tau_", tau,
                       "_acc_", acc,
                       "_HR_", theta,
                       "_medC_", medC,
                       "_kC_", kC,
                       "_lag2_", lag2,
                       "_lag1_", lag1,
                       "_", df.i, ".RData"))
  }
  rm(df.i, res.sim1)
  closeAllConnections()
}