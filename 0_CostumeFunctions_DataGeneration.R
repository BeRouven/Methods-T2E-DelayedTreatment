################################################################################
####           Costume functions needed for data generation                 ####
################################################################################

#-------------- Cumulative hazard function in experimental arm -----------------
#------------------ based on generalized linear lag model ----------------------
LambdaE <- function(t, theta, t1star, t2star, lambdaC, kC){
  ## Parameters:
  # t:            input parameter
  # theta:        full treatment effect
  # t1star:       changepoint
  # t2star:       delay
  # lambdaC:      scale parameter of the Weibull hazard in control arm
  # kC:           shape of the Weibull hazard in control arm
  if (t1star == t2star){
    if (t < t1star){
      res <- (lambdaC * t)^kC
    } else {
      res <- (1 - theta) * (lambdaC * t1star)^kC + theta * (lambdaC * t)^kC
    }
  } else {
    B <- (theta - 1) / (t2star - t1star)
    A <- 1 - B * t1star
    if (t < t1star){
      res <- (lambdaC * t)^kC
    } else if (t < t2star){
      res <- (A + B * t) * (lambdaC * t)^kC - B / (kC + 1) * (t * (lambdaC * t)^kC - t1star * (lambdaC * t1star)^kC)
    } else {
      res <- theta * (lambdaC * t)^kC - B / (kC + 1) * (t2star * (lambdaC * t2star)^kC - t1star * (lambdaC * t1star)^kC)
    }
  }
  return(res)
}
vLambdaE <- Vectorize(LambdaE, "t")

#----------- Inverse cumulative hazard function in experimental arm ------------
#------------------ based on generalized linear lag model ----------------------
invLambdaE <- function(z, theta, t1star, t2star, lambdaC, kC){
  ## Parameters:
  # z:            input parameter
  # theta:        full treatment effect
  # t1star:       changepoint
  # t2star:       delay
  # lambdaC:      scale parameter of the Weibull hazard in control arm
  # kC:           shape of the Weibull hazard in control arm
  if (t1star == t2star){
    if (z < (lambdaC * t1star)^kC){
      res <- 1 / lambdaC * z^(1/kC)
    } else {
      res <- 1 / lambdaC * ((z - (1 - theta) * (lambdaC * t1star)^kC)/theta)^(1/kC)
    }
  } else {
    B <- (theta - 1) / (t2star - t1star)
    A <- 1 - B * t1star
    CO1 <- (lambdaC * t1star)^kC
    CO2 <- theta * (lambdaC * t2star)^kC - B / (kC + 1) * (t2star * (lambdaC * t2star)^kC - t1star * (lambdaC * t1star)^kC)
    if (z < CO1){
      res <- 1 / lambdaC * z^(1/kC)
    } else if (z < CO2){
      res <- uniroot(function(x) LambdaE(x, theta = theta, t1star = t1star, t2star = t2star, lambdaC = lambdaC, kC = kC) - z,
                     lower = t1star, upper = t2star, extendInt = "no")$root
    } else {
      res <- 1 / lambdaC * ((z + B / (kC + 1) * (t2star * (lambdaC * t2star)^kC - t1star * (lambdaC * t1star)^kC)) / theta)^(1/kC)
    }
  }
  return(res)
}
vinvLambdaE <- Vectorize(invLambdaE, "z")