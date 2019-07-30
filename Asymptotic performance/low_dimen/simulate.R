### Figure 1 ###
### low-dimensional setting ###

library(ggplot2)
library(gridExtra)
library(cowplot)
library(MASS)


True_Y <- function(x,z){
  y <- sin(3*x) + z
  return(y)
}
True_Z <- function(x){
  z <- cos(2*x)
  return(z)
}

e_mu <- c(0,0)
e_sigma <- matrix(c(1,-0.5,-0.5,1),2,2)

Sim_Cov_YZ <- function(n, True_Y, True_Z, e_sigma, e_mu, iter, h_y, h_z, h_yz, true_Cov_YZ)
{
  Naive1 <- NULL
  Naive2 <- NULL
  Onestep <- NULL
  Onestep1.CI <- NULL
  Onestep2.CI <- NULL
  Naive1.CI <- NULL
  Naive2.CI <- NULL
  
  for(b in 1:iter){
    X <- runif(n, min = 0, max = 2)
    e <- mvrnorm(n, e_mu, e_sigma)
    mu_Z <- True_Z(X); Z <- mu_Z + e[,1]
    mu_Y <- True_Y(X,Z); Y <- mu_Y + e[,2]
    
    # Estimate E(Y|X), E(Z|X), and $E(YZ|X)$
    Hat.muY <- loess(Y ~ X, span = h_y, degree = 1, control=loess.control(surface="direct"))$fitted
    Hat.muZ <- loess(Z ~ X, span = h_z, degree = 1, control=loess.control(surface="direct"))$fitted
    Hat.muYZ <- loess(Y*Z ~ X, span = h_yz, degree = 1, control=loess.control(surface="direct"))$fitted
    
    naive1 <- mean(Hat.muYZ-Hat.muY*Hat.muZ)
    naive2 <- cov(Y,Z) - cov(Hat.muY, Hat.muZ)
    
    # bootstrap
    naive1.idx <- NULL
    naive2.idx <- NULL
    for (i in 1:400) {
      idx <- sample(1:n, n, replace = T)

      Hat.muY.idx <- loess(Y[idx] ~ X[idx], span = h_y, degree = 1, control=loess.control(surface="direct"))$fitted
      Hat.muZ.idx <- loess(Z[idx] ~ X[idx], span = h_z, degree = 1, control=loess.control(surface="direct"))$fitted
      Hat.muYZ.idx <- loess(Y[idx]*Z[idx] ~ X[idx], span = h_yz, degree = 1, control=loess.control(surface="direct"))$fitted

      naive1.idx <- c(naive1.idx, mean(Hat.muYZ.idx-Hat.muY.idx*Hat.muZ.idx))
      naive2.idx <- c(naive2.idx, cov(Y[idx],Z[idx]) - cov(Hat.muY.idx, Hat.muZ.idx))
    }
    
    # Naive estimator of marginal Covrelation
    Naive1 <- c(Naive1, naive1)
    Naive2 <- c(Naive2, naive2)
    
    naive1.ci <- ifelse(true_Cov_YZ<= quantile(naive1.idx, probs = 0.975) & true_Cov_YZ >= quantile(naive1.idx, probs = 0.025), 1, 0)
    naive2.ci <- ifelse(true_Cov_YZ<= quantile(naive2.idx, probs = 0.975) & true_Cov_YZ >= quantile(naive2.idx, probs = 0.025), 1, 0)

    Naive1.CI <- c(Naive1.CI, naive1.ci)
    Naive2.CI <- c(Naive2.CI, naive2.ci)
    
  
    # One step estimator of marginal Covariance
    Dp = Y*Z - Hat.muZ*Y - Hat.muY*Z + Hat.muY*Hat.muZ
    Dp1 <- mean((Dp-naive1)^2)

    onestep <- mean(Dp)
    Onestep <- c(Onestep, onestep)
    
    # Standard deviance of one step estimator
    onestep1.ci <- ifelse(true_Cov_YZ<= onestep + 1.96*sqrt(Dp1)/sqrt(n) & true_Cov_YZ >= onestep - 1.96*sqrt(Dp1)/sqrt(n), 1, 0)
    Onestep1.CI <- c(Onestep1.CI, onestep1.ci)
  }
  return(list(Naive1=Naive1, Naive2=Naive2, Onestep=Onestep, Naive1.CI=Naive1.CI, Naive2.CI=Naive2.CI, Onestep1.CI=Onestep1.CI))
}

size <- c(100, 500, 1000, 2000,3000, 6000)
iter = 400

# this span is pre-selected by cross-validation
span_size_y <- c(0.46, 0.46, 0.46, 0.46, 0.46, 0.46)
span_size_z <- c(0.67, 0.47, 0.37, 0.33, 0.30, 0.28)
span_size_yz <- c(0.54, 0.4, 0.37, 0.36, 0.36, 0.35)

Naive1.mat <- matrix(nrow = iter, ncol = length(size))
Naive2.mat <- matrix(nrow = iter, ncol = length(size))
Onestep.mat <- matrix(nrow = iter, ncol = length(size))
CI.Onestep1.mat <- matrix(nrow = iter, ncol = length(size))
CI.Naive1.mat <- matrix(nrow = iter, ncol = length(size))
CI.Naive2.mat <- matrix(nrow = iter, ncol = length(size))


for (i in 1:length(size)) 
{
  n <- size[i]
  Sim1 <- Sim_Cov_YZ(n, True_Y, True_Z, e_sigma, e_mu, iter = iter, h_y = span_size_y[i], h_z = span_size_z[i], h_yz = span_size_yz[i], true_Cov_YZ = 0.5)

  Naive1.mat[,i] <- Sim1$Naive1
  Naive2.mat[,i] <- Sim1$Naive2
  Onestep.mat[,i] <- Sim1$Onestep
  CI.Onestep1.mat[,i] <- Sim1$Onestep1.CI
  CI.Naive1.mat[,i] <- Sim1$Naive1.CI
  CI.Naive2.mat[,i] <- Sim1$Naive2.CI
  
}


write.table(Naive1.mat,".../low_dim/Naive1.txt")
write.table(Naive2.mat,"../low_dim/Naive2.txt")
write.table(Onestep.mat,"../low_dim/Onestep.txt")
write.table(CI.Onestep1.mat,"../low_dim/Onestep1_CI.txt")
write.table(CI.Naive1.mat,"../low_dim/Naive1_CI.txt")
write.table(CI.Naive2.mat,"../low_dim/Naive2_CI.txt")

