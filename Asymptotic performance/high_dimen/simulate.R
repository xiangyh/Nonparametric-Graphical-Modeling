### Figure 1 ###
### high-dimensional setting ###

library(MASS)
library(glmnet)
library(dplyr)

e_mu <- c(0,0)
e_sigma <- matrix(c(1,-0.5,-0.5,1),2,2)

beta_Y=c(rep(1,10),rep(0,4990))
beta_Z=c(rep(-1,10),rep(0,4990))

Sim_Cov_YZ <- function(n, beta_Y, beta_Z, e_sigma, e_mu, iter, true_Cov_YZ)
{
  Naive <- NULL
  Naive2 <- NULL
  Onestep <- NULL
  Onestep.CI <- NULL

  for(b in 1:iter){
    e <- mvrnorm(n, e_mu, e_sigma)
    X=matrix(rnorm(n*length(beta_Y)),n,length(beta_Y))
    Y <- X%*%beta_Y + e[,1]
    Z <- X%*%beta_Z + e[,2]
    YZ <- Y*Z
    
    # Estimate E(Y|X), E(Z|X), and $E(YZ|X)$
    # mean((Y-Hat.muY)^2), mean((Z-Hat.muZ)^2), mean((YZ-Hat.muYZ)^2)
    Hat.muY <- cv.glmnet(X, Y, family="gaussian", nfolds = 5)  %>% predict(s="lambda.min",newx=X) 
    Hat.muZ <- cv.glmnet(X, Z, family="gaussian", nfolds = 5) %>% predict(s="lambda.min",newx=X) 
    Hat.muYZ <- cv.glmnet(X, YZ, family="gaussian", nfolds = 5) %>% predict(s="lambda.min",newx=X) 
    # Naive estimator of marginal Covrelation
    naive <- mean(Hat.muYZ-Hat.muY*Hat.muZ)
    Naive <- c(Naive, naive)
    
    naive2 <- cov(Y,Z) - cov(Hat.muY, Hat.muZ)
    Naive2 <- c(Naive2, naive2)

    # One step estimator of marginal Covariance
    Dp = (Y-Hat.muY)*(Z - Hat.muZ)
    s2 <- mean((Dp-mean(Dp))^2)
    
    onestep <- mean(Dp)
    Onestep <- c(Onestep, onestep)
    
    # Standard deviance of one step estimator
    onestep.ci <- ifelse(true_Cov_YZ<= onestep + 1.96*sqrt(s2)/sqrt(n) & true_Cov_YZ >= onestep - 1.96*sqrt(s2)/sqrt(n), 1, 0)
    Onestep.CI <- c(Onestep.CI, onestep.ci)
  }
  return(list(Naive=Naive, Naive2=Naive2, Onestep=Onestep, Onestep.CI=Onestep.CI))
}

# size <- c(100,200)
size <- c(500, 1000, 2000, 3000, 4000, 5000, 6000)
iter <- 1000

Naive.mat <- matrix(nrow = iter, ncol = length(size))
Naive2.mat <- matrix(nrow = iter, ncol = length(size))
Onestep.mat <- matrix(nrow = iter, ncol = length(size))
CI.Onestep.mat <- matrix(nrow = iter, ncol = length(size))


for (i in 1:length(size)) 
{
  n <- size[i]
  Sim1 <- Sim_Cov_YZ(n, beta_Y, beta_Z, e_sigma, e_mu, iter, true_Cov_YZ=-0.5)
  
  Naive.mat[,i] <- Sim1$Naive
  Naive2.mat[,i] <- Sim1$Naive2
  Onestep.mat[,i] <- Sim1$Onestep
  CI.Onestep.mat[,i] <- Sim1$Onestep.CI
}


write.table(Naive.mat,".../high_dim/Naive.txt")
write.table(Naive2.mat,".../high_dim/Naive2.txt")
write.table(Onestep.mat,".../high_dim/Onestep.txt")
write.table(CI.Onestep.mat,".../high_dim/Onestep_CI.txt")

