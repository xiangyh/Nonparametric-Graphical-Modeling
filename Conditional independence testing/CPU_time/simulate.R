# CPU

require(CondIndTests)
library(MASS)
library(bootstrap)
library(ranger)
library(cdcsis)
library(dplyr)

dim <- 3:9
iter <- 1

KCI_time200 <- matrix(NA, ncol = length(dim), nrow = iter)
KCI_time400 <- matrix(NA, ncol = length(dim), nrow = iter)

CDI_time200 <- matrix(NA, ncol = length(dim), nrow = iter)
CDI_time400 <- matrix(NA, ncol = length(dim), nrow = iter)

Ecov_time200 <- matrix(NA, ncol = length(dim), nrow = iter)
Ecov_time400 <- matrix(NA, ncol = length(dim), nrow = iter)

n <- 200
for (i in 1:iter) {
  for (d in 1:length(dim)) {
    p <- dim[d]
    z <- matrix(rnorm(n*p),n,p)
    beta <- 0.9^(1:p)
    x <- z%*%beta
    y <- z%*%exp(beta)
    
    KCI.t1 <- Sys.time()
    KCI_p <- KCI(y, x, z)
    KCI.t2 <- Sys.time()
    KCI_time200[i,d] <- difftime(KCI.t2, KCI.t1, units = "secs") %>% as.numeric() %>% round(3)
    
    CDI.t1 <- Sys.time()
    CDI_p <- cdcov.test(y, x, z)
    CDI.t2 <- Sys.time()
    CDI_time200[i,d] <- difftime(CDI.t2, CDI.t1, units = "secs") %>% as.numeric() %>% round(3)
    
    Ecov.t1 <- Sys.time()
    Hat.muY <- ranger(y~., data=as.data.frame(cbind(y,z)), num.trees = 300)$predictions
    Hat.muX <- ranger(x~., data=as.data.frame(cbind(x,z)), num.trees = 300)$predictions
    ECov2 <-  mean((y-Hat.muY)*(x-Hat.muX))
    EVarY2 <- mean((y-Hat.muY)^2)
    EVarX2 <- mean((x-Hat.muX)^2)
    Dp = (y-Hat.muY)*(x-Hat.muX)/sqrt(EVarY2*EVarX2) - ECov2/(2*sqrt(EVarY2*EVarX2))*((y-Hat.muY)^2/EVarY2 + (x-Hat.muX)^2/EVarX2)
    s2 <- mean((Dp)^2)
    onestep <- ECov2/sqrt(EVarY2*EVarX2)
    test_stat <- sqrt(n)*onestep/sqrt(s2)
    Ecov_p <- 2*pnorm(-abs(test_stat))
    Ecov.t2 <- Sys.time()
    Ecov_time200[i,d] <- difftime(Ecov.t2, Ecov.t1, units = "secs") %>% as.numeric() %>% round(3)
  }
}


n <- 400
for (i in 1:iter) {
  for (d in 1:length(dim)) {
    p <- dim[d]
    z <- matrix(rnorm(n*p),n,p)
    beta <- 0.9^(1:p)
    x <- z%*%beta
    y <- z%*%exp(beta)
    
    KCI.t1 <- Sys.time()
    KCI_p <- KCI(y, x, z)
    KCI.t2 <- Sys.time()
    KCI_time400[i,d] <- difftime(KCI.t2, KCI.t1, units = "secs") %>% as.numeric() %>% round(3)
    
    CDI.t1 <- Sys.time()
    CDI_p <- cdcov.test(y, x, z)
    CDI.t2 <- Sys.time()
    CDI_time400[i,d] <- difftime(CDI.t2, CDI.t1, units = "secs") %>% as.numeric() %>% round(3)
    
    Ecov.t1 <- Sys.time()
    Hat.muY <- ranger(y~., data=as.data.frame(cbind(y,z)), num.trees = 300)$predictions
    Hat.muX <- ranger(x~., data=as.data.frame(cbind(x,z)), num.trees = 300)$predictions
    ECov2 <-  mean((y-Hat.muY)*(x-Hat.muX))
    EVarY2 <- mean((y-Hat.muY)^2)
    EVarX2 <- mean((x-Hat.muX)^2)
    Dp = (y-Hat.muY)*(x-Hat.muX)/sqrt(EVarY2*EVarX2) - ECov2/(2*sqrt(EVarY2*EVarX2))*((y-Hat.muY)^2/EVarY2 + (x-Hat.muX)^2/EVarX2)
    s2 <- mean((Dp)^2)
    onestep <- ECov2/sqrt(EVarY2*EVarX2)
    test_stat <- sqrt(n)*onestep/sqrt(s2)
    Ecov_p <- 2*pnorm(-abs(test_stat))
    Ecov.t2 <- Sys.time()
    Ecov_time400[i,d] <- difftime(Ecov.t2, Ecov.t1, units = "secs") %>% as.numeric() %>% round(3)
  }
}




save(KCI_time200,file=".../CPU_time/KCI_time200.rda")
save(KCI_time400,file=".../CPU_time/KCI_time400.rda")
save(CDI_time200,file=".../CPU_time//CDI_time200.rda")
save(CDI_time400,file=".../CPU_time//CDI_time400.rda")
save(Ecov_time200,file=".../CPU_time/Ecov_time200.rda")
save(Ecov_time400,file=".../CPU_time/Ecov_time400.rda")
