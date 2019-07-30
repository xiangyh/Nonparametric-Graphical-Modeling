# simulation: conditional independent testing
require(CondIndTests)
library(MASS)
library(bootstrap)
library(ranger)
library(cdcsis)
library(bootstrap)
library(dplyr)


Iter = 200
size = c(50,100,200,300,400,500)


# Type II error
KCI_type2_1 <- matrix(NA, ncol = length(size), nrow = Iter)
KCI_type2_5 <- matrix(NA, ncol = length(size), nrow = Iter)

CDI_type2_1 <- matrix(NA, ncol = length(size), nrow = Iter)
CDI_type2_5 <- matrix(NA, ncol = length(size), nrow = Iter)

Ecov_type2_1 <- matrix(NA, ncol = length(size), nrow = Iter)
Ecov_type2_5 <- matrix(NA, ncol = length(size), nrow = Iter)


for (i in 1:Iter) {
  for (j in 1:length(size)) {
    n <- size[j]
    
    # simulate data
    z <- rnorm(n)
    zz1 <- 0.7 *(z^3/5+z/2)
    zz2 = (z^3/4 + z)/3
    z2 <- rnorm(n)
    zz1_2 <- zz1/2 + z2
    zz1_2 <- zz1_2/2 + 0.7 * tanh(zz1_2)
    zz2_2 <- zz2/2 +z2
    zz2_2 <- zz2_2/2 + 0.7*tanh(zz2_2)
    z3 <- rnorm(n)
    zz1_3 <- zz1_2*2/3 + z3*5/6
    zz1_3 = zz1_3/2 + .7 * tanh(zz1_3)
    zz2_3 = zz2_2*2/3 + z3*5/6
    zz2_3 = zz2_3/2 + .7 * tanh(zz2_3)
    z4 <- rnorm(n)
    zz1_4 = zz1_3*2/3 + z4*5/6
    zz1_4 = zz1_4/2 + .7 * tanh(zz1_4)
    zz2_4 = zz2_3*2/3 + z4*5/6
    zz2_4 = zz2_4/2 + .7 * tanh(zz2_4)
    z5 <- rnorm(n)
    zz1_5 = zz1_4*2/3 + z5*5/6
    zz1_5 = zz1_5/2 + .7 * tanh(zz1_5)
    zz2_5 = zz2_4*2/3 + z5*5/6
    zz2_5 = zz2_5/2 + .7 * tanh(zz2_5)
    z = cbind(z,z2,z3,z4,z5)
    x <- rnorm(n)
    y <- rnorm(n)
    ff = rnorm(n,mean=1,sd=2)
    x = zz1_5 + tanh(x) + (x^3)/3 + tanh(x/3)/2 + cosh(ff)
    y = y + zz2_5 + tanh(y/3) + cosh(ff^2)
    
    
    # KCI
    KCI_p <- KCI(y, x, z)$pvalue
    KCI_type2_1[i,j] <- KCI_p > 0.01
    KCI_type2_5[i,j] <- KCI_p > 0.05
    
    # conditional independence testing
    CDI_p <- cdcov.test(y, x, z)$p.value
    CDI_type2_1[i,j] <-  CDI_p > 0.01
    CDI_type2_5[i,j] <- CDI_p > 0.05
    
    # Ecov
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
    
    Ecov_type2_1[i,j] <- Ecov_p > 0.01
    Ecov_type2_5[i,j] <- Ecov_p > 0.05
    
  }
}

save(KCI_type2_1,file=".../mod_dim_type2/KCI_type2_1.rda")
save(KCI_type2_5,file=".../mod_dim_type2/KCI_type2_5.rda")
save(CDI_type2_1,file=".../mod_dim_type2/CDI_type2_1.rda")
save(CDI_type2_5,file=".../mod_dim_type2/CDI_type2_5.rda")
save(Ecov_type2_1,file=".../mod_dim_type2/Ecov_type2_1.rda")
save(Ecov_type2_5,file=".../mod_dim_type2/Ecov_type2_5.rda")




