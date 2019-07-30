require(CondIndTests)
library(MASS)
library(bootstrap)
library(ranger)
library(cdcsis)
library(bootstrap)
library(dplyr)


Iter = 200
size = c(50,100,200,300,400,500)


# Type I error and power:
KCI_type1_1 <- matrix(NA, ncol = length(size), nrow = Iter)
KCI_type1_5 <- matrix(NA, ncol = length(size), nrow = Iter)

CDI_type1_1 <- matrix(NA, ncol = length(size), nrow = Iter)
CDI_type1_5 <- matrix(NA, ncol = length(size), nrow = Iter)

Ecov_type1_1 <- matrix(NA, ncol = length(size), nrow = Iter)
Ecov_type1_5 <- matrix(NA, ncol = length(size), nrow = Iter)

for (i in 1:Iter) {
  for (j in 1:length(size)) {
    n <- size[j]
    
    # simulate data
    x <- rnorm(n)
    y <- rnorm(n)
    z <- rnorm(n)
    zz1 <- 0.7 *(z^3/5+z/2)
    zz2 = (z^3/4 + z)/3
    z2 <- rnorm(n)
    zz1_2 <- zz1/2 + z2
    zz1_2 <- zz1_2/2 + 0.7 * tanh(zz1_2)
    zz2_2 <- zz2/2 +z2
    zz2_2 <- zz2_2/2 + 0.7 * tanh(zz2_2)
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
    
    x = zz1_5 + tanh(x) + (x^3)/3 + tanh(x/3)/2
    y = y + zz2_5 + tanh(y/3)
    
    z = cbind(z,z2,z3,z4,z5)
    
    # KCI
    KCI_p <- KCI(y, x, z)$pvalue
    KCI_type1_1[i,j] <- KCI_p <= 0.01
    KCI_type1_5[i,j] <- KCI_p <= 0.05
    
    # conditional independence testing
    CDI_p <- cdcov.test(y, x, z)$p.value
    CDI_type1_1[i,j] <-  CDI_p <= 0.01
    CDI_type1_5[i,j] <- CDI_p <= 0.05
    
    # Ecov
    Hat.muY <- ranger(y~., data=as.data.frame(cbind(y,z)), num.trees = 400)$predictions
    Hat.muX <- ranger(x~., data=as.data.frame(cbind(x,z)), num.trees = 400)$predictions
    
    Dp = (y-Hat.muY)*(x - Hat.muX)
    s2 <- mean((Dp-mean(Dp))^2)
    onestep <- mean(Dp)
    test_stat <- sqrt(n)*onestep/sqrt(s2)
    Ecov_p <- 2*pnorm(-abs(test_stat))
    
    Ecov_type1_1[i,j] <- Ecov_p <= 0.01
    Ecov_type1_5[i,j] <- Ecov_p <= 0.05
    
  }
}

save(KCI_type1_1,file=".../mod_dim_type1/KCI_type1_1.rda")
save(KCI_type1_5,file=".../mod_dim_type1/KCI_type1_5.rda")
save(CDI_type1_1,file=".../mod_dim_type1/CDI_type1_1.rda")
save(CDI_type1_5,file=".../mod_dim_type1/CDI_type1_5.rda")
save(Ecov_type1_1,file=".../mod_dim_type1/Ecov_type1_1.rda")
save(Ecov_type1_5,file=".../mod_dim_type1/Ecov_type1_5.rda")




