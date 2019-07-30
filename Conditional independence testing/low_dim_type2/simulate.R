# simulation: conditional independent testing
require(CondIndTests)
library(MASS)
library(bootstrap)
library(ranger)
library(cdcsis)
library(bootstrap)

loess_wrapper_extrapolate <- function (x, y, span.vals = seq(0.08, 0.8, by = 0.05), folds = 5){
  # Do model selection using mean squared error
  mean.sqr.error <- numeric(length(span.vals))
  
  # Quantify error for each span, using CV
  loess.model <- function(x, y, span){
    loess(y ~ x, span = span)
  }
  
  loess.predict <- function(fit, newdata) {
    predict(fit, newdata = newdata)
  }
  
  span.index <- 0
  for (each.span in span.vals) {
    span.index <- span.index + 1
    y.hat.cv <- crossval(x, y, theta.fit = loess.model, theta.predict = loess.predict, span = each.span, ngroup = folds)$cv.fit
    non.empty.indices <- !is.na(y.hat.cv)
    mean.sqr.error[span.index] <- mean((y[non.empty.indices] - y.hat.cv[non.empty.indices])^2)
  }
  
  # find the span which minimizes error
  best.span <- span.vals[which.min(mean.sqr.error)]
  
  # fit and return the best model
  best.model <- loess(y ~ x, span = best.span, control=loess.control(surface="direct"))
  return(list(model = best.model))
}

True_Y <- function(z){
  y <- sin(3*z)
  return(y)
}
True_X <- function(z){
  x <- cos(2*z)
  return(x)
}

e_mu <- c(0,0)
e_sigma <- matrix(c(1,-0.5,-0.5,1),2,2)
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
    z <- runif(n, 0, 2)
    e <- mvrnorm(n, e_mu, e_sigma)
    mu_X <- True_X(z); x <- mu_X + e[,1]
    mu_Y <- True_Y(z); y <- mu_Y + e[,2]
    
    # KCI
    KCI_p <- KCI(y, x, z)$pvalue
    KCI_type2_1[i,j] <- KCI_p > 0.01
    KCI_type2_5[i,j] <- KCI_p > 0.05
    
    # conditional independence testing
    CDI_p <- cdcov.test(y, x, z)$p.value
    CDI_type2_1[i,j] <-  CDI_p > 0.01
    CDI_type2_5[i,j] <- CDI_p > 0.05
    
    # Ecov
    Hat.muY <- loess_wrapper_extrapolate(z, y)$model$fitted
    Hat.muX <- loess_wrapper_extrapolate(z, x)$model$fitted
    
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


save(KCI_type2_1,file=".../low_dim_type2/KCI_type2_1.rda")
save(KCI_type2_5,file=".../low_dim_type2/KCI_type2_5.rda")
save(CDI_type2_1,file=".../low_dim_type2/CDI_type2_1.rda")
save(CDI_type2_5,file=".../low_dim_type2/CDI_type2_5.rda")
save(Ecov_type2_1,file=".../low_dim_type2/Ecov_type2_1.rda")
save(Ecov_type2_5,file=".../low_dim_type2/Ecov_type2_5.rda")




