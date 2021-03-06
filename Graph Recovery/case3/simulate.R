########################################
# Graph recovery simulation for Case 3 #
# transelliptical model                #
########################################

library(MASS)
library(copula)
library(dplyr)
library(cdcsis)
library(CondIndTests)
library(clime)
library(ROCR)
library(ranger)
library(psych)  # for pairs.panels()
library(huge)

# 
source(".../clime_kendall.R")

create_adj_cov <- function(n,p,r)
{
  
  L = huge.generator(n=n, d=p, u=0.7,v=0.3, prob = 0.5)
  adj_mat <- L$theta %>% as.matrix()
  cov_mat <- L$sigma
  
  Xi <- rchisq(n,df=p)^0.5
  ev <- eigen(cov_mat)
  V <- ev$vectors
  L <- ev$values[!is.na(ev$values)]
  A <- V%*% diag(sqrt(L))
  U <- matrix(rnorm(length(L)*n,sd=1), nrow = length(L), ncol = n)
  U <- apply(U, 2, function(x) x/sqrt(sum(x^2)))
  X_mat0 <- X_mat <- diag(Xi) %*% t(A%*%U)
  X_mat[,1] <- exp(X_mat0[,1])^0.5
  X_mat[,2] <- sign(X_mat0[,2])*abs(X_mat0[,2])^0.5
  X_mat[,3] <- X_mat0[,3]^3
  X_mat[,4] <- pnorm(X_mat0[,4])
  X_mat[,5] <- exp(X_mat0[,5])^0.5
  X_mat[,6] <- sign(X_mat0[,6])*abs(X_mat0[,6])^0.5
  X_mat[,7] <- X_mat0[,7]^3
  X_mat[,8] <- pnorm(X_mat0[,8])
  
  cov_inv_est_ggm <- cov(X_mat) %>% solve()
  re.clime <- clime2(X_mat,nlambda = 10)
  re.cv <- cv.clime2(re.clime)
  cov_inv_est_clime <- clime2(X_mat, re.cv$lambdaopt)$Omegalist[[1]]
  cov_inv_est_ecov <- matrix(rep(0,p^2),nrow = p)

  
  cb <- combn(p,2)
  for (i in 1:dim(cb)[2]) {
    idx1 <- cb[,i][1]
    idx2 <- cb[,i][2]
    Y <- X_mat[,idx1] 
    Z <- X_mat[,idx2]
    X <- X_mat[,-cb[,i]] %>% as.data.frame()
    Hat.muY <- ranger(Y~., data=cbind(Y,X), num.trees = 500)$predictions
    Hat.muZ <- ranger(Z~., data=cbind(Z,X), num.trees = 500)$predictions
    
    ECov2 <-  mean((Y-Hat.muY)*(Z-Hat.muZ))
    EVarY2 <- mean((Y-Hat.muY)^2)
    EVarZ2 <- mean((Z-Hat.muZ)^2)
    onestep <- ECov2/sqrt(EVarY2*EVarZ2)
    cov_inv_est_ecov[idx1,idx2] <- cov_inv_est_ecov[idx2,idx1] <- onestep
  }
  diag(adj_mat) <- NA
  diag(cov_inv_est_ggm) <- NA
  diag(cov_inv_est_ecov) <- NA
  diag(cov_inv_est_clime) <- NA
  
  adj_mat <- adj_mat[-seq(1,p^2,p+1)]
  cov_inv_est_ggm <- cov_inv_est_ggm[-seq(1,p^2,p+1)]
  cov_inv_est_ecov <- cov_inv_est_ecov[-seq(1,p^2,p+1)]
  cov_inv_est_clime <- cov_inv_est_clime[-seq(1,p^2,p+1)]
  
  return(list(adj_mat = adj_mat,
              cov_inv_est_ggm=cov_inv_est_ggm, 
              cov_inv_est_ecov=cov_inv_est_ecov, 
              cov_inv_est_clime=cov_inv_est_clime))
}

n=400
p=8
r=0.145
iter=20
Adj_Mat <- list()
Cov_inv_Est_GMM <- list()
Cov_inv_Est_Ecov <- list()
Cov_inv_Est_CLIME <- list()

for (k in 1:iter) {
  est <- create_adj_cov(n,p,r)
  Cov_inv_Est_GMM[[k]] <- abs(est$cov_inv_est_ggm)
  Cov_inv_Est_Ecov[[k]] <- abs(est$cov_inv_est_ecov)
  Cov_inv_Est_CLIME[[k]] <- abs(est$cov_inv_est_clime)
  Adj_Mat[[k]] <- est$adj_mat
}


pred <- prediction(Cov_inv_Est_Ecov, Adj_Mat)
perf <- performance(pred,"tpr","fpr")
pred2 <- prediction(Cov_inv_Est_GMM, Adj_Mat)
perf2 <- performance(pred2,"tpr","fpr")
pred3 <- prediction(Cov_inv_Est_CLIME, Adj_Mat)
perf3 <- performance(pred3,"tpr","fpr")


par(mfrow=c(1,1),mar=c(3.5, 3.5, 1.7, 0.8), mgp=c(2, 0.8, 0))
plot(perf,col="red",lty=1, lwd=2, avg="threshold", main="Case 3", xlab="False positive rate",ylab="True positive rate")
plot(perf2,lty=2, lwd=2, avg="threshold",add=TRUE)
plot(perf3,lty=2, lwd=2, add=TRUE, avg="threshold",col="green")
#legend("bottomright", c("ECov", "GGM","TGM"), col = c("red", 1,"green"), lty = c(1,2,2))
