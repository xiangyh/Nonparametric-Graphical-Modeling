##################################################
# Graph recovery simulation for Case 4           #
# Non-Gaussian, Non-Copulas, non-transelliptical #
##################################################

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
library(igraph)


# 
source("/home/yxiang/graph_recovery/clime_kendall.R")

dat <- read.table(text="A B 
1 2
                  1 3
                  1 4
                  1 5
                  1 6
                  1 7
                  1 8
                  1 9
                  1 10
                  6 7
                  6 8
                  6 9
                  6 10", header=TRUE)



adj_mat <- get.adjacency(graph.edgelist(as.matrix(dat), directed=FALSE)) %>% as.matrix()

create_adj_cov <- function(n,p,adj_mat, r=0.2)
{
  X_mat <- matrix(NA,n,p)
  X_mat[,2:p] <- matrix(rexp(n*(p-1), rate=2),n,p-1)
  X_mat[,6] <- sin(X_mat[,7]) + X_mat[,8]^2 + X_mat[,9] + X_mat[,10]^2
  X_mat[,1] <- X_mat[,2]^2+X_mat[,3]+2*X_mat[,4]+sin(X_mat[,5])+X_mat[,6]+exp(X_mat[,7]) + X_mat[,8] + X_mat[,9] + X_mat[,10]
  
  cov_inv_est_ggm <- cov(X_mat) %>% solve()
  re.clime <- clime2(X_mat,nlambda = 10)
  re.cv <- cv.clime2(re.clime)
  cov_inv_est_clime <- clime2(X_mat, re.cv$lambdaopt)$Omegalist[[1]]
  cov_inv_est_ecov <- matrix(rep(0,p^2),nrow = p)
  
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

    Dp <- (Y - Hat.muY)*(Z - Hat.muZ)
    onestep <- mean(Dp)
    sd <- sqrt(mean((Dp-mean(Dp))^2))/sqrt(n)
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
  
  return(list(adj_mat=adj_mat, cov_inv_est_ggm=cov_inv_est_ggm, cov_inv_est_ecov=cov_inv_est_ecov,cov_inv_est_clime=cov_inv_est_clime))
}



n=400
p=10
iter=10
Adj_Mat <- list()
Cov_inv_Est_GMM <- list()
Cov_inv_Est_Ecov <- list()
Cov_inv_Est_CLIME <- list()

for (k in 1:iter) {
  est <- create_adj_cov(n,p,adj_mat)
  
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
plot(perf,col="red",lty=1, lwd=2, avg="threshold", main="Case 4", xlab="False positive rate",ylab="True positive rate")
plot(perf2,lty=2, lwd=2, avg="threshold",add=TRUE)
plot(perf3,lty=2, lwd=2, add=TRUE, avg="threshold",col="green")
