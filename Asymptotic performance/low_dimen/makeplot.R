library(gridExtra)
library(cowplot)
library(ggplot2)
library(grid)


Naive <- read.table("/Users/yunhuaxiang/Desktop/independent study/Simulation/Covariance/7_10/Sim1/Naive1.txt")
Onestep <- read.table("/Users/yunhuaxiang/Desktop/independent study/Simulation/Covariance/7_10/Sim1/Onestep.txt")
Onestep_CI <- read.table("/Users/yunhuaxiang/Desktop/independent study/Simulation/Covariance/7_10/Sim1/Onestep1_CI.txt")
Naive_CI <- read.table("/Users/yunhuaxiang/Desktop/independent study/Simulation/Covariance/7_10/Sim1/Naive1_CI.txt")


Cover.ggplot <- function(plt, main="")
{
  plt + 
    geom_hline(yintercept = 0.95, linetype=2) + 
    geom_point(aes(color = Est), position=position_dodge(100), size=2.5) + # 21 is filled circle
    scale_color_manual(values = c("red","blue","green"))+
    xlab("Size") +
    ggtitle("The Effect of Vitamin C on\nTooth Growth in Guinea Pigs") +
    #expand_limits(y=0) +                        # Expand y range
    scale_y_continuous(breaks=seq(0,1,by=0.2)) +         # Set tick every 4
    theme_bw() +
    labs(title=main)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position="bottom")               # Position legend in bottom right
}



MSE.ggplot <- function(plt, main="")
{
  plt + 
    geom_errorbar(aes(ymin=value.l, ymax=value.u), colour="black", width=200, position=pd, size=0.2) +
    #geom_line(position=pd) +
    geom_hline(yintercept = 0, linetype=2) + 
    geom_point(aes(color = Est), position=pd, size=2.5) + # 21 is filled circle
    scale_color_manual(values = c("red","blue","green"))+
    xlab("Size") +
    #expand_limits(y=0) +                        # Expand y range
    #scale_y_continuous(breaks=seq(0,1,by=0.2)) +         # Set tick every 4
    theme_bw() +
    labs(title=main)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position="bottom")               # Position legend in bottom right
}

size <- c(100, 500, 1000, 2000,3000, 6000)

onestep_ci <- data.frame(Est = "Plug-in", size=size, CI = apply(Onestep_CI, 2, mean))
naive_ci <- data.frame(Est = "Naive", size=size, CI = apply(Naive_CI, 2, mean))

Onestep_ci <- rbind(naive_ci,onestep_ci)

pd <- position_dodge(100) # move them .05 to the left and right
ci.plt <- ggplot(Onestep_ci, aes(x=size, y=CI, group=Est))
cover.plt <- Cover.ggplot(ci.plt)+ylab("Coverage")


true_Cov_YZ=-0.5
N <- dim(Naive)[1]
# 
# Naive.bias <- data.frame(Est = "Naive", size = size, 
#                           value = sqrt(size)*(apply(Naive,2,mean)-true_Cov_YZ),
#                           value.l = sqrt(size)*(apply(Naive,2,mean) -true_Cov_YZ - apply(Naive, 2, sd)/sqrt(N)*1.96),
#                           value.u = sqrt(size)*(apply(Naive,2,mean) -true_Cov_YZ + apply(Naive, 2, sd)/sqrt(N)*1.96)
#                           )
# 
# Naive2.bias <- data.frame(Est = "Naive2", size = size, 
#                           value = sqrt(size)*(apply(Naive2,2,mean)-true_Cov_YZ),
#                           value.l = sqrt(size)*(apply(Naive2,2,mean) -true_Cov_YZ - apply(Naive2, 2, sd)/sqrt(N)*1.96),
#                           value.u = sqrt(size)*(apply(Naive2,2,mean) -true_Cov_YZ + apply(Naive2, 2, sd)/sqrt(N)*1.96)
# )
# 
# Onestep.bias <- data.frame(Est = "Onestep", size = size, 
#                            value = sqrt(size)*(apply(Onestep,2,mean)-true_Cov_YZ),
#                            value.l = sqrt(size)*(apply(Onestep,2,mean) -true_Cov_YZ - apply(Onestep, 2, sd)/sqrt(N)*1.96),
#                            value.u = sqrt(size)*(apply(Onestep,2,mean) -true_Cov_YZ + apply(Onestep, 2, sd)/sqrt(N)*1.96)
# )
# 
# Bias <- rbind(Naive.bias, Naive2.bias, Onestep.bias)
# bias.plt <- ggplot(Bias, aes(x=size, y=value, group=Est))
# Bias.plt <- MSE.ggplot(bias.plt)+ylab(expression(sqrt(n)~~x~~Bias))


## bootstrap
iter = 400
Boot = 1000
bias.Naive.bs <- matrix(ncol = length(size), nrow = Boot)
bias.Onestep.bs <- matrix(ncol = length(size), nrow = Boot)

var.Naive.bs <-  matrix(ncol = length(size), nrow = Boot)
var.Onestep.bs <- matrix(ncol = length(size), nrow = Boot)


for (b in 1:Boot) {
  idx <- sample(1:iter, iter, replace=T)
  Naive.bs <- Naive[idx,]
  Onestep.bs <- Onestep[idx,]
  
  # bias
  bias.Naive.bs[b,] <- apply(Naive.bs, 2, mean)-true_Cov_YZ
  bias.Onestep.bs[b,] <- apply(Onestep.bs, 2, mean)-true_Cov_YZ
  
  # mse
  var.Naive.bs[b,] <- apply(Naive.bs, 2, var)
  var.Onestep.bs[b,] <- apply(Onestep.bs, 2, var)
}

# bias

Naive.bias <- data.frame(Est = "Naive", size = size, 
                          value = sqrt(size)*(apply(Naive,2,mean)-true_Cov_YZ),
                          value.u = sqrt(size)*apply(bias.Naive.bs, 2, function(x) quantile(x, 0.975)),
                          value.l = sqrt(size)*apply(bias.Naive.bs, 2, function(x) quantile(x, 0.025)))


Onestep.bias <- data.frame(Est = "Plug-in", size = size, 
                           value = sqrt(size)*(apply(Onestep,2,mean)-true_Cov_YZ),
                           value.u = sqrt(size)*apply(bias.Onestep.bs, 2, function(x) quantile(x, 0.975)),
                           value.l = sqrt(size)*apply(bias.Onestep.bs, 2, function(x) quantile(x, 0.025)))

Bias <- rbind(Naive.bias, Onestep.bias)
bias.plt <- ggplot(Bias, aes(x=size, y=value, group=Est))
Bias.plt <- MSE.ggplot(bias.plt)+ylab(expression(sqrt(n)~~x~~Bias))

# var
Naive.var <- data.frame(Est = "Naive", size = size, 
                         value = size*apply(Naive,2,var),
                         value.u = size*apply(var.Naive.bs, 2, function(x) quantile(x, 0.975)),
                         value.l = size*apply(var.Naive.bs, 2, function(x) quantile(x, 0.025)))


Onestep.var <- data.frame(Est = "Plug-in", size = size, 
                          value = size*apply(Onestep,2,var),
                          value.u = size*apply(var.Onestep.bs, 2, function(x) quantile(x, 0.975)),
                          value.l = size*apply(var.Onestep.bs, 2, function(x) quantile(x, 0.025)))


# Naive.var <- data.frame(Est = "Naive", size = size, 
#                          value = size*apply(Naive,2,var))
# 
# Naive2.var <- data.frame(Est = "Naive2", size = size, 
#                          value = size*apply(Naive2,2,var))
# 
# Onestep.var <- data.frame(Est = "Onestep", size = size, 
#                           value = size*apply(Onestep,2,var))


Var <- rbind(Naive.var, Onestep.var)
var.plt <- ggplot(Var, aes(x=size, y=value, group=Est))
Var.plt <- MSE.ggplot(var.plt)+ylab(expression(n~~x~~Variance))



# Naive.mse <- data.frame(Est = "Naive", size = size, 
#                          value = sqrt(size)*apply((Naive-true_Cov_YZ)^2,2,mean))
# 
# Naive2.mse <- data.frame(Est = "Naive2", size = size, 
#                          value = sqrt(size)*apply((Naive2-true_Cov_YZ)^2,2,mean))
# 
# Onestep.mse <- data.frame(Est = "Onestep", size = size, 
#                           value = sqrt(size)*apply((Onestep-true_Cov_YZ)^2,2,mean))
# 
# 
# Mse <- rbind(Naive.mse, Naive2.mse, Onestep.mse)
# mse.plt <- ggplot(Mse, aes(x=size, y=value, group=Est))
# Mse.plt <- MSE.ggplot(mse.plt)+ylab("MSE")
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(Var.plt)

pdf("/Users/yunhuaxiang/Desktop/independent study/Simulation/Covariance/7_10/plot3.pdf", height = 3, width = 5.5)
grid.arrange(arrangeGrob(Bias.plt+theme(legend.position="none"), 
                         cover.plt+theme(legend.position="none"), nrow = 1), 
             mylegend, nrow=2,heights=c(10, 1),
             top = textGrob("Low-dimensional setting",gp=gpar(fontsize=15), vjust = 1.3))
dev.off()
 
