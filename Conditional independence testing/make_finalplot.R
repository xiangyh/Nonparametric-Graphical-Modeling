### Figure 2 ###

library(gridExtra)
library(cowplot)
library(ggplot2)

TypeI.plt <- function(plt, main="")
{
  plt + 
    #geom_hline(yintercept = 0.05,alpha=0.3) + 
    #geom_hline(yintercept = 0.01,alpha=0.3) + 
    geom_line(aes(linetype=Alpha, color=Test),size=0.8) + 
    geom_point(aes(shape = Test,color=Test),size=2.3) +
    facet_grid(rows = vars(Type),scales = "free_y") + 
    #facet_wrap(~Type, ncol=1, scales = "free_y") +
    xlab("Size") + ylab("") + 
    #expand_limits(y=ylt) +                        # Expand y range
    scale_y_continuous(breaks=seq(0,1,by=0.1)) +         # Set tick every 4
    theme_bw() +
    labs(title=main)+
    theme(panel.grid.minor = element_blank(), 
          legend.position="bottom")               # Position legend in bottom right
}

# Low dimensional case
load("../low_dim_type1/KCI_type1_1.rda")
load("../low_dim_type1/KCI_type1_5.rda")

load("../low_dim_type1/CDI_type1_1.rda")
load("../low_dim_type1/CDI_type1_5.rda")

load("../low_dim_type1/Ecov_type1_1.rda")
load("../low_dim_type1/Ecov_type1_5.rda")


size = c(50,100,200,300,400,500)


KCI_type1_1 <- data.frame(Test = "KCI-test",
                          Alpha = "0.01",
                          Type = "Type I error",
                          Size = size,
                          value = apply(KCI_type1_1,2,mean))


KCI_type1_5 <- data.frame(Test = "KCI-test",
                          Alpha = "0.05",
                          Type = "Type I error",
                          Size = size,
                          value = apply(KCI_type1_5,2,mean))



CDI_type1_1 <- data.frame(Test = "CDI-test",
                          Alpha = "0.01",
                          Type = "Type I error",
                          Size = size,
                          value = apply(CDI_type1_1,2,mean))


CDI_type1_5 <- data.frame(Test = "CDI-test",
                          Alpha = "0.05",
                          Type = "Type I error",
                          Size = size,
                          value = apply(CDI_type1_5,2,mean))


Ecov_type1_1 <- data.frame(Test = "scaled-Ecov",
                           Alpha = "0.01",
                           Type = "Type I error",
                           Size = size,
                           value = apply(Ecov_type1_1,2,mean))


Ecov_type1_5 <- data.frame(Test = "scaled-Ecov",
                           Alpha = "0.05",
                           Type = "Type I error",
                           Size = size,
                           value = apply(Ecov_type1_5,2,mean))

TypeI.low.dim <- rbind(Ecov_type1_1,Ecov_type1_5,KCI_type1_1,KCI_type1_5,CDI_type1_1,CDI_type1_5)



# power
load("../low_dim_type2/KCI_type2_1.rda")
load("../low_dim_type2/KCI_type2_5.rda")

load("../low_dim_type2/CDI_type2_1.rda")
load("../low_dim_type2/CDI_type2_5.rda")

load("../low_dim_type2/Ecov_type2_1.rda")
load("../low_dim_type2/Ecov_type2_5.rda")

KCI_power_1 <- data.frame(Test = "KCI-test",
                          Alpha = "0.01",
                          Type = "Power",
                          Size = size,
                          value = 1-apply(KCI_type2_1,2,mean))


KCI_power_5 <- data.frame(Test = "KCI-test",
                          Alpha = "0.05",
                          Type = "Power",
                          Size = size,
                          value = 1-apply(KCI_type2_5,2,mean))



CDI_power_1 <- data.frame(Test = "CDI-test",
                          Alpha = "0.01",
                          Type = "Power",
                          Size = size,
                          value = 1-apply(CDI_type2_1,2,mean))


CDI_power_5 <- data.frame(Test = "CDI-test",
                          Alpha = "0.05",
                          Type = "Power",
                          Size = size,
                          value = 1-apply(CDI_type2_5,2,mean))


Ecov_power_1 <- data.frame(Test = "scaled-Ecov",
                           Alpha = "0.01",
                           Type = "Power",
                           Size = size,
                           value = 1-apply(Ecov_type2_1,2,mean))


Ecov_power_5 <- data.frame(Test = "scaled-Ecov",
                           Alpha = "0.05",
                           Type = "Power",
                           Size = size,
                           value = 1-apply(Ecov_type2_5,2,mean))

TypeII.low.dim <- rbind(TypeI.low.dim,Ecov_power_1,Ecov_power_5,KCI_power_1,KCI_power_5,CDI_power_1,CDI_power_5)



# moderate dimensional case

load("../mod_dim_type1/KCI_type1_1.rda")
load("../mod_dim_type1/KCI_type1_5.rda")

load("../mod_dim_type1/CDI_type1_1.rda")
load("../mod_dim_type1/CDI_type1_5.rda")

load("../mod_dim_type1/Ecov_type1_1.rda")
load("../mod_dim_type1/Ecov_type1_5.rda")


KCI_type1_1 <- data.frame(Test = "KCI-test",
                          Alpha = "0.01",
                          Type = "Type I error",
                          Size = size,
                          value = apply(KCI_type1_1,2,mean))


KCI_type1_5 <- data.frame(Test = "KCI-test",
                          Alpha = "0.05",
                          Type = "Type I error",
                          Size = size,
                          value = apply(KCI_type1_5,2,mean))



CDI_type1_1 <- data.frame(Test = "CDI-test",
                          Alpha = "0.01",
                          Type = "Type I error",
                          Size = size,
                          value = apply(CDI_type1_1,2,mean))


CDI_type1_5 <- data.frame(Test = "CDI-test",
                          Alpha = "0.05",
                          Type = "Type I error",
                          Size = size,
                          value = apply(CDI_type1_5,2,mean))


Ecov_type1_1 <- data.frame(Test = "scaled-Ecov",
                           Alpha = "0.01",
                           Type = "Type I error",
                           Size = size,
                           value = apply(Ecov_type1_1,2,mean))


Ecov_type1_5 <- data.frame(Test = "scaled-Ecov",
                           Alpha = "0.05",
                           Type = "Type I error",
                           Size = size,
                           value = apply(Ecov_type1_5,2,mean))

TypeI.mod.dim <- rbind(Ecov_type1_1,Ecov_type1_5,KCI_type1_1,KCI_type1_5,CDI_type1_1,CDI_type1_5)



# power
load("../mod_dim_type2/KCI_type2_1.rda")
load("../mod_dim_type2/KCI_type2_5.rda")

load("../mod_dim_type2/CDI_type2_1.rda")
load("../mod_dim_type2/CDI_type2_5.rda")

load("../mod_dim_type2/Ecov_type2_1.rda")
load("../mod_dim_type2/Ecov_type2_5.rda")


KCI_power_1 <- data.frame(Test = "KCI-test",
                          Alpha = "0.01",
                          Type = "Power",
                          Size = size,
                          value = 1-apply(KCI_type2_1,2,mean))


KCI_power_5 <- data.frame(Test = "KCI-test",
                          Alpha = "0.05",
                          Type = "Power",
                          Size = size,
                          value = 1-apply(KCI_type2_5,2,mean))



CDI_power_1 <- data.frame(Test = "CDI-test",
                          Alpha = "0.01",
                          Type = "Power",
                          Size = size,
                          value = 1-apply(CDI_type2_1,2,mean))


CDI_power_5 <- data.frame(Test = "CDI-test",
                          Alpha = "0.05",
                          Type = "Power",
                          Size = size,
                          value = 1-apply(CDI_type2_5,2,mean))


Ecov_power_1 <- data.frame(Test = "scaled-Ecov",
                           Alpha = "0.01",
                           Type = "Power",
                           Size = size,
                           value = 1-apply(Ecov_type2_1,2,mean))


Ecov_power_5 <- data.frame(Test = "scaled-Ecov",
                           Alpha = "0.05",
                           Type = "Power",
                           Size = size,
                           value = 1-apply(Ecov_type2_5,2,mean))

TypeII.mod.dim <- rbind(TypeI.mod.dim,Ecov_power_1,Ecov_power_5,KCI_power_1,KCI_power_5,CDI_power_1,CDI_power_5)


low.dim <- TypeI.plt(ggplot(TypeII.low.dim, aes(x=Size, y=value))) + labs(title = "Low-dimensional setting",linetype=expression(alpha))
mod.dim <- TypeI.plt(ggplot(TypeII.mod.dim, aes(x=Size, y=value))) + labs(title = "Moderate-dimensional setting",linetype=expression(alpha))




## time
load("../CPU_time/KCI_time200.rda")
load("../CPU_time/KCI_time400.rda")
load("../CPU_time/CDI_time200.rda")
load("../CPU_time/CDI_time400.rda")
load("../CPU_time/Ecov_time200.rda")
load("../CPU_time/Ecov_time400.rda")

KCI_time200 <- data.frame(Test = "KCI_test",
                          Size = "200",
                          Dimension = 3:9,
                          value = apply(KCI_time200, 2, mean))

KCI_time400 <- data.frame(Test = "KCI_test",
                          Size = "400",
                          Dimension = 3:9,
                          value = apply(KCI_time400, 2, mean))



CDI_time200 <- data.frame(Test = "CDI_test",
                          Size = "200",
                          Dimension = 3:9,
                          value = apply(CDI_time200, 2, mean))

CDI_time400 <- data.frame(Test = "CDI_test",
                          Size = "400",
                          Dimension = 3:9,
                          value = apply(CDI_time400*60, 2, mean)) # turn minutes to seconds



Ecov_time200 <- data.frame(Test = "scaled-Ecov",
                          Size = "200",
                          Dimension = 3:9,
                          value = apply(Ecov_time200, 2, mean))

Ecov_time400 <- data.frame(Test = "scaled-Ecov",
                          Size = "400",
                          Dimension = 3:9,
                          value = apply(Ecov_time400, 2, mean))

Time <- rbind(Ecov_time200,Ecov_time400,KCI_time200,KCI_time400,CDI_time200,CDI_time400)


time.plt <- ggplot(Time, aes(x=Dimension, y=value)) +
  geom_line(aes(linetype=Size, color=Test),size=0.8) + 
  geom_point(aes(shape = Test,color=Test),size=2.3,alpha = 0.8) +
  xlab("Dimension") + ylab("Time (sec, log-scale)") + 
  #expand_limits(y=ylt) +                        # Expand y range
  scale_y_log10() +         # Set tick every 4
  theme_bw() +
  labs(title="CPU Time")+
  theme(legend.position="bottom")+
  guides(color=FALSE,shape=FALSE)


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend1<-g_legend(low.dim)
mylegend2<-g_legend(time.plt)


grid.arrange(arrangeGrob(low.dim+theme(legend.position="none"), 
                         mod.dim+theme(legend.position="none"), nrow = 1), 
             time.plt+theme(legend.position="none"),
             mylegend1, mylegend2, nrow=2,heights=c(12, 1), widths=c(2,1))
