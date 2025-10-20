library(latex2exp)
library(gridExtra)
theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )

# simulate data on Laplace scale
# select first two variables
Ylap <- as.data.frame(Yboot$Y[1:6000,1:2])
names(Ylap) <- c("Y1","Y2")
# including a threshold
v <- 0.95
vL <- quantile(Ylap$Y1,v)
tmp <- Ylap %>% mutate(above_thres= as.character(Y1>vL))
vL1 <- quantile(Ylap$Y2,v)
tmp1 <- Ylap %>% mutate(above_thres= as.character(Y2>vL1))
p2 <- ggplot(tmp) + geom_point(aes(x=Y1,y=Y2),size=0.5,alpha=0.5) +
         #           scale_color_manual(values = c("FALSE"="black","TRUE" = "#009ADA")) +
                    geom_vline(xintercept=vL,color="#009ADA",linetype="dashed") +
                    xlab(TeX("$Y_1$")) + ylab(TeX("$Y_2$"))+
                    theme(legend.position = "none") + geom_rect(xmax=20,xmin=vL,ymin=-20,ymax=20,colour="#009ADA",alpha=0.5)

p3 <-     p2 + geom_point(data=tmp1,aes(x=Y1,y=Y2),size=0.5,alpha=0.5) +
                 #   scale_color_manual(values = c("FALSE"="black","TRUE" = "#C11432")) +
                    geom_hline(yintercept=vL1,color="#C11432",linetype="dashed") +
                    xlab(TeX("$Y_1$")) + ylab(TeX("$Y_2$")) +
                    geom_rect(xmax=vL1,xmin=vL,ymin=-20,ymax=20,col="#C11432",alpha=0.5)



# series of plots

# p1 data on laplace scale
p1 <- ggplot(Ylap) + geom_point(aes(x=Y1,y=Y2),size=0.5) + xlab(TeX("$Y_1$")) + ylab(TeX("$Y_2$")) + scale_color_manual(values = c("black", "#009ADA")) + theme(legend.position="none") + coord_fixed() + xlim(-10,10) + ylim(-10,10) 

# p2 data with Y1 threshold

# p3 data with Y2 threshold

# p4 data with simulated body

# p5 data with simulated Y1

ggsave(p1,filename="../Documents/DataChallenge/p1.png")
ggsave(p2,filename="../Documents/DataChallenge/p2.png")
ggsave(p3,filename="../Documents/DataChallenge/p3.png")
ggsave(p4,filename="../Documents/DataChallenge/p4.png")
ggsave(p5,filename="../Documents/DataChallenge/p5.png")

