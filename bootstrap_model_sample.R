# libraries ---------
library(parallel)
library(LaplacesDemon)
library(dplyr)
library(texmex)
library(evd)
library(rvinecopulib)
library(tidyverse)
library(VGAM)
library(gtools)
library(purrr)
library(quantreg)

# function scripts ------
source('conditional_model_helpers.R')
source('conditional_model_fit.R')
source('target_quantity_helpers.R')
source('model_diagnostics.R')
source("NonSta_GPD_Functions.R")
source("NonSta_GPD_Fitting.R")

# load data -----
load("runs/run1.RData")
load("runs/run2.RData")
load("runs/run3.RData")
load("runs/run4.RData")


# obtain bootstrap and refit model -----
run_number <- 1
  source("Block_Bootstrapping.R")
    start_time <- Sys.time()
    Y_Bootstrapped <- NonSta_GPD_to_Lapalce(df = data,
                                            run = run_number,
                                            tau = 0.98,
                                            B = 1,
                                            l = 2)
    end_time <- Sys.time()
    end_time - start_time
    Yboot <- Y_Bootstrapped[[1]]
    names(Yboot$Y) <- paste0("Y",1:25)
    Yboot$Y <- as.data.frame(Yboot$Y)
    q = 0.95 # quantile threshold for fitting the conditional model
    cond_model_fit_wrapper = function(i){
      return(fn(site_index=i,v=q,Yrun1Lap =Yboot$Y,res_dist = "AGG_vinecopula"))
    }
    cond_modelvc1 <- sapply(1:25,cond_model_fit_wrapper,simplify=FALSE)
    
Nrun <- 1
sims_vc <- sim_cond_model(Yrun=Yboot$Y,q=q,cond_model=cond_modelvc1,res_dist = "AGG_vinecopula",Y_boot=Yboot)
sims_emp <- sim_cond_model(Yrun=Yboot$Y,q=q,cond_model=cond_modelvc1,res_dist = "empirical",Y_boot=Yboot)
    
save(Yboot,cond_modelvc1,sims_vc,sims_emp,file="../run1_bootstrap_resfit.RData")

# load to also save on other margins
load("../run1_bootstrap_resfit.RData")
set.seed(1234)
sims_vc <- sim_cond_model(seed=1234,Yrun=Yboot$Y,cond_model=cond_modelvc1,res_dist = "AGG_vinecopula",Y_boot=Yboot,include_temporal=TRUE,return_Laplace=FALSE)
# exclude temporal reordering and run again
set.seed(1234)
sims_notemp <- sim_cond_model(seed=1234,Yrun=Yboot$Y,cond_model=cond_modelvc1,res_dist = "AGG_vinecopula",Y_boot=Yboot,include_temporal = FALSE,return_Laplace = FALSE)

# return Laplace
set.seed(1234)
sims_vc_Laplace <-  sim_cond_model(seed=1234,Yrun=Yboot$Y,cond_model=cond_modelvc1,res_dist = "AGG_vinecopula",Y_boot=Yboot, include_temporal = TRUE, return_Laplace = TRUE)
# comment out temporal reordering and return Laplace
set.seed(1234)
sims_notemp_Laplace <-  sim_cond_model(seed=1234,Yrun=Yboot$Y,cond_model=cond_modelvc1,res_dist = "AGG_vinecopula",Y_boot=Yboot, include_temporal = FALSE, return_Laplace = TRUE)
# save
save(Yboot,cond_modelvc1,sims_vc,sims_notemp,sims_vc_Laplace,sims_notemp_Laplace,file="../run1_bootstrap_resfit.RData")

#compare vine copula fit and Gaussian copula fit
# source Aiden's functions
source("MVAGG_Functions_2.R")
q = 0.95 # quantile threshold for fitting the conditional model
cond_model_fit_wrapper = function(i){
  return(fn(site_index=i,v=q,Yrun1Lap =Yboot$Y,res_dist = "AGG_Gauscopula"))
}
cond_model_Gauscop1 <- sapply(1:25,cond_model_fit_wrapper,simplify=FALSE)

set.seed(1234)
sims_gauscop <- sim_cond_model(seed=1234,Yrun=Yboot$Y,cond_model=cond_model_Gauscop1,res_dist = "AGG_Gauscopula",Y_boot=Yboot,include_temporal=TRUE,return_Laplace=FALSE)

theme_set(theme_bw())
theme_replace(
  panel.spacing = unit(2, "lines"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA) )


# subset where two variables are extreme
j <- 1
nv <- dim(sims_vc)[1]*(1-q)
x <- sort(as.data.frame(Yboot$X) %>% pull(all_of(1)),decreasing=TRUE)[1:nv]
y1 <- sort(sims_vc %>% pull(all_of(1)),decreasing=TRUE)[1:nv]
y2 <- sort(sims_gauscop %>% pull(all_of(1)),decreasing=TRUE)[1:nv]
tmp <- rbind(data.frame("data"=x,"model" = y1,"Z_dependence"="vine_copula"),data.frame("data"=x,"model" = y2,"Z_dependence"="gaussian_copula"))
ggplot(tmp) + geom_point(aes(x=data,y=model,col=Z_dependence)) + geom_abline(slope=1,linetype="dashed")


# try overlapping
q <- 0.99
tmp <- data.frame("data"=numeric(),"model"=numeric(),"Z_residual"=character(),"Site"=numeric())
for (j in 1:dim(sims_vc)[2]) {
  nv <- dim(sims_vc)[1]*(1-q)
x <- sort(as.data.frame(Yboot$X) %>% pull(all_of(j)),decreasing=TRUE)[1:nv]
y1 <- sort(sims_vc %>% pull(all_of(j)),decreasing=TRUE)[1:nv]
y2 <- sort(sims_gauscop %>% pull(all_of(j)),decreasing=TRUE)[1:nv]
tmp <- rbind(tmp,rbind(data.frame("data"=x,"model" = y1,"Z_dependence"="vine_copula","Site"=j),data.frame("data"=x,"model" = y2,"Z_dependence"="gaussian_copula","Site"=j)))
}
p <- ggplot(tmp) + geom_point(aes(x=data,y=model,col=Z_dependence,alpha=Site)) + geom_abline(slope=1,linetype="dashed")
ggsave(p,filename="QQ_plot_simulated_original.pdf",width=7,height=7)

# try also comparing bias
library(Metrics)
bias_vc <- bias_gauscop <- c()
for (j in 1:dim(sims_vc)[2]) {
  nv <- dim(sims_vc)[1]*(1-q)
  x <- sort(as.data.frame(Yboot$X) %>% pull(all_of(j)),decreasing=TRUE)[1:nv]
  y1 <- sort(sims_vc %>% pull(all_of(j)),decreasing=TRUE)[1:nv]
  y2 <- sort(sims_gauscop %>% pull(all_of(j)),decreasing=TRUE)[1:nv]
  bias_vc <- append(bias_vc,mean(y1-x) )
  bias_gauscop <- append(bias_gauscop,mean(y2-x) )
}

A <- matrix(bias_vc, nrow = 5, ncol = 5)  
print(A)

# Plot a heatmap 
heatmap(A,Rowv=NA,Colv=NA,col=heat.colors(3))
tmp <- data.frame("Vine.copula"=bias_vc,"Gaussian.copula"=bias_gauscop,"Bias.difference" = bias_vc-bias_gauscop,"latitude"=rep(1:5,n=5),"longitude" = rep(1:5,each=5)) %>% pivot_longer(cols = c("Vine.copula","Gaussian.copula","Bias.difference"),names_to ="Z_dependence",values_to = "bias") %>% mutate("Z_dependence"=factor(Z_dependence,levels=c("Vine.copula","Gaussian.copula","Bias.difference")))
levels(tmp$Z_dependence) <- c("vine copula","Gaussian copula","bias difference")

p <- ggplot(tmp, aes(x = longitude, y = latitude, fill = bias)) +
  geom_tile() +
  scale_fill_gradientn(colours=c( "#009ADA", "#fffffe",
                     "#C11432")
                      ) +
  geom_text(aes(label = round(bias,2)), color = "black", size = 5) +
  theme(text = element_text(size = 15)) + facet_wrap("Z_dependence") + coord_fixed()

ggsave(p,filename="bias_gauscop_vc.pdf",width=15,height=5)

for (i in 1:25) {
  print(summary(cond_modelvc1[[i]][[4]]))
}

# compare for a range of thresholds
u <- c(0.9,0.95,0.99)
thres_est <- list()
for (ui in 1:length(u)) {
  cond_model_fit_wrapper = function(i){
    return(fn(site_index=i,v=u[ui],Yrun1Lap =Yboot$Y,res_dist = "AGG_Gauscopula"))
  }
  thres_est[[ui]] <- sapply(1:25,cond_model_fit_wrapper,simplify=FALSE)

}

# compare alpha and beta values
thres_est[[1]][[1]][[1
                     ]]

# compare for alpha and beta
tmpall <- data.frame("threshold"=numeric(),"a"=numeric(),"b"=numeric(),"cond_site"=numeric(),"res_site" = numeric())
for (site in 1:25) {
  a90 <- thres_est[[1]][[site]][[1]]$a
  a95 <- thres_est[[2]][[site]][[1]]$a
  a99 <- thres_est[[3]][[site]][[1]]$a
  b90 <- thres_est[[1]][[site]][[1]]$b
  b95 <- thres_est[[2]][[site]][[1]]$b
  b99 <- thres_est[[3]][[site]][[1]]$b
  tmpa <- data.frame(a90,a95,a99)
  tmpb <- data.frame(b90,b95,b99)
  print(summary(tmpa))
  print(summary(tmpb))
  tmpa <- tmpa %>% pivot_longer(cols=everything(),names_to = "threshold",values_to = "a") %>% mutate(threshold=factor(recode(threshold,a90=0.90,a95=0.95,a99=0.99),levels=c(0.90,0.95,0.99)))
  tmpb <- tmpb %>% pivot_longer(cols=everything(),names_to = "threshold",values_to = "b") %>% mutate(threshold=factor(recode(threshold,b90=0.90,b95=0.95,b99=0.99),levels=c(0.90,0.95,0.99)))
  tmp <- cbind(tmpa,tmpb %>% select(b)) %>% mutate("cond_site"=site) %>% mutate("res_site" = rep(c(1:25)[-site],each=3))
 p <-  ggplot(tmp) + geom_point(aes(x=a,y=b,colour=threshold))
 print(p)
 tmpall <- rbind(tmpall,tmp)
}
library(latex2exp)
library(gridExtra)
tmpall <- tmpall %>% mutate("s_s0"=rep(1:(nrow(tmpall)/3),each=3)) 
p <-  ggplot(tmpall) + geom_line(aes(x=a,y=b,group=s_s0),alpha=0.5,linewidth=0.3)+ geom_point(aes(x=a,y=b,colour=threshold),size=0.6) + scale_colour_manual(values = c("0.9"="#009ADA","0.95"="#FDD10A","0.99"= "#C11432")) + xlab(TeX("$\\alpha_{s|s_0}$")) + ylab(TeX("$\\beta_{s|s_0}$")) + coord_equal()
ggsave(p,filename="../ab_threshold_compare.pdf",width=6,height=5)

# select high alpha
p1 <-  ggplot(tmpall %>% filter(a>0.8)) + geom_line(aes(x=a,y=b,group=s_s0),alpha=0.5)+ geom_point(aes(x=a,y=b,colour=threshold),size=1.2) + scale_colour_manual(values = c("0.9"="#009ADA","0.95"="#FDD10A","0.99"= "#C11432")) + xlab(TeX("$\\alpha_{s|s_0}$")) + ylab(TeX("$\\beta_{s|s_0}$")) + coord_fixed()
p1
ggsave(grid.arrange(p,p1,ncol=2),filename="../ab_threshold_compare2.pdf",width=10,height=5)

