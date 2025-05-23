# libraries ---------
library(parallel)
library(dplyr)
library(texmex)
library(evd)
library(rvinecopulib)
library(tidyverse)
library(VGAM)

# function scripts ------
source('conditional_model_helpers.R')
source('conditional_model_fit.R')
source('target_quantity_helpers.R')
source('model_diagnostics.R')

### loading data file
load('Yrun.RData')
Yrun1Lap <- as.data.frame((Yrun1 %>% apply(c(2),FUN=row_number))/(nrow(Yrun1)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()

# marginal and conditional fits and data -----
load('cond_modelvc.RData')

# obtain bootstrap using model refitting -----
model_refit <- function(k) {
  Yboot <- sim_cond_model(Yrun=Yrun1,cond_model=cond_modelvc1,res_dist = "AGG_vinecopula")
  Ybootemp <- sim_cond_model(Yrun=Yrun1,cond_model=cond_modelvc1,res_dist = "empirical")

  cond_model_fit_wrapper = function(i){
    return(fn(site_index=i,v=q,Yrun=Yboot,res_dist = "AGG_vinecopula"))
  }
  cond_modelvc1 <- sapply(1:25,cond_model_fit_wrapper,simplify=FALSE)

cond_model_fit_wrapper = function(i){
    return(fn(site_index=i,v=q,Yrun=Ybootemp,res_dist = "empirical"))
  }
  cond_modelvc1e <- sapply(1:25,cond_model_fit_wrapper,simplify=FALSE)
  
  Nrun <- 50
   boot <- replicate(n = Nrun, expr = sim_cond_model(Yrun=Yboot,cond_model=cond_modelvc1,res_dist = "AGG_vinecopula"), simplify = FALSE)

  bootemp <- replicate(n = Nrun, expr = sim_cond_model(Yrun=Ybootemp,cond_model=cond_modelvc1e,res_dist = "empirical"), simplify = FALSE)

  # evaluate target quantities

  tq <- do.call(rbind,lapply(bootemp,FUN = Qeval))
  # sum and divide by Nrun
  temp <- apply(tq,MARGIN=2,FUN=sum)/Nrun
  
  # evaluate target quantities
  tq <- do.call(rbind,lapply(boot,FUN = Qeval))
  # sum and divide by Nrun
  tvc <- apply(tq,MARGIN=2,FUN=sum)/Nrun

  
  # print in latex table
  tq_table <- cbind(data.frame("Residual method" = c("empirical","vine copula")),
                    rbind(temp,tvc))

  # try model diagnostics
  diag1 <- Q1diag(sims1=boot,sims2=bootemp)
  diag2 <- Q2diag(sims1=boot,sims2=bootemp)
  diag3 <- Q3diag(sims1=boot,sims2=bootemp)
  return(list(tq_table,diag1,diag2,diag3))
}

d = 25
q = 0.95
Nrun = 100
# model_refit(1)
param_bootstrap = mclapply(1:Nrun, model_refit, mc.cores = 100)

# save results -----
save(param_bootstrap, file = "param_bootstrap3.RData")
