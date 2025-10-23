# libraries ---------
library(parallel)
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
    sims_vc <- sim_cond_model(Yrun=Yboot$Y,cond_model=cond_modelvc1,res_dist = "AGG_vinecopula",Y_boot=Yboot)
    sims_emp <- sim_cond_model(Yrun=Yboot$Y,cond_model=cond_modelvc1,res_dist = "empirical",Y_boot=Yboot)
    
    save(Yboot,cond_modelvc1,sims_vc,sims_emp,file="../run1_bootstrap_resfit.RData")
    