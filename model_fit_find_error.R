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


# obtain bootstrap using model refitting -----
model_refit <- function(k) {
  source("Block_Bootstrapping.R")
  tq_table_run <- list()
  for (run_number in 1:4) {
    set.seed(k*100+run_number)
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
    
    cond_model_fit_wrapper = function(i){
      return(fn(site_index=i,v=q,Yrun1Lap =Yboot$Y,res_dist = "empirical"))
    }
    start_time <- Sys.time()
    cond_modelvc1 <- sapply(1:25,cond_model_fit_wrapper,simplify=FALSE)
    end_time <- Sys.time()
    end_time - start_time 
    Nrun <- 1
    # # simulate from cond. model with AGG_vinecopula residuals
    # boot <- replicate(n = Nrun, expr = sim_cond_model(Yrun=Yboot$Y,cond_model=cond_modelvc1,res_dist = "AGG_vinecopula",Y_boot=Yboot), simplify = FALSE)
    # # evaluate target quantities for all 6 questions (preliminary+target)
    # tq <- do.call(rbind,lapply(bootemp,FUN = Qeval))
    # # sum and divide by Nrun
    # temp <- apply(tq,MARGIN=2,FUN=sum)/Nrun
    
    # simulate and evaluate in one function
    simulate_evaluate <- function(x="AGG_vinecopula") {
      sims <- sim_cond_model(Yrun=Yboot$Y,cond_model=cond_modelvc1,res_dist = x,Y_boot=Yboot)
      y <- Qeval(sims)
      # v <- seq(1, 1.7, by = 0.1)
      # q1 <- rbind(sapply(X = v, FUN = Q1, Yrun = sims))
      # 
      # v <- seq(2.1,5.7,by=0.6)
      # q2 <- rbind(sapply(X = v, FUN = Q2, Yrun = sims))
      # 
      # v <- seq(2,5,by=0.5)
      # q3 <- rbind(sapply(X = v, FUN = Q3, Yrun = sims))
      
      return(list(y))
    }
    
    # evaluate target quantities for all 6 questions (preliminary+target)
    start_time <- Sys.time()
    tq <- replicate(Nrun,expr = simulate_evaluate(x="empirical"), simplify=TRUE)
    end_time <- Sys.time()
    end_time - start_time
    
    tq_matrix <- unlist(tq) %>% matrix(nrow=Nrun, byrow=TRUE)
    
    # sum and divide by Nrun
    tq_means <- apply(tq_matrix,MARGIN=2,FUN=sum)/Nrun
    tvc <- data.frame(matrix(tq_means, 1))
    names(tvc) <- names(tq[[1]])
    
    # # simulate from cond.model with empirical residuals
    # bootemp <- replicate(n = Nrun, expr = sim_cond_model(Yrun=Yboot$Y,cond_model=cond_modelvc1e,res_dist = "empirical",Y_boot=Yboot), simplify = FALSE)
    # # evaluate target quantities for all 6 questions (preliminary+target)
    # tq <- do.call(rbind,lapply(boot,FUN = Qeval))
    # # sum and divide by Nrun
    # tvc <- apply(tq,MARGIN=2,FUN=sum)/Nrun
    
    # print in latex table
    tq_table <- cbind(data.frame("Residual method" = c("empirical","vine copula")),
                      rbind(temp,tvc))
    tq_table_run[[run_number]] <- tq_table 
  }
  # return a list of lists
  return(tq_table_run)
  
}

d = 25
q = 0.95 # quantile threshold for fitting the conditional model
Nboot = 10
model_refit(1)
#param_bootstrap = mclapply(1:Nboot, model_refit, mc.cores = 2)


