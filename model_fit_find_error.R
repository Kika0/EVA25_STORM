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

d = 25
q = 0.95 # quantile threshold for fitting the conditional model
Nboot = 50
k <- 58

fn <- function(site_index=1,v=0.95,Yrun1Lap=Yboot$Y,res_dist="empirical") {
  j <- site_index
  pe_cond1 <- par_est(df=Yrun1Lap,v=v,given=j,keef_constraints = c(1,2))
  # calculate a vector of observed residuals
  Y_given1extreme <- Yrun1Lap %>% dplyr::filter(Yrun1Lap[,j]>quantile(Yrun1Lap[,j],v))
  Z <- as.data.frame(matrix(nrow=nrow(Y_given1extreme),ncol=0))
  Y1 <- Y_given1extreme[,j]
  res <- c(1:ncol(Y_given1extreme))[-j]
  for (i in 2:ncol(Y_given1extreme)) {
    Y2 <- Y_given1extreme[,res[i-1]]
    Z <- cbind(Z,data.frame("Z"= (Y2-pe_cond1$a[i-1]*Y1)/(Y1^pe_cond1$b[i-1]) ) )
  }
  names(Z) <- paste0("Z",res)
  
  if (res_dist=="empirical") {
    return(list(pe_cond1,Z))
  }
  if (res_dist=="AGG_vinecopula") {
    # 3a. estimate parameters for AGG residual margins
    res_margin <- res_margin_par_est(obs_res=Z,method = "AGG")
    # transform to AGG margins
    pAGG_wrapper <- function (i) {
      pAGG(x=Z[i],theta=as.numeric(unlist(res_margin[i,2:6])))
    }
    Z_AGG <- sapply(1:ncol(Z),FUN=pAGG_wrapper)
    # 3b. fit a vine copula
    # res_vc <- rvinecopulib::vinecop( Z_AGG,selcrit = "mbicv")  
    return(list(pe_cond1,Z,res_margin))
  }
}


# obtain bootstrap using model refitting -----
  source("Block_Bootstrapping.R")
  tq_table_run <- list()
  cond_modelvc1 <- list()
  for (run_number in 1:4) {
    set.seed(k*100+run_number)
    start_time <- Sys.time()
    Y_Bootstrapped <- NonSta_GPD_to_Lapalce(df = data,
                                            run = run_number,
                                            tau = 0.98,
                                            B = 1,
                                            l = 2)
    end_time <- Sys.time()
  print(end_time - start_time )
    Yboot <- Y_Bootstrapped[[1]]
    names(Yboot$Y) <- paste0("Y",1:25)
    Yboot$Y <- as.data.frame(Yboot$Y)
    
    cond_model_fit_wrapper = function(i){
      return(fn(site_index=i,v=q,Yrun1Lap =Yboot$Y,res_dist = "AGG_vinecopula"))
    }
    start_time <- Sys.time()
    cond_modelvc1[[run_number]] <- sapply(1:25,cond_model_fit_wrapper,simplify=FALSE)
    end_time <- Sys.time()
  print(  end_time - start_time )
}

