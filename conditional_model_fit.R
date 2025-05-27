fn <- function(site_index=1,v=0.95,Yrun=Yrun1,res_dist="empirical") {
  # transform to Laplace margins
  Yrun1Lap <- as.data.frame((Yrun %>% apply(c(2),FUN=row_number))/(nrow(Yrun)+1)) %>% apply(c(1,2),FUN=unif_laplace_pit) %>% as.data.frame()
  j <- site_index
  pe_cond1 <- par_est(df=Yrun1Lap,v=v,given=j,keef_constraints = c(1,2))
  # calculate a vector of observed residuals
  Y_given1extreme <- Yrun1Lap %>% filter(Yrun1Lap[,j]>quantile(Yrun1Lap[,j],v))
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
    # 3b. fit a vine copula
    res_vc <- rvinecopulib::vinecop( (Z %>% apply(c(2),FUN=row_number))/(nrow(Yrun1Lap)*(1-v)+1),selcrit = "mbicv")  
    return(list(pe_cond1,Z,res_margin,res_vc))
  }
}

sim_cond_model <- function(Nrun=1,Yrun=Yrun1,q=0.95,res_dist="empirical",cond_model=cond_model,Y_boot=Yboot) {
  ## Reading in required packages
  v <- qlaplace(q)
  Nsim <- sum(apply(Yrun1Lap, 1, max) > v)*Nrun
  n_sim <- Nsim*10
  required_pckgs <- c("Matrix", "LaplacesDemon")
  t(t(sapply(required_pckgs, require, character.only = TRUE)))
  
  ## Get the conditioning sites
  cond_sites <- sample(x = 1:d, size = n_sim, replace = TRUE)
  cond_sites_unique <- sort(unique(cond_sites))
  n_cond_sites <- as.numeric(table(cond_sites))
  ## Convert the threshold onto Laplace margins and simulate from the conditioning site
  
  large_Laplace <- v +  rexp(n = n_sim, rate = 1)
  ## Sample residuals from the fitted model
  if (res_dist=="empirical") {
    Resid <- lapply(cond_model, function(x){x[[2]]})
    # 2. simulate residuals from an empirical distribution
    Z_sample <- lapply(cond_sites_unique, function(i){
      unlist(unname(as.matrix(Resid[[i]])[sample(1:nrow(Resid[[i]]), n_cond_sites[i], replace = TRUE),]))})
  } else if (res_dist=="AGG_vinecopula") {
    # AGG parameters
    AGGpar <- lapply(cond_model, function(x){x[[3]]})
    # vine copula objects
    Resid_vc <- lapply(cond_model, function(x){x[[4]]})
    # 2. simulate residuals from a vine copula
    Zsim_unif <- lapply(cond_sites_unique, function(i){
      rvinecopulib::rvinecop(n=n_cond_sites[i],vine=Resid_vc[[i]])})
    # 3. convert to original margins
    # rewrite using apply
    AGG_inverse_wrapper <- function(i,j){
      return(qAGG(p=as.numeric(unlist(Zsim_unif[[i]][,j])),theta=as.numeric(unlist(AGGpar[[i]][j,2:6]))))
    }
    # transform a vector of residuals for each of the 25 models  
    Z_sample <- lapply(cond_sites_unique, function(i){
      sapply(1:ncol(Zsim_unif[[i]]), AGG_inverse_wrapper,i=i)})
  } # end of simulating from AGG_vinecopula residual distribution
  
  ## Convert the simulated data onto Laplace margins
  a_cond_sites <- lapply(cond_sites_unique, function(i){cond_model[[i]][[1]]$a})
  b_cond_sites <- lapply(cond_sites_unique, function(i){cond_model[[i]][[1]]$b})
  
  Laplace_Samples <- vector("list", d)
  for(i in 1:d){
    large_Laplace_cond <- large_Laplace[cond_sites == i]
    Y_Yi_large <- sapply(a_cond_sites[[i]], function(a){a*large_Laplace_cond}) + 
      sapply(b_cond_sites[[i]], function(b){large_Laplace_cond^b})*Z_sample[[i]]
    
    ## Now combine Y{-i}|Y_{i} > u_{Y_{i}} and Y{-i}|Y_{i} > u_{Y_{i}}
    if(i == 1){
      Laplace_Samples[[i]] <- cbind(large_Laplace_cond, Y_Yi_large)
    }
    else if(i == d){
      Laplace_Samples[[i]] <- cbind(Y_Yi_large, large_Laplace_cond)
    }
    else{
      Laplace_Samples[[i]] <- cbind(Y_Yi_large[,1:(i-1)], large_Laplace_cond, Y_Yi_large[,i:(d-1)])
    }
  }
  Laplace_Samples_Full <- do.call(rbind, Laplace_Samples)
  
  ## Now use importance sampling to resample the simulated data
  ## Idea is to up-weight those close to the boundary and down-weight those in the
  ## centre of the spatial process
  weights <- d/apply(Laplace_Samples_Full, 1, function(x){length(which(x > v))})
  Index <- sample(x = 1:nrow(Laplace_Samples_Full), size = Nsim, prob = weights/sum(weights),replace=TRUE)
  ## Now we convert the data onto the original scale
  Final_Laplace_Samples <- Laplace_Samples_Full[Index,]
  ## Above simulated Y(s) | max(Y(s)) > v
  ## We want to get the unconditioned process
  ## First get the data points where we have no extremes
  data_Laplace <- Yrun1Lap
  Index_No_Extremes <- which(apply(data_Laplace, 1, max) < v)
  Data_Body_Laplace <- data_Laplace[Index_No_Extremes,]
  
  ## then get the probability we will draw from the joint tail
  #p_tail <- mean(apply(data_Laplace, 1, max) > v)
  #p_accecpt <- runif(n = nrow(Yrun))
  
  ## Get the index of body/tail
  Index_Body <- sample(x = 1:length(Index_No_Extremes), size = nrow(Yrun)*Nrun-Nsim, replace = TRUE)
  #Index_Tail <- which(p_accecpt < p_tail)
  
  Final_Laplace_Samples <- as.data.frame(Final_Laplace_Samples)
  names(Final_Laplace_Samples) <- names(Data_Body_Laplace)
  ## Get the final data sets
  Data_Final_Laplace_Margins <- rbind(Data_Body_Laplace[Index_Body,],
                                      Final_Laplace_Samples)
  
  final_uniform_data <- apply(X=Data_Final_Laplace_Margins,MARGIN=2,plaplace)
  final_uniform_data <- as.data.frame(final_uniform_data)
  
# reorder the data
ordering <- rank(rowSums(apply(Y_boot$X, 2, rank)) ## Y_boot$X blocked bootstrapped data

final_uniform_data <- final_uniform_data[ordering,] ## reorder for temporal ordering

## Obtain the data on the original margins
Data_orig_margins <- sapply(1:25, function(i){
    qspliced_nonsta(p = final_uniform_data[,i],
                      t = 1:nrow(final_uniform_data,
                      x = Y_boot$X[,i],
                      gpd_par = list(u = Y_boot$par[i,1:3],
                                     scale = Y_boot$par[i,4:6],
                                     shape = Y_boot$par[i,7])
  
  names(Data_orig_margin) <- names(Data_Final_Laplace_Margins)
  
  return(Data_orig_margin)
  #  return(final_uniform_data)
}

cond_model_fit_wrapper = function(i){
  return(fn(site_index=i,v=q,Yrun=Yrun1))
}
