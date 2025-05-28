## Function to fit threshold wrt time
Model_Year <- function(t, beta_1, beta_2, beta_3){
  if(!is.numeric(t) | min(t) < 0){stop("t must be a positive real number")}
  if(!is.numeric(beta_1)){stop("beta_1 must be a real number")}
  if(!is.numeric(beta_2)){stop("beta_2 must be a real number")}
  if(!is.numeric(beta_3)){stop("beta_3 must be a real number")}
  return(beta_1 + beta_2 * cos((2 * pi * t / 365) - beta_3))
}

Fit_GPD_Seasonal_Threshold <- function(df, tau){
  ## Check df is a data frame
  if(!is.data.frame(df)){
    stop("df must be a data frame")
  }
  ## Check tau is a single number in the region (0, 1)
  if(!is.numeric(tau) | length(tau) != 1 | min(tau) <= 0 | max(tau) >= 1){
    stop("tau must be a single positive number in the region (0,1)")
  }
  
  ## Fit the model
  fit <- try(nlrq(formula = Value ~ Model_Year(t = Time, beta_1 = beta_1, beta_2 = beta_2, beta_3 = beta_3),
                  data = df, 
                  start = list(beta_1 = 0, beta_2 = 1, beta_3 = 0), 
                  tau = tau),
             silent = TRUE)
  
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Fit_GPD_Seasonal_Threshold")
    out <- list()
    out$beta <- rep(NA, 3)
    out$u <- rep(NA, nrow(df))
  }
  else{
    ## Extract the output
    pars <- coef(fit)
    
    out <- list()
    out$beta <- pars
    out$u <- Model_Year(t = df$Time, beta_1 = pars[1], beta_2 = pars[2], beta_3 = pars[3])
  }
  return(out)
}

## Function to fit the non-stationary GPD wrt to time
Loglike_GPD_Seaonsal_Scale <- function(par, x, t, negative = FALSE){
  
  ## Extract parameter estimates
  scale <- exp(Model_Year(t = t, beta_1 = par[1], beta_2 = par[2], beta_3 = par[3]))
  shape <- par[4]
  
  if(shape < 0){
    if(any(1 + shape*x/scale <= 0)){
      return((-10^10)*(-1)^negative)
    }
  }
  
  if(abs(shape) <= 1e-10){
    z <- sum(-log(scale) - x/scale)
  }
  else{
    ## Calculate the log-likelihood
    z <- sum(-log(scale) - (1 + 1/shape)*log(pmax(0, 1 + shape*(x/scale)))) 
  }
  
  ## Check the output
  if(is.na(z) | is.nan(z) | is.infinite(z)){
    return((-10^10)*(-1)^negative)
  }
  else{
    return(z*(-1)^negative)
  }
}

## Function to fit the non-stationary GPD wrt respect to time
Fit_GPD_Seasonal_Scale <- function(df){
  ## Check df is a data frame
  if(!is.data.frame(df)){
    stop("df must be a data frame")
  }
  
  ## Fit the model
  fit <- try(optim(par = c(0, 1, 0, 0.1),
                   fn = Loglike_GPD_Seaonsal_Scale,
                   x = df$Excess,
                   t = df$Time,
                   negative = TRUE,
                   method = "BFGS"),
             silent = TRUE)
  
  ## Return the output
  if(inherits(fit, "try-error")){
    warning("Error in optim call from Fit_GPD_Seasonal_Scale")
    out <- list()
    out$par <- list(beta = rep(NA, 3), scale = rep(NA, nrow(df)), shape = NA)
    out$value <- NA
    out$convergence <- NA
    out$counts <- NA
    out$message <- NA
  }
  else if(fit$convergence != 0 | fit$value == 1e+10){
    warning("Non-convergence from Fit_GPD_Seasonal_Scale")
    out <- list()
    out$par <- list(beta = rep(NA, 3), scale = rep(NA, nrow(df)), shape = NA)
    out$value <- NA
    out$convergence <- NA
    out$counts <- NA
    out$message <- NA
  }
  else if(!is.na(fit$par[1])){
    out <- list()
    out$par <- list(beta = fit$par[-4],
                    scale = exp(Model_Year(t = df$Time, beta_1 = fit$par[1], beta_2 = fit$par[2], beta_3 = fit$par[3])), 
                    shape = fit$par[4])
    out$llh <- -fit$value
    out$convergence <- fit$convergence
    out$counts <- fit$counts
    out$message <- fit$message
  }
  else{
    warning("Unknown error from Fit_GPD_Seasonal_Scale")
    out <- list()
    out$par <- list(beta = rep(NA, 3), scale = rep(NA, nrow(df)), shape = NA)
    out$value <- NA
    out$convergence <- NA
    out$counts <- NA
    out$message <- NA
  }
  return(out)
}

## Fit the threshold and scale parameter 
Fit_GPD_Seasonal_Threshold_Scale <- function(df, tau){
  ## Check df is a data frame
  if(!is.data.frame(df)){
    stop("df must be a data frame")
  }
  ## Check tau is a single number in the region (0, 1)
  if(!is.numeric(tau) | length(tau) != 1 | min(tau) <= 0 | max(tau) >= 1){
    stop("tau must be a single positive number in the region (0,1)")
  }
  
  ## First fit the Threshold
  fit_u <- Fit_GPD_Seasonal_Threshold(df = df, tau = tau)
  if(any(is.na(fit_u$beta))){
    warning("Error in optim call from Fit_GPD_Seasonal_Threshold")
    
    ## Obtain the output
    out <- list()
    out$par <- list(beta_u = rep(NA, 3),
                    beta_scale = rep(NA, 3),
                    u = rep(NA, nrow(df)),
                    scale = rep(NA, nrow(df)),
                    shape = NA)
    out$llh <- NA
    out$convergence <- NA
    out$counts <- NA
    out$message <- NA
    out$excess_data <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(df)))
    colnames(out$excess_data) <- colnames(df)
  }
  else{
    ## Add the threshold to the data frame and filter for positive excesses
    df <- filter(mutate(df, 
                        u = fit_u$u, 
                        Excess = Value - u), 
                 Excess > 0)
    
    ## Fit the seasonal scale
    fit_scale <- Fit_GPD_Seasonal_Scale(df)
    if(is.na(fit_scale$convergence)){
      warning("Error in optim call from Fit_GPD_Seasonal_Scale")
      
      ## Obtain the output
      out <- list()
      out$par <- list(beta_u = rep(NA, 3),
                      beta_scale = rep(NA, 3),
                      u = rep(NA, nrow(df)),
                      scale = rep(NA, nrow(df)),
                      shape = NA)
      out$llh <- NA
      out$convergence <- NA
      out$counts <- NA
      out$message <- NA
      out$excess_data <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(df)))
      colnames(out$excess_data) <- colnames(df)
      
    }
    else{
      
      ## Obtain the output
      out <- list()
      out$par <- list(beta_u = fit_u$beta,
                      beta_scale = fit_scale$par$beta,
                      u = fit_u$u,
                      scale = fit_scale$par$scale,
                      shape = fit_scale$par$shape)
      out$llh <- fit_scale$llh
      out$convergence <- fit_scale$convergence
      out$counts <- fit_scale$counts
      out$message <- fit_scale$message
      out$excess_data <- df
    }
  }
  
  ## Return the output
  return(out)
}

QQ_plot_Seasonal_Threshold_Scale <- function(fit, B = 100){
  ## check inputs
  if(!is.list(fit)){
    stop("fit must be a list of output from Threshold_Selection_Seasonal_Threshold_Scale")
  }
  if(B <= 0 | B%%1 != 0){
    stop("Number of bootstrapped samples must be a positive integer")
  }
  
  ## Obtain the empirical qunatiles
  m <- nrow(fit$excess_data)
  p <- (1:m)/(m+1)
  q_exp <- transform_to_exp(y = fit$excess_data$Excess,
                            sig = fit$par$scale,
                            xi = fit$par$shape)
  q_emprical <- quantile(q_exp, probs = p)
  
  ## Obtain the model quantiels
  q_model <- qexp(p, rate = 1)
  
  ## Obtain the confidence intervale for the model quantiles
  
  ## Bootstrap the data
  extreme_data_boot <- replicate(n = B,
                                 expr = slice_sample(fit$excess_data, 
                                                     n = m, 
                                                     replace = TRUE),
                                 simplify = FALSE)
  
  ## Fit the model to the bootstrapped data
  boot_fits <- lapply(FUN = Fit_GPD_Seasonal_Scale,
                      X = extreme_data_boot)
  
  ## Extract the MLE
  scale_mle <- lapply(boot_fits, function(x){x$par$scale})
  shape_mle <- lapply(boot_fits, function(x){x$par$shape})
  
  ## Transform the data onto standard exponential margins
  y_exp <- do.call(cbind, pmap(.l = list(y = extreme_data_boot,
                                         sig = scale_mle,
                                         xi = shape_mle),
                               .f = function(y, sig, xi){
                                 transform_to_exp(y$Excess, sig, xi)}))
  
  ## Compute bootstrap confidence intervals
  y_exp <- apply(y_exp, 2, sort)
  boot_ci <- t(apply(y_exp, 1, quantile, probs = c(0.025, 0.975)))
  
  ## Create the plot
  plot(x = 1, type = "n",
       xlab = "Model Quantiles",
       ylab = "Empirical Quantiles",
       main = sprintf("Location (%d, %d) - Run %d", 
                      unique(fit$excess_data$Grid_1), unique(fit$excess_data$Grid_2), unique(fit$excess_data$Run)),
       xlim = c(0, max(q_emprical)),
       ylim = c(0, max(boot_ci) + 0.1))
  
  polygon(x = c(q_model, rev(q_model)),
          y = c(boot_ci[, 1], rev(boot_ci[, 2])),
          col = 'grey80', border = NA)
  
  points(q_model, q_emprical, pch = 19)
  
  abline(0, 1, col = 2, lwd = 2)
}