#' Transform from Uniform to Laplace margins
#'
#' @param x A number sampled from U(0,1) distribution.
#'
#' @return x transformed to Laplace margin
#' @export
#'
#' @examples unif_laplace_pit(x=runif(1))
unif_laplace_pit <-  function(x) {
  if (x<0.5) { y <-log(2*x) }
  else { y <- -log(2*(1-x)) }
  return(y)
}

#' Transform from Laplace to uniform margins
#'
#' @param x A number sampled from U(0,1) distribution.
#'
#' @return x transformed to Laplace margin
#' @export
#'
#' @examples laplace_unif_pit(x=runif(1))
laplace_unif_pit <-  function(x) {
  if (x<0.5) { y <- exp(x/2) }
  else { y <- 1-exp(-x)/2 }
  return(y)
}

#' CDF of AGG margins for residuals
#'
#' @param x A vector on AGG margins.
#'
#' @return y A vector of probabilities from cdf of AGG distribution.
#' @export
#'
#' @examples F_AGG(x=runif(1))
pAGG <- function(x,theta) {
  y <- numeric(length(x))
  mu <- theta[1]
  sigl <- theta[2]
  sigu <- theta[3]
  deltal <- theta[4]
  deltau <- theta[5]
  C_AGG <-  (sigl/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
  y[x<mu] <- C_AGG*sigl/deltal*gamma(1/deltal) *( 1-pgamma(q=  ((mu-x[x<mu])/sigl)^deltal ,shape=1/deltal,rate=1)  ) 
  y[x>=mu] <- C_AGG*sigl/deltal*gamma(1/deltal) + C_AGG*sigu/deltau*gamma(1/deltau)*(pgamma(((x[x>=mu]-mu)/sigu)^deltau,shape=1/deltau,rate=1))
  return(y)
}

#' Inverse CDF of AGG margins for residuals
#'
#' @param x A vector of probabilities.
#'
#' @return y A vector on AGG margins.
#' @export
#'
#' @examples qAGG(p=seq(0,1,length.out=100),theta=c(0,1,1,2,2))
qAGG <- function(p,theta) {
  y <- numeric(length(p))
  mu <- theta[1]
  sigl <- theta[2]
  sigu <- theta[3]
  deltal <- theta[4]
  deltau <- theta[5]
  C_AGG <-  (sigl/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
  pmu <- C_AGG*sigl/deltal*gamma(1/deltal)
  y[p<pmu] <- -qgamma(1-(p[p<pmu]*deltal/(C_AGG*sigl*gamma(1/deltal))),shape = 1/deltal,rate = 1)^(1/deltal)*sigl+mu
  y[p>=pmu] <- qgamma( (p[p>=pmu]/C_AGG - sigl/deltal*gamma(1/deltal))*deltau/(sigu*gamma(1/deltau)),shape=1/deltau,rate=1)^(1/deltau)*sigu+mu
  return(y)
}

#' Calculate negative log-likelihood of Normal regression
#' 
#' Conditioning on Y_1 being extreme to model Y_2 (conditional model in bivariate case)
#'
#' @param theta A set of 4 parameters: a,b,mu,sig.
#' @param df A dataset with column names of paste0("Y",number).
#' @param given A numeric specifying column name of cond. variable Y1.
#' @param sim A numeric specifying column name of other variable Y2.
#' @return A numeric negative log-likelihood.
#' @export
#'
#' @examples
Y_NLL <- function(theta,df=Y_given1extreme,given=1,sim=2,a_hat=NULL,b_hat=NULL,b_max=1) {
  if (is.null(a_hat)==FALSE) {
    a <- a_hat
  } else {a <- theta[1]}
  if (is.null(b_hat)==FALSE) {
    b <- b_hat
  } else {b <- theta[length(theta)-2]}
  mu <- theta[length(theta)-1]
  sig <- theta[length(theta)]
  Y1 <- df %>% dplyr::select(paste0("Y",given)) %>% pull()
  Y2 <- df %>% dplyr::select(paste0("Y",sim)) %>% pull()
  if (a<(-1) | a>1 | b<0 | b>=b_max | sig<0) {
    log_lik <- (10^6) # low log-likelihood outside Keef bounds
  }
  else {
    log_lik <- sum(log(Y1^b *sig*sqrt(2*pi)) + ((Y2-a*Y1-mu*Y1^b)^2/(2*(Y1^b*sig)^2))  )
  }
  return(log_lik)
}

#' NLL for AGG with different both scale and shape for lower and upper tail
#'
#' @param x A numerical vector of data.
#' @param theta A vector of parameters c(mu,sigl,sigu,deltal,deltau).
#'
#' @return A number of negative log-likelihood.
#' @export
#'
#' @examples NLL_AGG(x=rnorm(50),theta=c(0,1,1,2,2))
NLL_AGG <- function(x,theta) {
  mu <- theta[1]
  sigl <- theta[2]
  sigu <- theta[3]
  deltal <- theta[4]
  deltau <- theta[5]
  z <- c()
  if(sigl<=0 | sigu<=0 | deltal<=0 |deltau<=0 ){return(10e10)}
  C_AGG <-  (sigl/deltal*gamma(1/deltal) + sigu/deltau*gamma(1/deltau)  )^(-1)
  for (i in 1:length(x)) {
    if (x[i]<mu) {
      #   z[i] <- C_AGG*exp(-abs((x[i]-mu)/sigl)^deltal)
      z[i] <- log(C_AGG)-((mu-x[i])/sigl)^deltal 
    }
    else 
      #   z[i] <- C_AGG*exp(-abs((x[i]-mu)/sigu)^deltau)
    {z[i] <- log(C_AGG)-((x[i]-mu)/sigu)^deltau }
  }
  return(-sum(z))
}


# keef constraints for beta
keef_constraint1 <- function(b,a,Y1,Y2) {
  if (b<0 | b>1) {return(10^6)}
  v <- max(Y1)
  ZmAI <- max((Y2-a*Y1)/(Y1^b))
  ZmAD <- max(Y2-Y1)
  return((1-b*ZmAI*v^(b-1)  -a)^2)
}

keef_constraint2 <- function(b,a,Y1,Y2) {
  if (b<0 | b>1) {return(10^6)}
  v <- max(Y1)
  ZmAI <- max((Y2-a*Y1)/(Y1^b))
  ZmAD <- max(Y2-Y1)
  return((1-v^(b-1)*ZmAI+v^(-1)*ZmAD  -a)^2)
}


# generate a table of parameter estimates conditional on (given) each of the specified vector of variables
par_est <- function(df=sims,v=0.99,given=c(1),margin="Normal",method="sequential2", a=NULL, keef_constraints=0) {
  lika <- likb <- likmusig <- a_hat <- b_hat <- bmax <- mu_hat <- sig_hat <- res_var <- c()
  names(df) <- paste0("Y",1:ncol(df))
  d <- ncol(df)
  for (j in given) {
    Y_given1extreme <- df %>% filter(df[,j]>quantile(df[,j],v))
    res <- c(1:d)[-j]
    init_par <- c()
    init_lik <- c()
    for (i in 2:d) {
      # optimise using the initial parameters
      Y1 <- Y_given1extreme[,j]
      Y2 <- Y_given1extreme[,res[i-1]]
      if (method=="sequential2") {
        init_para <- c(0.8,0,1)
        opta <- optim(par=init_para,fn = Y_NLL,df=Y_given1extreme,given=j,sim=res[i-1],b_hat=0,control = list(maxit=2000))
        a_hat <- append(a_hat,opta$par[1])
        lika <- append(lika,opta$value)
        if (1 %in% keef_constraints) {
          b_max1 <- optim(par=0.8,fn = keef_constraint1,a=a_hat[length(a_hat)],Y1=Y1,Y2=Y2,control = list(maxit=2000),lower=0,upper=1, method = "Brent")$par
        } else {b_max1 <- 1}
        if (2 %in% keef_constraints) {
          b_max2 <- optim(par=0.8,fn = keef_constraint2,a=a_hat[length(a_hat)],Y1=Y1,Y2=Y2,control = list(maxit=2000),lower=0, upper=1, method = "Brent")$par
        } else {b_max2 <- 1} 
        b_max <- min(b_max1,b_max2)
        bmax <- append(bmax,b_max)
        init_parb <- c(b_max/2,0,1)
        optb <- optim(par=init_parb,fn = Y_NLL,df=Y_given1extreme,given=j,sim=res[i-1],a_hat=opta$par[1],b_max=b_max,control = list(maxit=2000), method = "Nelder-Mead")
        b_hat <- append(b_hat,optb$par[length(optb$par)-2])
   #     optmusig <- optim(par=init_parb,fn = Y_NLL,df=Y_given1extreme,given=j,sim=res[i-1],a_hat=opta$par[1],b_hat=optb$par[1],control = list(maxit=2000))
        mu_hat <- append(mu_hat,optb$par[length(optb$par)-1])
        sig_hat <- append(sig_hat,optb$par[length(optb$par)])  
        likmusig <- append(likmusig,optb$value)
        likb <- append(likb,optb$value)
      }
      res_var <- append(res_var,res[i-1])
    }
  }
  if (margin=="Normal") {
     if (method %in% c("sequential2")) {
      par_sum <- data.frame("lika" = lika,"likb"=likb,"likmusig"=likmusig,
                            "a" = a_hat, "b" = b_hat,"b_max"=bmax,
                            "mu" = mu_hat,
                            "sig" = sig_hat,
                            "given" = rep(given,each=(d-1)), "res" = res_var)  }
  }

  
  return(par_sum)
}

# function which speeds up and vectorises qspliced ------------------------
qspliced_v_fast <- function(p, x, q, gpd_par) {  
  #' Vectorised quantile function for GPD spliced with empirical
  #' 
  #' @author Wanchen Yue
  #' @param p vector of probabilities between [0,1] for which to find the corresponding quantile
  #' @param x vector of sample used to fit the distribution
  #' @param q float between [0,1], exceedance quantile
  #' @param gpd_par vector of gpd scale and shape parameters
  #' 
  #' @returns vector of quantiles, estimated pth quantile for data x with gpd tail
  #' @export
  
  # Precompute the threshold   
  u <- quantile(x, q, names = FALSE)      
  # Allocate output vector   out <- numeric(length(p))      
  # Split indices based on threshold q   
  below_q <- p <= q   
  above_q <- !below_q      
  # Empirical quantiles for p <= q   
  out <- c()
  out[below_q] <- quantile(x, p[below_q], names = FALSE)      
  # GPD quantiles for p > q   
  if (any(above_q)) {     
    gpd_prob <- 1 - (1 - p[above_q]) / (1 - q)     
    gpd_quant <- evd::qgpd(gpd_prob, scale = gpd_par[1], shape = gpd_par[2])     
    out[above_q] <- u + gpd_quant
  }      
  return(out) 
}

# estimate residual margins
res_margin_par_est <- function(obs_res,method="AGG") {
  if (method=="AGG") {
    tmp <- data.frame("likres"=numeric(),"mu_agg"=numeric(),"sigl"=numeric(),"sigu"=numeric(),"deltal"=numeric(),deltau=numeric())
  }  
  for (i in 1:ncol(obs_res)) {
    Z2 <- as.numeric(unlist(obs_res[,i]))
    if (method=="AGG") {
      opt <- optim(fn=NLL_AGG,x=Z2,par=c(mean(Z2),sd(Z2),sd(Z2),1.2,1.8),control=list(maxit=2000),method = "Nelder-Mead")
      tmp[nrow(tmp)+1,] <- c(opt$value,opt$par)
    }
  }
  return(tmp)
}

