################################################################################
## Distribution Function, Quantile Function, and Random number generator for non-stationary GPD
pspliced_nonsta <- function(x, t, gpd_par) {
  # x: sample used for empirical CDF 
  # gpd_par: list with:
  #    - scale: vector for x0 > u only
  #    - shape: single shared value
  stopifnot(length(x) == length(t),
            length(gpd_par$u) == 3,
            length(gpd_par$scale) == 3,
            length(gpd_par$shape) == 1)
  
  
  n_data <- length(x)
  #-----------------------------------------------------------------------------
  ## GPD parameters
  
  ## Obtain the threshold
  u <- Model_Year(t = t, 
                  beta_1 = gpd_par$u[1],
                  beta_2 = gpd_par$u[2],
                  beta_3 = gpd_par$u[3])
  
  ## Obtain the scale parameter
  scale <- exp(Model_Year(t = t, 
                          beta_1 = gpd_par$scale[1],
                          beta_2 = gpd_par$scale[2],
                          beta_3 = gpd_par$scale[3]))
  
  #-----------------------------------------------------------------------------
  
  ## Empirical CDF
  F_all <- ecdf(x)
  F_body <- ecdf(x[x <=  u])
  
  ## GPD component
  phi_t <- 1 - F_all(u)
  cdf_vals <- (1 - phi_t) + phi_t*(1 - pmax(0, 1 + gpd_par$shape*(x - u)/scale)^(-1/gpd_par$shape))
  
  ## Empirical part  
  below_idx <- which(x <= u)
  if (length(below_idx) > 0){
    constant <- (1 - phi_t[below_idx])/F_body(u[below_idx])
    cdf_vals[below_idx] <- constant*F_body(x[below_idx])
  }
  
  return(cdf_vals)
}

qspliced_nonsta <- function(p, x, t, gpd_par) {
  # p: vector of probabilities
  # x: empirical sample
  # gpd_par:
  #   - scale: vector, same length as sum(p > q)
  #   - shape: scalar
  
  stopifnot(length(p) == length(x),
            length(p) == length(t),
            length(gpd_par$u) == 3,
            length(gpd_par$scale) == 3,
            length(gpd_par$shape) == 1)
  
  ## Number of probabilities to evaluate
  n_data <- length(p)
  
  #-----------------------------------------------------------------------------
  ## GPD parameters
  
  ## Obtain the threshold
  u <- Model_Year(t = t, 
                  beta_1 = gpd_par$u[1],
                  beta_2 = gpd_par$u[2],
                  beta_3 = gpd_par$u[3])
  
  ## Obtain the scale parameter
  scale <- exp(Model_Year(t = t, 
                          beta_1 = gpd_par$scale[1],
                          beta_2 = gpd_par$scale[2],
                          beta_3 = gpd_par$scale[3]))
  
  #-----------------------------------------------------------------------------
  
  ## Empirical CDFs
  F_all <- ecdf(x)
  F_body <- ecdf(x[x <= u])
  
  #-----------------------------------------------------------------------------
  
  ## Obtain the tail component 
  u_unif <- F_all(u)
  phi_t <- 1 - u_unif
  
  ## Obtain the tail part
  gpd_prob <- (p - (1 - phi_t))/phi_t
  out <- u + scale*((1 - gpd_prob)^(-gpd_par$shape) - 1)/gpd_par$shape
  
  ## Obtain the body
  below_idx <- p <= u_unif
  if (any(below_idx)) {
    p_body <- p[below_idx]*F_body(u[below_idx])/(1 - phi_t[below_idx])
    out[below_idx] <- quantile(x[x <= u], p_body, type = 1, names = FALSE)
  }
  
  return(out)
}

rspliced_nonsta <- function(x, gpd_par) {
  # x: empirical sample
  # gpd_par:
  #   - scale: vector, same length as sum(p > q)
  #   - shape: scalar
  
  stopifnot(is.numeric(x),
            length(gpd_par$u) == 3,
            length(gpd_par$scale) == 3,
            length(gpd_par$shape) == 1)
  
  ## Probabilities to evaluate
  p <- runif(length(x))
  
  ## Obtain the random sample
  out <- qspliced_nonsta(p = p,
                         x = x,
                         t = 1:lenghth(x),
                         gpd_par = gpd_par)
  return(out)
}