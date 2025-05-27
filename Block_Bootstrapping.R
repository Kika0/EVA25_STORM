## Combine the runs into a single data frame
d <- dim(run1)[1]
n_grid <- d*d
n_data <- dim(run1)[3]
n_runs <- 4

grid_locs <- permutations(n = d, r = 2, v = 1:d, repeats.allowed = TRUE)
n_years <- 165
n_months <- 12
n_days <- 365
n_seasons <- 4

days_per_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
days_per_season <- c(sum(days_per_month[1:3]),
                     sum(days_per_month[4:6]),
                     sum(days_per_month[7:9]),
                     sum(days_per_month[10:12]))

data <- data.frame(Run = rep(1:n_runs, each = n_days*n_years*n_grid),
                   Grid_1 = rep(rep(grid_locs[,1], each = n_days*n_years), n_runs),
                   Grid_2 = rep(rep(grid_locs[,2], each = n_days*n_years), n_runs),
                   Time = rep(1:(n_days*n_years), times = n_grid),
                   Day = rep(rep(rep(do.call(c, sapply(days_per_month, function(i){1:i})), times = n_years), times = n_grid), n_runs),
                   Month = rep(rep(rep(do.call(c, sapply(1:n_months, function(i){rep(i, times = days_per_month[i])})), times = n_years), times = n_grid), n_runs),
                   Year = rep(rep(rep(1:n_years, each = n_days), times = n_grid), n_runs),
                   Season = rep(rep(rep(do.call(c, sapply(1:n_seasons, function(i){rep(i, times = days_per_season[i])})), times = n_years), times = n_grid), n_runs),
                   Value = c(c(apply(grid_locs, 1, function(x){c(run1[x[1],x[2],])})),
                             c(apply(grid_locs, 1, function(x){c(run2[x[1],x[2],])})),
                             c(apply(grid_locs, 1, function(x){c(run3[x[1],x[2],])})),
                             c(apply(grid_locs, 1, function(x){c(run4[x[1],x[2],])}))))

data$Run <- factor(data$Run, levels = 1:n_runs, labels = 1:n_runs)
data$Grid_1 <- factor(data$Grid_1, levels = 1:d, labels = 1:d)
data$Grid_2 <- factor(data$Grid_2, levels = 1:d, labels = 1:d)
# data$Value <- as.numeric(data$Value)
################################################################################
## Function to block bootstrap from the non-stationary GPD
NonSta_GPD_to_Lapalce <- function(df, run, tau, B, l){
  
  ## Filter the data to the run of choice
  data <- filter(df, Run == run)
  
  ## Obtain the marginal fits
  grid_list <- lapply(apply(expand.grid(Grid_2 = 1:d, Grid_1 = 1:d, Run = run), 1, list),
                      function(x){x[[1]]})
  Marginal_Fits <- map(.x = grid_list,
                       .f = function(x){
                         Fit_GPD_Seasonal_Threshold_Scale(df = filter(data,
                                                                      Run == x[3],
                                                                      Grid_1 == x[2],
                                                                      Grid_2 == x[1]),
                                                          tau = tau)})
  
  ## Obtain the data on uniform margins
  Value_unif <- lapply(seq_along(grid_list), function(i){
    pspliced_nonsta(x = unname(unlist(select(filter(data,
                                                    Run == grid_list[[i]][3],
                                                    Grid_1 == grid_list[[i]][2],
                                                    Grid_2 == grid_list[[i]][1]),
                                             Value))),
                    t = unname(unlist(select(filter(data,
                                                    Run == grid_list[[i]][3],
                                                    Grid_1 == grid_list[[i]][2],
                                                    Grid_2 == grid_list[[i]][1]),
                                             Time))),
                    gpd_par = list(u = Marginal_Fits[[i]]$par$beta_u,
                                   scale = Marginal_Fits[[i]]$par$beta_scale,
                                   shape = Marginal_Fits[[i]]$par$shape))
  })
  data$Value_unif <- do.call(c, Value_unif)
  
  ## Convert the data to wide format to allow for easy bootstrapping
  data <- data %>%
    mutate(ID = paste0("Run", Run, "_G", Grid_1, "_", Grid_2))
  data_wide <- data %>%
    select(Time, ID, Value_unif) %>%
    pivot_wider(names_from = ID, values_from = Value_unif)
  
  
  ## Function to bootstrap the above sample
  Boot_Refit <- function(n_data, l, data_wide, data, grid_list){
    
    ## Create the blocks to bootstrap the data
    ## Upper bound on the block to not exceed the time series length
    blocks <- c()
    while(sum(blocks) < n_data){
      blocks <- c(blocks, min(n_data - sum(blocks), rgeom(n = 1, prob = 1/l)))
    }
    blocks <- blocks[blocks > 0]
    n_blocks <- length(blocks)
    
    ## Obtain the time stamps with those blocks
    ## Upper bound the time so that block does not exceeed time series length
    time_samples <- sapply(blocks, function(x){sample(1:(n_data - x), 1)})
    time_samples_end <- time_samples + blocks - 1
    
    ## Order the time in ascending order 
    ## (not necessary for dependence model but) good for diagnostics
    ordering_time <- order(time_samples)
    time_samples <- time_samples[ordering_time]
    time_samples_end <- time_samples_end[ordering_time]
    
    ## bootstrap the data
    data_boot <- do.call(rbind, pmap(.f = function(x, y){data_wide[x:y,]},
                                     .l = list(x = time_samples, y = time_samples_end)))
    
    ## Convert back to long format
    data_boot_long <- data_boot %>%
      pivot_longer(
        cols = matches("^Run\\d+_G\\d+_\\d+$"),  # Select all the value columns
        names_to = c("Run", "Grid_1", "Grid_2"),
        names_pattern = "Run(\\d+)_G(\\d+)_(\\d+)",
        values_to = "Value_Unif"
      ) %>%
      mutate(across(c(Run, Grid_1, Grid_2), as.integer))
    
    data_boot_long <- arrange(data_boot_long, Run, Grid_1, Grid_2)
    
    ## Transform the data back onto the original margins
    data_boot_long$Value <- c(sapply(seq_along(grid_list), function(i){
      qspliced_nonsta(p = unname(unlist(select(filter(data_boot_long,
                                                      Run == grid_list[[i]][3],
                                                      Grid_1 == grid_list[[i]][2],
                                                      Grid_2 == grid_list[[i]][1]),
                                               Value_Unif))),
                      t = unname(unlist(select(filter(data_boot_long,
                                                      Run == grid_list[[i]][3],
                                                      Grid_1 == grid_list[[i]][2],
                                                      Grid_2 == grid_list[[i]][1]),
                                               Time))),
                      x = unname(unlist(select(filter(data,
                                                      Run == grid_list[[i]][3],
                                                      Grid_1 == grid_list[[i]][2],
                                                      Grid_2 == grid_list[[i]][1]),
                                               Value))),
                      gpd_par = list(u = Marginal_Fits[[i]]$par$beta_u,
                                     scale = Marginal_Fits[[i]]$par$beta_scale,
                                     shape = Marginal_Fits[[i]]$par$shape))
    }))
    
    ## Refit the model
    Marginal_Fits_Boot <- map(.x = grid_list,
                              .f = function(x){
                                Fit_GPD_Seasonal_Threshold_Scale(df = filter(data_boot_long,
                                                                             Run == x[3],
                                                                             Grid_1 == x[2],
                                                                             Grid_2 == x[1]),
                                                                 tau = tau)})
    
    ## Transform onto uniform margins
    U_Boot <- sapply(seq_along(grid_list), function(i){
      pspliced_nonsta(x = unname(unlist(select(filter(data_boot_long,
                                                      Run == grid_list[[i]][3],
                                                      Grid_1 == grid_list[[i]][2],
                                                      Grid_2 == grid_list[[i]][1]),
                                               Value))),
                      t = unname(unlist(select(filter(data_boot_long,
                                                      Run == grid_list[[i]][3],
                                                      Grid_1 == grid_list[[i]][2],
                                                      Grid_2 == grid_list[[i]][1]),
                                               Time))),
                      gpd_par = list(u = Marginal_Fits_Boot[[i]]$par$beta_u,
                                     scale = Marginal_Fits_Boot[[i]]$par$beta_scale,
                                     shape = Marginal_Fits_Boot[[i]]$par$shape))
    })
    
    ## Transform onto standard Laplace margins
    Y_Boot <- apply(U_Boot, 2, qlaplace)
    
    Par_Boot <- t(sapply(Marginal_Fits_Boot, function(x){
      c(x$par$beta_u, x$par$beta_scale, x$par$shape)})) 
    
    return(list(X = matrix(data_boot_long$Value, ncol = length(grid_list)),
                Y = Y_Boot,
                t = data_boot$Time,
                par = Par_Boot))
  }
  
  ## Obtain the B bootstraps
  Boot_Out <- replicate(n = B,
                        expr = Boot_Refit(n_data = nrow(data_wide),
                                          l = l,
                                          data_wide = data_wide,
                                          data = data,
                                          grid_list = grid_list),
                        simplify = FALSE)
  return(Boot_Out)
}

start_time <- Sys.time()
Y_Bootstrapped <- NonSta_GPD_to_Lapalce(df = data,
                                        run = 1,
                                        tau = 0.98,
                                        B = 1,
                                        l = 2)
end_time <- Sys.time()
end_time - start_time


