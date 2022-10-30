generate_fixed_parameters <- function() {
#============= Parameter Sampling ===========
nparams = 8

#sigma parameters
sigma_mu = c(1.4, 1.2, 1, 0.8)
sigma_sd = rep(0.2, 4)
sigma_bounds = c(0.3, 5)

#boundary parameters
boundaries_mu = c(2.7, 1.9, 1.2, 0.3)
boundaries_sd = rep(0.2, 4)
boundaries_bounds = c(-10, 10) #limits on boundary position

repeat {
  # sampling sigma first
  sigmas <- (rnorm(length(sigma_mu), 0, 1)*sigma_sd)+sigma_mu
  
  #if sampling sigma - need them to be ordered consecutively
  if ((sum(diff(sigmas) < 0) == (length(sigmas)-1)) & #order check
      all(sigmas > sigma_bounds[1]) & #lower bound check
      all(sigmas < sigma_bounds[2])) { #upper bound check
    break
  }
}

#set a timer because sometimes this will run on an infinite loop
start_time = Sys.time()
time_out = FALSE

repeat {    
  #if sigmas are consecutive - try to find some boundary parameters
  #sample some randomly and see if they meet conditions 
  bounds <- (rnorm(length(boundaries_mu), 0, 1)*boundaries_sd)+boundaries_mu
  
  run_time = Sys.time()
  if (all(bounds > boundaries_bounds[1]) & #lower bound check
      all(bounds > boundaries_bounds[1]) &
      all(bounds > boundaries_bounds[1]) &
      all(bounds > boundaries_bounds[1]) &
      all(bounds < boundaries_bounds[2]) & #upper bound check
      all(bounds < boundaries_bounds[2]) &
      all(bounds < boundaries_bounds[2]) &
      all(bounds < boundaries_bounds[2]) &
      (sum(diff((bounds)) < 0) == 3) & #order check
      (sum(diff((bounds)) < 0) == 3) &
      (sum(diff((bounds)) < 0) == 3) &
      (sum(diff((bounds)) < 0) == 3)) {
    break
  }
  
  if (difftime(run_time, start_time, units = "mins") > 30) {
    time_out = TRUE
    break
  }
}

if (time_out == TRUE) {
  generating_parameters <- rep(0, nparams)
} else {
#use resulting parameters to generate some simulated data
generating_parameters <- c(bounds, sigmas)
}
return(generating_parameters)
}