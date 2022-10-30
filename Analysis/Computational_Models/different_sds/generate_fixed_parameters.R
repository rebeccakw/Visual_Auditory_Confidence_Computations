generate_fixed_parameters <- function() {
#============================
#Rebecca West 2022
#This function randomly generates parameters for the unscaled evidence strength model class.
#Models in this class: Distance 
#Parameters include: ks + ms (boundaries), sigmas
#Parameters are sampled from a normal distribution with a defined mean (referred to with mu) 
#and standard deviation (referred to as sd). All parameters also have a lower and upper bound. 
  
#Function continues to resample parameters from the define distribution until all parameters meet 
#the requirements of that model
#============= Define Distributions to Sample Parameters ===========
nparams = 11 #number of parameters for this model 

#sigma parameters
sigma_mu = c(1.37, 1.03, 0.69, 0.34)
sigma_sd = rep(1, 4)
sigma_bounds = c(0, 5)

#boundary parameters
boundaries_mu = c(0.11, 0.23, 0.54, 1.32, 1.70, 2.34, 2.71)
boundaries_sd = rep(1, 7)
boundaries_bounds = c(0, 10) #limits on boundary position

#=========== Sample Sigma Parameters ===============
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

#=========== Sample Boundary Parameters ===============
#set a timer because sometimes this will run on an infinite loop
start_time = Sys.time()
time_out = FALSE

repeat {    
  #if sigmas are consecutive - try to find some boundary parameters
  #sample some randomly and see if they meet conditions 
  bounds <- (rnorm(length(boundaries_mu), 0, 1)*boundaries_sd)+boundaries_mu
  
  run_time = Sys.time()
  if (all(bounds > boundaries_bounds[1]) & #lower bound check
      all(bounds < boundaries_bounds[2]) & #upper bound check
      all(sort(bounds) == bounds)) {
    break
  }
  
  if (difftime(run_time, start_time, units = "mins") > 30) {
    time_out = TRUE
    break
  }
}

if (time_out == TRUE) {
  generating_parameters <- rep(NA, nparams)
} else {
#use resulting parameters to generate some simulated data
generating_parameters <- c(bounds, sigmas)
}
return(generating_parameters)
}