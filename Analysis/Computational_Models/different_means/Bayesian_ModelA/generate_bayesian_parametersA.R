generate_bayesian_parameters <- function(modality = 'visual') {
  #============= Parameter Sampling ===========
  nparams = 8
  
  #sigma parameters
  sigma_mu = c(2.5, 1.5, 1, 0.9)
  sigma_sd = rep(1.5, 4)
  sigma_bounds = c(0.1, 10)
  
  #boundary parameters
  boundaries_mu = c(-1.4, -1, -0.6, -0.2)
  boundaries_sd = rep(1.5, 4)
  boundaries_bounds = c(-30, 30) #limits on boundary position in orientation space
  
  if (modality == 'visual') {
    cat_mean = -0.6247561
    cat_sigma = 0.7809452
  } 
  if (modality == 'auditory') {
    cat_mean = -0.6449836
    cat_sigma = 0.7608733
  }
  
  #write function to convert from d space to measurement space for each boundary 
  measurement <- function(d, sigma)  {
    measurement = (d*((sigma^2)+(cat_sigma^2)))/(2*cat_mean)
    return(measurement)  
  }
  
  repeat {
    # sampling sigma first
    sigmas <- (rnorm(length(sigma_mu), 0, 1)*sigma_sd)+sigma_mu
    
    #if sampling sigma - need them to be ordered consecutively
    if ((sum(diff(sigmas) <= 0) == (length(sigmas)-1)) & #order check
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
    
    rel1 = c(-measurement(bounds[1:3], sigmas[1]), measurement(bounds[4], sigmas[1]), rev(measurement(bounds[1:3], sigmas[1])))
    rel2 = c(-measurement(bounds[1:3], sigmas[2]), measurement(bounds[4], sigmas[2]), rev(measurement(bounds[1:3], sigmas[2])))
    rel3 = c(-measurement(bounds[1:3], sigmas[3]), measurement(bounds[4], sigmas[3]), rev(measurement(bounds[1:3], sigmas[3])))
    rel4 = c(-measurement(bounds[1:3], sigmas[4]), measurement(bounds[4], sigmas[4]), rev(measurement(bounds[1:3], sigmas[4])))
    boundaries_all = rbind(rel1, rel2, rel3, rel4) 
    
    run_time = Sys.time()
    if (all(boundaries_all > boundaries_bounds[1]) & #lower bound check
        all(boundaries_all < boundaries_bounds[2]) & #upper bound check
        all(t(apply(boundaries_all, 1, sort)) == boundaries_all)) { #order check
      break
    }
    
    if (difftime(run_time, start_time, units = "mins") > 10) {
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