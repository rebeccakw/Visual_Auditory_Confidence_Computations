generate_exp_parameters <- function() {
  #============= Parameter Sampling ===========
  nparams = 14
  
  #sigma parameters
  sigma_mu = c(2.84,   1.85,   1.25,    0.8)
  #sigma_mu = c(3,	2,	1,	0.8)
  sigma_sd = rep(1, 4)
  sigma_bounds = c(0.1, 10)
  
  #boundary parameters
  kboundaries_mu = c(1.27, 0.61, 0.17, -0.04)
  #kboundaries_mu = c(1.078, 0.003, 0.002, 0.025)
  kboundaries_sd = rep(1, 4)
  
  mboundaries_mu = c(0.73, 0.7, 0.46, 0.16)
  #mboundaries_mu = c(0.959, 0.726, 0.622, 0.078)
  mboundaries_sd = rep(1, 4)
  
  exp_mu = 2
  exp_sd = 1
  
  scale_mu = 1.2
  scale_sd = 1
  
  boundaries_bounds = c(0, 30) #limits on boundary position
  
  boundary <- function(k, m, sigma, exp) {
    b = k + (m*(sigma^exp))
    return(b)
  }
  
  repeat {
    # sampling sigma first
    sigmas <- (rnorm(length(sigma_mu), 0, 1)*sigma_sd)+sigma_mu
    exp <- rnorm(1, exp_mu, exp_sd)
    scaler <- rnorm(1, scale_mu, scale_sd)
    
    #if sampling sigma - need them to be ordered consecutively
    if ((sum(diff(sigmas) <= 0) == (length(sigmas)-1)) & #order check
        all(sigmas > sigma_bounds[1]) & #lower bound check
        all(sigmas < sigma_bounds[2]) & 
        exp < 10 & 
        exp > 1 & 
        scaler < 3 & 
        scaler > 0) { #upper bound check
      break
    }
  }
  
  #set a timer because sometimes this will run on an infinite loop
  start_time = Sys.time()
  time_out = FALSE
  boundaries_all = matrix(0, 4, 8)
  
  repeat {    
    #if sigmas are consecutive - try to find some boundary parameters
    #sample some randomly and see if they meet conditions 
    kbounds <- (rnorm(length(kboundaries_mu), 0, 1)*kboundaries_sd)+kboundaries_mu
    mbounds <- (rnorm(length(mboundaries_mu), 0, 1)*mboundaries_sd)+mboundaries_mu
    sigmas_all = c(sigmas, (scaler*sigmas))
                   
    for (s in 1:length(sigmas_all)) {
      boundaries_all[,s] <- boundary(kbounds, mbounds, sigmas_all[s], exp)
    }

    run_time = Sys.time()
    if (all(boundaries_all < boundaries_bounds[2]) & #lower bound check
        all(boundaries_all[2:4,] > boundaries_bounds[1]) &
        all(abs(boundaries_all[1,]) < abs(boundaries_all[2,])) &
        all(apply(boundaries_all, 2, sort) == (boundaries_all))) { #order check
      break
    }
    
    if (difftime(run_time, start_time, units = "mins") > 10) {
      time_out = TRUE
      break
    }
  }
  
  if (time_out == TRUE) {
    generating_parameters <- rep(NA, nparams)
  } else {
    #use resulting parameters to generate some simulated data
    generating_parameters <- c(kbounds, mbounds, sigmas, exp, scaler)
  }
  return(generating_parameters)
}