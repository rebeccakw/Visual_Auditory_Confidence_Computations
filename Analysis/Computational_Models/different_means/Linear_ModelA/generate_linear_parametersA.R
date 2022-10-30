generate_linear_parameters <- function() {
  #============= Parameter Sampling ===========
  nparams = 12
  
  #sigma parameters
  sigma_mu = c(3,	2,	1,	0.8)
  sigma_sd = rep(0.5, 4)
  sigma_bounds = c(0.1, 10)
  
  #boundary parameters
  kboundaries_mu = c(-1.7,	-2.1,	-1.8,	-0.2)
  kboundaries_sd = rep(0.5, 4)
  
  mboundaries_mu = c(3.7,	2.8, 2.4,	0.3)
  mboundaries_sd = rep(0.5, 4)
  
  boundaries_bounds = c(-30, 30) #limits on boundary position
  
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
    kbounds <- (rnorm(length(kboundaries_mu), 0, 1)*kboundaries_sd)+kboundaries_mu
    mbounds <- (rnorm(length(mboundaries_mu), 0, 1)*mboundaries_sd)+mboundaries_mu
    
    rel1 = c(-(kbounds[1:3] + (mbounds[1:3]*sigmas[1])), (kbounds[4] + mbounds[4]*sigmas[1]), rev(kbounds[1:3] + (mbounds[1:3]*sigmas[1])))
    rel2 = c(-(kbounds[1:3] + (mbounds[1:3]*sigmas[2])), (kbounds[4] + mbounds[4]*sigmas[2]), rev(kbounds[1:3] + (mbounds[1:3]*sigmas[2])))
    rel3 = c(-(kbounds[1:3] + (mbounds[1:3]*sigmas[3])), (kbounds[4] + mbounds[4]*sigmas[3]), rev(kbounds[1:3] + (mbounds[1:3]*sigmas[3])))
    rel4 = c(-(kbounds[1:3] + (mbounds[1:3]*sigmas[4])), (kbounds[4] + mbounds[4]*sigmas[4]), rev(kbounds[1:3] + (mbounds[1:3]*sigmas[4])))
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
    generating_parameters <- c(kbounds, mbounds, sigmas)
  }
  return(generating_parameters)
}