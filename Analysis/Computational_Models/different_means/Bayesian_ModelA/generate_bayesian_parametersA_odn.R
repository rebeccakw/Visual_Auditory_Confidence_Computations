generate_bayesian_parameters_odn <- function(modality = 'visual') {
  #============= Parameter Sampling ===========
  nparams = 9
  
  #sigma parameters
  sigma_mu = c(13,  7.5,  4.5,  2.5) #from group fits 
  #sigma_mu = c(14.50,  8.10,  7,  5)
  sigma_sd = rep(1, 4) #4
  sigma_bounds = c(1, 50)
  
  #boundary parameters
  boundaries_mu = c(-0.2, -0.6, -1, -1.5)
  #boundaries_mu = c(-1.5, -1, -0.6, -0.2) #from group fits 
  boundaries_sd = rep(2, 4) #3 
  boundaries_bounds = c(-60, 60) #limits on boundary position in orientation space
  
  odn_mu = 5 #from group fits 
  #odn_mu = 2
  odn_sd = 2
  odn_bounds = c(0, 20)
  
  if (modality == 'visual') {
    cat_mean = -4
    cat_sigma = 5
  } 
  if (modality == 'auditory') {
    cat_mean = 2300
    cat_sigma = 475
  }
  
  #write function to convert from d space to measurement space for each boundary 
  measurement <- function(d, sigma)  {
    measurement = (d*((sigma^2)+(cat_sigma^2)))/(2*cat_mean)
    return(measurement)  
  }
  
  start_time = Sys.time()
  time_out = FALSE
  
  repeat {
    # sampling sigma first
    sigmas <- (rnorm(length(sigma_mu), 0, 1)*sigma_sd)+sigma_mu
    odn <- rnorm(1, odn_mu, odn_sd)
    
    #if sampling sigma - need them to be ordered consecutively
    if (all(sigmas == sort(sigmas, decreasing = TRUE)) & #order check
        all(sigmas > sigma_bounds[1]) & #lower bound check
        all(sigmas < sigma_bounds[2]) & 
        (odn > odn_bounds[1]) & 
        (odn < odn_bounds[2])) { #upper bound check
      break
    }
    
    run_time = Sys.time()
    if (difftime(run_time, start_time, units = "mins") > 5) {
      time_out = TRUE
      break
    }
  }
  
  #set a timer because sometimes this will run on an infinite loop
  if (time_out == FALSE) {
  start_time = Sys.time()
  min_ori = 0
  max_ori = 50
  noise_min = sigmas + (odn*(abs(sin((min_ori*pi)/90))))
  noise_max = sigmas + (odn*(abs(sin((max_ori*pi)/90))))
  
  repeat {    
    #if sigmas are consecutive - try to find some boundary parameters
    #sample some randomly and see if they meet conditions 
    bounds <- (rnorm(length(boundaries_mu), 0, 1)*boundaries_sd)+boundaries_mu
    bounds <- cumsum(bounds)
    
    rel1_min = c(-measurement(bounds[1:3], noise_min[1]), measurement(bounds[4], noise_min[1]), rev(measurement(bounds[1:3], noise_min[1])))
    rel2_min = c(-measurement(bounds[1:3], noise_min[2]), measurement(bounds[4], noise_min[2]), rev(measurement(bounds[1:3], noise_min[2])))
    rel3_min = c(-measurement(bounds[1:3], noise_min[3]), measurement(bounds[4], noise_min[3]), rev(measurement(bounds[1:3], noise_min[3])))
    rel4_min = c(-measurement(bounds[1:3], noise_min[4]), measurement(bounds[4], noise_min[4]), rev(measurement(bounds[1:3], noise_min[4])))
    rel1_max = c(-measurement(bounds[1:3], noise_max[1]), measurement(bounds[4], noise_max[1]), rev(measurement(bounds[1:3], noise_max[1])))
    rel2_max = c(-measurement(bounds[1:3], noise_max[2]), measurement(bounds[4], noise_max[2]), rev(measurement(bounds[1:3], noise_max[2])))
    rel3_max = c(-measurement(bounds[1:3], noise_max[3]), measurement(bounds[4], noise_max[3]), rev(measurement(bounds[1:3], noise_max[3])))
    rel4_max = c(-measurement(bounds[1:3], noise_max[4]), measurement(bounds[4], noise_max[4]), rev(measurement(bounds[1:3], noise_max[4])))
    boundaries_all = rbind(rel1_min, rel2_min, rel3_min, rel4_min,
                           rel1_max, rel2_max, rel3_max, rel4_max) 
    
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
  }

  if (time_out == TRUE) {
    generating_parameters <- rep(0, nparams)
  } else {
    #use resulting parameters to generate some simulated data
    generating_parameters <- c(bounds, sigmas, odn)
  }
  return(generating_parameters)
}