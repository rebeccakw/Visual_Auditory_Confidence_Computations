generate_exp_parameters <- function(confidencedata) {
  #============= Parameter Sampling ===========
  #required for ODN 
  confidencedata$unscaled_stim_orientation = confidencedata$stim_orientation
  
  sample_mean = mean(confidencedata$stim_orientation)
  sample_sd = sd(confidencedata$stim_orientation)
  #use this for odn noise 
  z_score <- function(value) {
    z = (value - sample_mean)/sample_sd
    return(z)
  }
  
  boundary_function <- function(ks, ms, sig, exp) {
    bound = ks + (ms*(sig^exp))
    return(bound)
  }
  
  confidencedata$stim_orientation = scale(confidencedata$stim_orientation)
  
  nparams = 14
  
  #sigma parameters
  sigma_mu = c(2.84,   1.85,   1.25,    0.8)
  #sigma_mu = c(3,	2,	1,	0.8)
  sigma_sd = rep(1, 4)
  sigma_bounds = c(0.1, 10)
  
  #boundary parameters
  kboundaries_mu = c(1.27, 0.61, 0.17, -0.04)
  #kboundaries_mu = c(1.078, 0.003, 0.002, 0.025)
  kboundaries_sd = rep(2, 4)
  
  mboundaries_mu = c(0.73, 0.7, 0.46, 0.16)
  #mboundaries_mu = c(0.959, 0.726, 0.622, 0.078)
  mboundaries_sd = rep(2, 4)
  
  exp_mu = 2
  exp_sd = 1
  
  odn_mu = 1
  odn_sd = 0.5
  odn_bounds = c(0, 20)
  
  boundary <- function(k, m, sigma, exp) {
    b = k + (m*(sigma^exp))
    return(b)
  }
  
  zscore <- function(value){
  sample_mean = mean(confidencedata$stim_orientation)
  sample_sd = sd(confidencedata$stim_orientation)
  z = (value-sample_mean)/sample_sd
  return(z)
  }
  
  repeat {
    # sampling sigma first
    sigmas <- (rnorm(length(sigma_mu), 0, 1)*sigma_sd)+sigma_mu
    exp <- rnorm(1, exp_mu, exp_sd)
    odn <- rnorm(1, odn_mu, odn_sd)
    
    
    #if sampling sigma - need them to be ordered consecutively
    if ((sum(diff(sigmas) <= 0) == (length(sigmas)-1)) & #order check
        all(sigmas > sigma_bounds[1]) & #lower bound check
        all(sigmas < sigma_bounds[2]) & 
        exp < 10 & 
        exp > 1 & 
        (odn > odn_bounds[1]) & 
        (odn < odn_bounds[2])) { #upper bound check
      break
    }
  }
  
  #set a timer because sometimes this will run on an infinite loop
  start_time = Sys.time()
  time_out = FALSE
  
  baseline_sigma = sigmas[confidencedata$stim_reliability_level] 
  all_odn = z_score(odn*(abs(sin((confidencedata$unscaled_stim_orientation*pi)/90))))
  all_sigma = baseline_sigma + all_odn
  
  repeat {    
    #if sigmas are consecutive - try to find some boundary parameters
    #sample some randomly and see if they meet conditions 
    kbounds <- (rnorm(length(kboundaries_mu), 0, 1)*kboundaries_sd)+kboundaries_mu
    mbounds <- (rnorm(length(mboundaries_mu), 0, 1)*mboundaries_sd)+mboundaries_mu
    
    boundaries = t(sapply(all_sigma, function(x) boundary_function(kbounds, mbounds, x, exp)))
    all_boundaries = cbind(-boundaries[,1:3], #make first 3 boundaries negative
                           boundaries[,4], #category boundary
                           t(apply(boundaries[,1:3], 1, function(x) rev(x)))) #positive boundaries
    run_time = Sys.time()
    if (all(t(apply(all_boundaries, 1, function(x) sort(x))) == all_boundaries)) {
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
    generating_parameters <- c(kbounds, mbounds, sigmas, exp, odn)
  }
  return(generating_parameters)
}