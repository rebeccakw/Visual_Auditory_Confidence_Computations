generate_bayesian_parameters_odn <- function(mus, sigmas, confidencedata, modality = 'visual') {
  #mus input 
  #1-7: boundaries
  #8-11: sigma4 (largest), sigma3, sigma2, sigma1 (smallest)
  #12: psi (if it exists)
  #13: dnoise
  #sigma_mu = c(4,	3,	2,	1)
  #sigma_sd = rep(6, 4)
  #boundaries_mu = c(0.5, -0.1, -0.2, -0.6, -1.0, -1.6, -1.7)
  #boundaries_sd = rep(4, 7)
  #d_noise_mu = 0.24
  #d_noise_sd = 1
  
  #================================
  # Set Distribution of Parameters To Be Sampled
  #================================
  #boundary parameters
  boundaries_mu = mus[1:7]
  boundaries_sd = sigmas[1:7]
  boundaries_bounds = c(0, 20)
  
  #sigma parameters
  sigma_mu = mus[8:11]
  sigma_sd =  sigmas[8:11]
  sigma_bounds = c(0, 10)
  
  odn_mu =  mus[12]
  odn_sd = 0.5
  odn_bounds = c(0.2, 20)
  
  if (modality == 'visual') {
    cat_one_sigma = 0.3427148
    cat_two_sigma = 1.370859
  } else if (modality == 'auditory') {
    cat_one_sigma = 0.3429777
    cat_two_sigma = 1.371911
  }
  
  #write function to convert from d space to measurement space for each boundary 
  d_to_ori <- function(d, s) {
    a <- (1/2*log(( (s^2) + (cat_two_sigma^2) )/( (s^2) + (cat_one_sigma^2) )))
    b <- ((( (cat_two_sigma^2) - (cat_one_sigma^2) )/( 2*( (s^2) + (cat_one_sigma^2) ) *( (s^2) + (cat_two_sigma^2) ) )))
    prod <- (a - d)/b
    
    #find where values are less than 0 
    index <- (prod < 0)
    
    #take the absolute value of negative values, so we can find the square root
    prod[index] <- abs(prod[index])
    prod <- sqrt(prod)
    
    #returns the original negative values as negative values 
    prod[index] <- - (prod[index])
    return(prod)  
  }
  
  #required for ODN 
  confidencedata$unscaled_stim_orientation = confidencedata$stim_orientation
  
  sample_mean = mean(confidencedata$stim_orientation)
  sample_sd = sd(confidencedata$stim_orientation)
  
  #use this for odn noise 
  z_score <- function(value) {
    z = (value - sample_mean)/sample_sd
    return(z)
  }
  confidencedata$stim_orientation = scale(confidencedata$stim_orientation)
  
  #set a timer because sometimes this will run on an infinite loop
  start_time = Sys.time()
  time_out = FALSE
  
  #================================
  # Sample Sigma Parameters
  #================================
    boundary_draws <-  matrix(0, 7, 4)
    repeat {
      # sampling sigma first
      sigmas <- (rnorm(length(sigma_mu), 0, 1)*sigma_sd)+sigma_mu
      bounds <- (rnorm(length(boundaries_mu), 0, 1)*boundaries_sd)+boundaries_mu
      boundaries = cumsum(bounds[1:7])
      odn <- rnorm(1, odn_mu, odn_sd)
     
      
      #check if sampled sigma parameters meet conditions
      if ((sum(diff(sigmas) < 0) == (length(sigmas)-1)) & #order check
          all(sigmas > sigma_bounds[1]) & #lower bound check
          all(sigmas < sigma_bounds[2]) & #upper bound check 
          all(diff(round(sigmas, 1)) < 0) &
          (odn > odn_bounds[1]) & 
          (odn < odn_bounds[2])) { 
        
        #if sampled sigma parameters meet condition, check boundary parameters
        #have to check sigmas first because if sigma is less than 0, this will throw an error
        baseline_sigma = sigmas[confidencedata$stim_reliability_level] 
        all_odn = z_score(odn*(abs(sin((confidencedata$unscaled_stim_orientation*pi)/90))))
        all_sigma = baseline_sigma + all_odn
        
        #get all possible boundary draws for each intensity level
        all_boundaries = t(sapply(all_sigma, function(x) d_to_ori(boundaries, x)))

        #check that boundary draws meet conditions
        if (all(all_boundaries > boundaries_bounds[1]) & #lower bound check
            all(all_boundaries < boundaries_bounds[2]) & #upper bound check - probably not necessary
            all(t(apply(all_boundaries, 1, function(x) sort(x))) == all_boundaries)) { #check that boundaries are monotonically increasing 
          break   #stop sampling if conditions have been met  
        }
      }
      
      #if sampling is taking too long (which it shouldn't), stop the sampling
      run_time = Sys.time()
      if (difftime(run_time, start_time, units = "mins") > 15) {
        time_out = TRUE
        break
      }
      
    }

#just so this function doesn't throw an error, if sampling has stopped
if (time_out == TRUE) {
  generating_parameters <- rep(NA, length(mus))
} else {
  #return parameters
generating_parameters <- c(bounds, sigmas, odn)
}
return(generating_parameters)
}
