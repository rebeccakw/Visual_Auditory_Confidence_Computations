generate_bayesian_parameters <- function(mus, sigmas, modality = 'visual') {
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
  boundaries_bounds = c(0, 100)
  
  #sigma parameters
  sigma_mu = mus[8:11]
  sigma_sd =  sigmas[8:11]
  sigma_bounds = c(0, 40)
  
  if (length(mus) > 12) {
    psi_mu = mus[12]
    psi_sd = sigmas[12]
    d_noise_mu = mus[13]
    d_noise_sd = sigmas[13]
  } else {
    d_noise_mu = mus[12]
    d_noise_sd = sigmas[12]  
  }
  
  psi_bounds = c(0, 40)
  d_noise_bounds = c(0, 10)
  
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
  
  
  get_draws <- function(d_bound, d_sd) {
    spaced_draws <- seq((qnorm(0.01, mean = d_bound, sd = d_sd)), (qnorm(0.99, mean = d_bound, sd = d_sd)), length.out = 101)
    return(spaced_draws)
  }
  
  odn_function <- function(psi, ori) {
    odn <- psi*(abs(sin((ori*pi)/90)))
    return(odn)
  }
  
  #set a timer because sometimes this will run on an infinite loop
  start_time = Sys.time()
  time_out = FALSE
  
  #================================
  # Sample Sigma Parameters
  #================================
  # doing the sampling seperately when psi is involved (because it's much more annoying)
  if (length(mus) == 12) {
    boundary_draws <-  matrix(0, 101, 28)
    repeat {
      # sampling sigma first
      sigmas <- (rnorm(length(sigma_mu), 0, 1)*sigma_sd)+sigma_mu
      d_noise <- (rnorm(length(d_noise_mu), 0, 1)*d_noise_sd)+d_noise_mu
      bounds <- (rnorm(length(boundaries_mu), 0, 1)*boundaries_sd)+boundaries_mu
      boundaries = cumsum(bounds[1:7])
      
      #check if sampled sigma parameters meet conditions
      if ((sum(diff(sigmas) < 0) == (length(sigmas)-1)) & #order check
          all(sigmas > sigma_bounds[1]) & #lower bound check
          all(sigmas < sigma_bounds[2]) & #upper bound check 
          all(diff(round(sigmas, 1)) < 0) & #make sure that sigmas are at least 1 decimal point different from one another
          d_noise > d_noise_bounds[1] &
          d_noise < d_noise_bounds[2]) { 
        
        #if sampled sigma parameters meet condition, check boundary parameters
        #have to check sigmas first because if sigma is less than 0, this will throw an error
        
        #get all possible boundary draws for each intensity level
        for (bound in 1:7) {
          boundary_draws[, (bound)] <- d_to_ori(get_draws(boundaries[bound], d_noise), sigmas[1]) 
          boundary_draws[, (bound+7)] <- d_to_ori(get_draws(boundaries[bound], d_noise), sigmas[2])
          boundary_draws[, (bound+14)] <- d_to_ori(get_draws(boundaries[bound], d_noise), sigmas[3])
          boundary_draws[, (bound+21)] <- d_to_ori(get_draws(boundaries[bound], d_noise), sigmas[4])
        } 
        
        
        #check that boundary draws meet conditions
        if (all(boundary_draws > boundaries_bounds[1]) & #lower bound check
            all(boundary_draws < boundaries_bounds[2]) & #upper bound check - probably not necessary
            all((diff(boundary_draws[51,])[-seq(7, 22, 7)]) > 0)) { #check that boundaries are monotonically increasing 
          break   #stop sampling if conditions have been met  
        }
      }
      
      #if sampling is taking too long (which it shouldn't), stop the sampling
      run_time = Sys.time()
      if (difftime(run_time, start_time, units = "mins") > 10) {
        time_out = TRUE
        break
      }
      
    }
    
  } else {
    
    low_noise <-  matrix(0, 101, 28)
    high_noise <-  matrix(0, 101, 28)
    repeat {
      # sampling sigma first
      sigmas <- (rnorm(length(sigma_mu), 0, 1)*sigma_sd)+sigma_mu
      bounds <- (rnorm(length(boundaries_mu), 0, 1)*boundaries_sd)+boundaries_mu
      boundaries = cumsum(bounds[1:7])
      
      psi <- rnorm(1, psi_mu, psi_sd)
      d_noise <- rnorm(1, d_noise_mu, d_noise_sd)
      
      #check if sampled sigma parameters meet conditions
      if ((sum(diff(sigmas) < 0) == (length(sigmas)-1)) & #order check
          all(sigmas > sigma_bounds[1]) & #lower bound check
          all(sigmas < sigma_bounds[2]) & #upper bound check 
          all(diff(round(sigmas, 1)) < 0) & #make sure that sigmas are at least 1 decimal point different from one another
          d_noise > d_noise_bounds[1] &
          d_noise < d_noise_bounds[2] &
          psi > psi_bounds[1] &
          psi < psi_bounds[2]) { 
        
        #I ACTUALLY HAVE NO IDEA IF IT MAKES SENSE TO DO IT LIKE THIS 
        #get all possible boundary draws for each intensity level
        #I am testing the smallest resonable orientation here (0.1) and largest (40)
        for (bound in 1:7) {
          low_noise[, (bound)] <- d_to_ori(get_draws(boundaries[bound], d_noise), sigmas[1] + odn_function(psi, 1)) 
          low_noise[, (bound+7)] <- d_to_ori(get_draws(boundaries[bound], d_noise), sigmas[2] + odn_function(psi, 1))
          low_noise[, (bound+14)] <- d_to_ori(get_draws(boundaries[bound], d_noise), sigmas[3] + odn_function(psi, 1))
          low_noise[, (bound+21)] <- d_to_ori(get_draws(boundaries[bound], d_noise), sigmas[4] + odn_function(psi, 1))
          high_noise[, (bound)] <- d_to_ori(get_draws(boundaries[bound], d_noise), sigmas[1] + odn_function(psi, 46)) 
          high_noise[, (bound+7)] <- d_to_ori(get_draws(boundaries[bound], d_noise), sigmas[2] + odn_function(psi, 46))
          high_noise[, (bound+14)] <- d_to_ori(get_draws(boundaries[bound], d_noise), sigmas[3] + odn_function(psi, 46))
          high_noise[, (bound+21)] <- d_to_ori(get_draws(boundaries[bound], d_noise), sigmas[4] + odn_function(psi, 46))
        } 
        
        #check that boundary draws meet conditions
        if (all(low_noise > boundaries_bounds[1]) & #lower bound check
            all(low_noise < boundaries_bounds[2]) & #upper bound check - probably not necessary
            all(high_noise > boundaries_bounds[1]) & #lower bound check
            all(high_noise < boundaries_bounds[2]) & 
            all((diff(low_noise[51,])[-seq(7, 22, 7)]) > 0)) { #check that boundaries are monotonically increasing 
          break   #stop sampling if conditions have been met  
        }
      }
      
      #if sampling is taking too long (which it shouldn't), stop the sampling
      run_time = Sys.time()
      if (difftime(run_time, start_time, units = "mins") > 10) {
        time_out = TRUE
        break
      }
      
    }
}




#just so this function doesn't throw an error, if sampling has stopped
if (time_out == TRUE) {
  generating_parameters <- rep(NA, length(mus))
} else {
  #return parameters
  if (length(mus) == 12) {
    generating_parameters <- c(bounds, sigmas, d_noise)
  } else {
    generating_parameters <- c(bounds, sigmas, psi, d_noise) 
  }
}

return(generating_parameters)
}
