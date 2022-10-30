generate_bayesian_parameters <- function(mus, sigmas, confidencedata, modality = 'visual', fit_prior = FALSE) {
  #============================
  #Rebecca West 2022
  #This function randomly generates parameters for the Bayesian model class
  #Parameters include: ks + ms (boundaries), sigmas
  #Parameters are sampled from a normal distribution with a defined mean (referred to with mu) 
  #and standard deviation (referred to as sd). All parameters also have a lower and upper bound. 
  
  #Function continues to resample parameters from the define distribution until all parameters meet 
  #the requirements of that model
  
  #Note: unlike the other sampling functions for this project, I pass means and sds of the 
  #distributions to sample parameter values from as arguments 
  #================== Function Arguments ===================
  #1. mus  
      #1-7: boundaries
      #8-11: sigma4 (largest), sigma3, sigma2, sigma1 (smallest)
      #12: psi (if it exists - orientation dependent noise parameter)
      #13: dnoise (decision noise parameters)
  #2. sigmas 
      #as above but for standard deviations 
  #3. confidencedata = data (needed for psi parameter)
  #4. modality 
    #either 'visual' or 'auditory' depending on data 
  #5. fit_prior
    #For the Bayesian model with free category distribution parameters 
    #Gets suitable parameters for the category distributions 
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
  sigma_bounds = c(0, 3.5)
  
  #true category distribution SDs
  if (modality == 'visual') {
    cat_one_sigma_mu = 0.3427148
    cat_two_sigma_mu = 1.370859
  } else if (modality == 'auditory') {
    cat_one_sigma_mu = 0.3429777
    cat_two_sigma_mu = 1.371911
  }
  cat_sd = sigmas[12]
  cat_bounds = c(0, 3)
  
  #define function to convert from d space to perceptual space for each boundary 
  d_to_ori <- function(d, s, cat_one_sigma, cat_two_sigma) {
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
  
  
  #set a timer because sometimes this will run on an infinite loop
  start_time = Sys.time()
  time_out = FALSE
  
  #================================
  # Sample Boundary & Sigma Parameters
  #================================
    boundary_draws <-  matrix(0, 7, 4)
    repeat {
      # sampling
      sigmas <- (rnorm(length(sigma_mu), 0, 1)*sigma_sd)+sigma_mu
      bounds <- (rnorm(length(boundaries_mu), 0, 1)*boundaries_sd)+boundaries_mu
      boundaries = cumsum(bounds[1:7])
      
      ## find good category parameters
      if (fit_prior == TRUE & (modality == 'visual' | modality == 'auditory')) {
      repeat { 
      cat_one_sigma <- (rnorm(length(cat_one_sigma_mu), 0, 1)*cat_sd)+cat_one_sigma_mu
      cat_two_sigma <- (rnorm(length(cat_two_sigma_mu), 0, 1)*cat_sd)+cat_two_sigma_mu
      values = c(cat_one_sigma, cat_two_sigma)
      if (all(values > cat_bounds[1]) & all(values < cat_bounds[2])) {
        break
      }
      }
        #if we aren't fitting category distribution, use true values 
      } else if (fit_prior == FALSE & modality == 'visual' | modality == 'auditory') {
        cat_one_sigma <- cat_one_sigma_mu
        cat_two_sigma <- cat_two_sigma_mu
      } 
      
      #check if sampled sigma parameters meet conditions
      if ((sum(diff(sigmas) < 0) == (length(sigmas)-1)) & #order check
          all(sigmas > sigma_bounds[1]) & #lower bound check
          all(sigmas < sigma_bounds[2]) & #upper bound check 
          all(diff(round(sigmas, 1)) < 0)) { #make sure that sigmas are at least 1 decimal point different from one another
        
        #if sampled sigma parameters meet condition, check boundary parameters
        #have to check sigmas first because if sigma is less than 0, this will throw an error
        
        #get all possible boundary draws for each intensity level
        for (rel in 1:4)  {
        boundary_draws[,rel] <- d_to_ori(boundaries, sigmas[rel], cat_one_sigma, cat_two_sigma) 
        }
        
        #check that boundary draws meet conditions
        if (all(boundary_draws > boundaries_bounds[1]) & #lower bound check
            all(boundary_draws < boundaries_bounds[2]) & #upper bound check - probably not necessary
            all(t(apply(boundary_draws, 1, function(x) sort(x))) == boundary_draws)) { #check that boundaries are monotonically increasing 
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

#just so this function doesn't throw an error, if sampling has stopped
if (time_out == TRUE) {
  if (fit_prior == TRUE & (modality == 'visual' | modality == 'auditory')) {
    generating_parameters <- rep(NA, 13)
 } else {
    generating_parameters <- rep(NA, 11)
  }
} else {
  if (fit_prior == TRUE & (modality == 'visual' | modality == 'auditory')) {
    generating_parameters <- c(bounds, sigmas, cat_one_sigma, cat_two_sigma)
  } else {
    generating_parameters <- c(bounds, sigmas)
  }
}
return(generating_parameters)
}
