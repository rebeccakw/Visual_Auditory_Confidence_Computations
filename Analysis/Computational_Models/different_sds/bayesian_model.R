bayesian_model_function <- function(confidencedata, sub, reps = 100, starting_parameters, 
                                    modality = 'visual', fit_prior = FALSE) {
  #Rebecca West 2022
  #This function estimates parameters for the bayesian model class
  #only used for standard bayesian model and bayesian model with free category distribution parameters
  #============= Function Arguments ==============
  #1. confidencedata = data that we are fitting
  #2. sub = subject number
  #3. Number of times to run optimisation 
  #each repetition uses the values returned by optim in the previous iteration as starting values 
  #4. Starting values for first repetition
  #5. Modality 
  #either 'visual' or 'auditory' depending on data 
  #6. fit_prior
  #For the Bayesian model with free category distribution parameters 
  #Gets suitable parameters for the category distributions 
  #==================================================
  ##Standardise Data (if required) 
  #confidencedata$stim_orientation = scale(confidencedata$stim_orientation)
  if (all(confidencedata$stim_orientation == 0)) {
    confidencedata$stim_orientation = confidencedata$stim_frequency #hacky way of dealing with auditory data 
  }
  #=================== Required Functions =================
  #define function to convert from d space to perceptual space
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
  #====================================  
  n_params = length(starting_parameters)
  n_intensities = length(unique(confidencedata$stim_reliability_level))
  
  #create dataframe to store parameter estimates
  parameters <- data.frame(matrix(0, reps, n_params))
  parameter_names <- c("bound1", "bound2", "bound3","bound4", "bound5", "bound6", "bound7", 
                       "sigma1", "sigma2", "sigma3","sigma4")
  if (fit_prior == TRUE) {
    parameter_names = c(parameter_names, "cat_one_sigma", "cat_two_sigma")
  } 
  colnames(parameters) <- parameter_names
  
  if (length(sub) > 1) {
    parameters$subject <- 0
  } else {
    parameters$subject <- sub  
    confidencedata = confidencedata[confidencedata$subject_name == sub,]
  }
  
  #create dataframe to store likelihoods in each rep
  n_trials = nrow(confidencedata)
  b_high <- matrix(Inf, nrow(confidencedata), 1)
  b_low <- matrix(0, nrow(confidencedata), 1)
  
  #need to set an initial value for this to work
  log.likelihoods <- 0
  
  for (rep in 1:reps) {   
    
    bayesian_model <- function(par){
      
      #get relevant category distribution parameters 
      if (fit_prior == TRUE & unimodal == TRUE) {
        cat_one_sigma = par[12]
        cat_two_sigma = par[13]
      } else if (fit_prior == FALSE) {
        if (modality == 'visual') {
          cat_one_sigma = 0.3427148
          cat_two_sigma = 1.370859
        } else if (modality == 'auditory') {
          cat_one_sigma = 0.3429777
          cat_two_sigma = 1.371911
        }
      }
      
      # test for any undesirable parameter values 
      if (any(par[8:length(par)] <= 0) | #none of the sigma parameters can be 0 
          #any(par[8:11] > 3) |
          #any(cumsum(par[1:7]) < -10) |
          any(diff(round(par[8:11], 1)) >= 0)) { #make sure sigma parameters are at least one decimal point different
        log.likelihoods <- Inf
      }
      
      if (log.likelihoods < Inf) { 
        sigma = par[(7+confidencedata$stim_reliability_level)] #get sigma for each intensity level
          
        #calculate boundary parameters for each value of sigma   
        bs <- t(sapply(par[8:11], function(x) d_to_ori(cumsum(par[1:7]), x, cat_one_sigma, cat_two_sigma)))
          
          if (any(bs < 0) | #make sure lower boundary positions are not less than 0 
              any(t(apply(bs, 1, function(x) sort(x))) != bs)) { #order check
            log.likelihoods <- Inf
          }
          
          if (log.likelihoods < Inf) { 
          
          boundaries = c(0, cumsum(par[1:7]), Inf)
          lower_boundaries <- boundaries[confidencedata$r]
          upper_boundaries <- boundaries[(confidencedata$r + 1)]
          
          #transform boundary parameters to orientation space using sigma for that trial
          b_high[confidencedata$r != 8,] = d_to_ori(upper_boundaries[confidencedata$r != 8], sigma[confidencedata$r != 8],
                                                    cat_one_sigma, cat_two_sigma)
          b_low[confidencedata$r != 1,] = d_to_ori(lower_boundaries[confidencedata$r != 1], sigma[confidencedata$r != 1],
                                                   cat_one_sigma, cat_two_sigma)
          
          orientations <- confidencedata$stim_orientation
          
          #get likelihood for upper positive boundary 
          likelihoods = (((pnorm(b_high, mean = orientations, 
                                         sd = sigma)) - 
                                    #get likelihood for lower positive boundary 
                                    (pnorm(b_low, mean = orientations, 
                                           sd = sigma))) + 
                                   #get likelihood for upper negative boundary
                                   ((pnorm(-b_low, mean = orientations, 
                                           sd = sigma)) -
                                      #get likelihood for lower positive boundary 
                                      (pnorm(-b_high, mean = orientations, 
                                             sd = sigma))))
          #am I allowed to do this?
          likelihoods[likelihoods == 0] <- min(likelihoods[likelihoods>0])
          
          log.likelihoods = -sum(log(likelihoods))
          }
      }
      
      return(log.likelihoods)
  }
    
    # =================== Optimisation =================
    
    #if first rep - use starting_parameters
    if (rep == 1) { 
      parameter.fits <- optim(par = starting_parameters,
                              fn = bayesian_model)
    } 
    
    # if rep is greater than one - use parameter values from previous iteration
    else {
      parameter.fits <- optim(par = as.double(parameters[(rep - 1), 1:n_params]),
                              fn = bayesian_model) 
    }
    
    
    parameters[rep, 1:n_params] <- parameter.fits$par
    parameters$log_likelihood[rep] <- parameter.fits$value
    
    #calculate and save BIC & AIC
    bic <- (-2*(-1*parameter.fits$value)) + length(parameter.fits$par) * log(nrow(confidencedata))
    aic <- (-2*(-1*parameter.fits$value)) + (2*length(parameter.fits$par))
    parameters$AIC[rep] <- aic
    parameters$BIC[rep] <- bic
    
    print(rep) 
    
    rep = rep + 1 
    
  }
  
  best_estimate = which.min(parameters$log_likelihood) #get best fitting parameters 
  estimated_parameters <- data.frame(parameters[best_estimate,])
  
  return(estimated_parameters)
  
}


