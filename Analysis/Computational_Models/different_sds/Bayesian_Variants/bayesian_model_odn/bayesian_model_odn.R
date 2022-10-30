bayesian_model_odn <- function(confidencedata, sub, reps = 100, starting_parameters, modality = 'visual') {
  
  #=================== Required Functions =================
  if (modality == 'visual') {
    cat_one_sigma = 0.3427148
    cat_two_sigma = 1.370859
  } else if (modality == 'auditory') {
    cat_one_sigma = 0.3429777
    cat_two_sigma = 1.371911
  }
  
  confidencedata$unscaled_stim_orientation <- confidencedata$stim_orientation
  confidencedata$stim_orientation = scale(confidencedata$stim_orientation)
  sample_mean = mean(confidencedata$stim_orientation)
  sample_sd = sd(confidencedata$stim_orientation)
  
  #use this for odn noise 
  z_score <- function(value) {
    z = (value - sample_mean)/sample_sd
    return(z)
  }
  
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
  
  #====================================  
  n_params = length(starting_parameters)
  n_intensities = length(unique(confidencedata$stim_reliability_level))
  
  #create dataframe to store parameter estimates
  parameters <- data.frame(matrix(0, reps, n_params))
  parameter_names <- c("bound1", "bound2", "bound3","bound4", "bound5", "bound6", "bound7", 
                       "sigma1", "sigma2", "sigma3","sigma4", "psi")
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
      
      # || is true if at least one of the conditions is True
      if (any(par[8:12] <= 0) | #none of the sigma parameters can be 0 
          #any(par[8:11] > 30) |
          any(diff(round(par[8:11], 1)) >= 0)) { #make sure sigma parameters are at least one decimal point different
        log.likelihoods <- Inf
      }
      
      if (log.likelihoods < Inf) {
      
        baseline_sigma = par[(7+confidencedata$stim_reliability_level)] #get baseline sigma for each reliability level
        odn = z_score(par[12]*(abs(sin((confidencedata$unscaled_stim_orientation*pi)/90))))
        sigma = baseline_sigma + odn
      
        boundaries = t(sapply(sigma, function(x) c(0, d_to_ori(cumsum(par[1:7]), x), Inf)))
       
        if (any(boundaries[,2:8] <= 0)|
            any(t(apply(boundaries, 1, function(x) sort(x))) != boundaries)) {
         log.likelihoods <- Inf 
       }
      
      if (log.likelihoods < Inf) { 
        #get upper and lower boundary draws for each trial, based on response
        lower_boundaries <- boundaries[cbind(1:nrow(confidencedata),confidencedata$r)]
        upper_boundaries <- boundaries[cbind(1:nrow(confidencedata),(confidencedata$r + 1))]
          
          #get density of boundary draws -- using boundary 1 
          #density are normalised (divided by total density) so use density values from boundary 1 for all boundary parameters.
          #different boundaries are just shifts in the mean and d_noise sd stays the same
          orientations <- confidencedata$stim_orientation
          
          #get likelihood for upper positive boundary 
          likelihoods = (((pnorm(upper_boundaries, mean = orientations, 
                                         sd = sigma)) - 
                                    #get likelihood for lower positive boundary 
                                    (pnorm(lower_boundaries, mean = orientations, 
                                           sd = sigma))) + 
                                   #get likelihood for upper negative boundary
                                   ((pnorm(-lower_boundaries, mean = orientations, 
                                           sd = sigma)) -
                                      #get likelihood for lower positive boundary 
                                      (pnorm(-upper_boundaries, mean = orientations, 
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
    parameters$AIC[rep] <- bic
    parameters$BIC[rep] <- aic
    
    print(rep) 
    
    rep = rep + 1 
    
  }
  
  best_estimate = which.min(parameters$log_likelihood) #get best fitting parameters 
  estimated_parameters <- data.frame(parameters[best_estimate,])
  
  #save all parameters as a back-up
  #with data and time
  #write.csv(parameters, paste0(wd, "/all_bayesian_model_parameters_", format(Sys.time(), "%d%B%Y-%X"), ".csv"))
  
  return(estimated_parameters)
  
}


