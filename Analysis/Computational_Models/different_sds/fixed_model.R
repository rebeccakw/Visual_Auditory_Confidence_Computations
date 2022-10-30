fixed_model_function <- function(confidencedata, sub, reps = 100, starting_parameters) {
  #Rebecca West 2022
  #This function estimates parameters for the fixed model (unscaled evidence strength class)
  #============= Function Arguments ==============
  #1. confidencedata = data that we are fitting
  #2. sub = subject number
  #3. Number of times to run optimisation 
    #each repetition uses the values returned by optim in the previous iteration as starting values 
  #4. Starting values for first repetition
  
  #============= Standardise Data ======================
  #confidencedata$stim_orientation <- scale(confidencedata$stim_orientation)
  #confidencedata$stim_frequency <- scale(confidencedata$stim_frequency)
  #hacky way of dealing with auditory data
  if (all(confidencedata$stim_orientation == 0)) {
    confidencedata$stim_orientation = confidencedata$stim_frequency
  }
  #=================== Set Inital Values =================
  parameter_names <- c("bound1", "bound2", "bound3","bound4", "bound5", "bound6", "bound7", 
                       "sigma1", "sigma2", "sigma3","sigma4")
  n_params = length(parameter_names)
  n_intensities = length(unique(confidencedata$stim_reliability_level))
  
  #=================== Create Storage =================
  #create dataframe to store parameter estimates
  parameters <- data.frame(matrix(0, reps, (n_params)))
  colnames(parameters) <- parameter_names
  
  #for fitting group model
  if (length(sub) > 1) {
    parameters$subject <- 0
  } else {
  #filter for individual subject if fitting individual
    parameters$subject <- sub 
    confidencedata <- confidencedata[confidencedata$subject_name == sub,] 
  }
  
  #create dataframe to store likelihoods in each rep
  likelihoods <- matrix(0, nrow(confidencedata), 1)
  
  #need to set an initial value for this to work
  log.likelihoods <- 0
  
  for (rep in 1:reps) {   
    
    #likelihood function 
    fixed_model <- function(par){
      
      #test for any undesirable parameter values 
      if (any(par[1:11] <= 0) | #none of the boundaries or sigmas can be less than 0 
        (sum(diff(par[1:7]) > 0) != 6) | #boundaries are monotonically increasing
        (sum(diff(par[8:11]) <= 0) != 3)) { #{ #make sure sigma is monotnically increasing
        #any(par[8:11] > 10) | #upper limit on sigma values
        #any(par[1:7] > 5)) { #upper limit on boundary parameters
        log.likelihoods <- Inf
      }
      
      if (log.likelihoods < Inf) {
        
        sigma = par[(7+confidencedata$stim_reliability_level)] #get baseline sigma for each reliability level
        boundaries = c(0, par[1:7], Inf)
        
        #get likelihood for upper postive boundary
        likelihoods = ((pnorm((boundaries[(confidencedata$r + 1)]), mean = confidencedata$stim_orientation, 
                              sd = sigma)) - 
                         #get likelihood for lower positive boundary 
                         (pnorm((boundaries[confidencedata$r]), mean = confidencedata$stim_orientation, 
                                sd = sigma))) + 
          #get likelihood for upper negative boundary
          ((pnorm(-(boundaries[confidencedata$r]), mean = confidencedata$stim_orientation, 
                  sd = sigma)) -
             #get likelihood for lower positive boundary 
             (pnorm(-(boundaries[(confidencedata$r + 1)]), mean = confidencedata$stim_orientation, 
                    sd = sigma)))
        
        likelihoods[likelihoods <= 0] = min(likelihoods[likelihoods > 0])
        
        #calculate 
        log.likelihoods <- -sum(log(likelihoods))
      }
      
      return(log.likelihoods)
    }
    
    
    
    # =================== Optimisation =================
    
    #if first rep - use starting_parameters
    if (rep == 1) { 
      parameter.fits <- optim(par = starting_parameters,
                              fn = fixed_model)
    } 
    
    # if rep is greater than one - use parameter values from previous iteration
    else {
      parameter.fits <- optim(par = as.double(parameters[(rep - 1), 1:n_params]),
                              fn = fixed_model) 
    }
    
    
    #save parameters and log likelihood for fit
    parameters[rep, 1:n_params] <- parameter.fits$par
    parameters$log_likelihood[rep] <- parameter.fits$value
    
    #calculate and save model comparison metrics (AIC and BIC)
    bic <- (-2*(-1*parameter.fits$value)) + length(parameter.fits$par) * log(nrow(confidencedata))
    aic <- (-2*(-1*parameter.fits$value)) + (2*length(parameter.fits$par))
    parameters$AIC[rep] <- aic
    parameters$BIC[rep] <- bic
    
    rep = rep + 1 
    
    print(rep)
  }
  
  #save parameters from the best fit
  best_estimate = which.min(parameters$log_likelihood)
  parameters <- data.frame(parameters[best_estimate,])
  
  return(parameters)
}
