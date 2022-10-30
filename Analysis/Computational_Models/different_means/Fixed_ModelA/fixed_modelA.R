fixed_model_function <- function(sub, reps = 50, starting_parameters, modality = 'visual') {
  #============= Converting to Z Scores ======================
  #setwd("/home/s4323621/confidence")
  setwd("/Users/s4323621/Dropbox/Documents/Cluster_Files")
  
  if (modality == 'visual') {
    confidencedata <- read.csv(file = 'visual_taskA_data.csv')
  } 
  if (modality == 'auditory') {
    confidencedata <- read.csv(file = 'visual_taskA_auditory_data.csv')
  }
  confidencedata$stim_orientation <- scale(confidencedata$stim_orientation)
  
  #filter for that subject's data
  name <- unique(confidencedata$subject_name)[sub]
  confidencedata <- confidencedata[confidencedata$subject_name == name,] 
  #=================== Set Inital Values =================
  #starting_parameters = c(0.50, 0.36, 0.20, 0.06, 0.6, 0.5, 0.4, 0.3)
  parameter_names <- c("b1", "b2", "b3","b4", 
                       "sigma1", "sigma2", "sigma3","sigma4")
  
  #=================== Flexible Coding  =================
  n_params = length(starting_parameters)
  n_intensities = length(unique(confidencedata$stim_reliability_level))
  
  #=================== Create Storage =================
  
  #create dataframe to store parameter estimates
  parameters <- data.frame(matrix(0, reps, (n_params)))
  colnames(parameters) <- parameter_names
  
  #for fitting group model
  if (length(sub) > 1) {
    parameters$subject <- 0
  } else {
    parameters$subject <- sub  
  }
  
  #create dataframe to store likelihoods in each rep
  likelihoods <- matrix(0, nrow(confidencedata), 1)
  
  #need to set an initial value for this to work
  log.likelihoods <- 0
  
  for (rep in 1:reps) {   
    
    #likelihood function 
    fixed_model <- function(par){
      
      sigma = par[(4+confidencedata$stim_reliability_level)] #get baseline sigma for each reliability level
      
      b_params = c(-Inf, -par[1:3], par[4], rev(par[1:3]), Inf)
      
      if (any(par[5:8] <= 0) | #none of the sigma parameters can be less than 0 
           any(sort(b_params) != b_params)| #make sure boundaries are monotonically increasing
           (sum(diff(par[5:8]) <= 0) != 3)) { #make sure sigma is monotnically increasing
            log.likelihoods <- Inf
          }
          
          if (log.likelihoods < Inf) {
            
            #get likelihood for upper boundary
            likelihoods = (pnorm(b_params[(confidencedata$r + 1)], mean = confidencedata$stim_orientation, 
                                  sd = sigma)) - 
                           #get likelihood for lower boundary 
                           (pnorm(b_params[confidencedata$r], mean = confidencedata$stim_orientation, 
                                  sd = sigma))
            
            #where likelihoods have been rounded to 0, set to minimum value
            likelihoods[likelihoods == 0] <- min(likelihoods[likelihoods>0])
            
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
    parameters$AIC[rep] <- bic
    parameters$BIC[rep] <- aic
    
    rep = rep + 1 
    
    print(rep)
  }
  
  #save parameters from the best fit
  best_estimate = which.min(parameters$log_likelihood)
  parameters <- data.frame(parameters[best_estimate,])
  
  return(parameters)
}
