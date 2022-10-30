exp_model_function <- function(sub, reps = 50, starting_parameters) {
  #============= Converting to Z Scores ======================
  #setwd("/home/s4323621/confidence")
  setwd("/Users/rebeccawest/Dropbox/Documents/Cluster_Files")
  visual <- read.csv(file = 'visual_taskA_data.csv')
  auditory <- read.csv(file = 'visual_taskA_auditory_data.csv')
  
  visual$stim_orientation <- scale(visual$stim_orientation)
  auditory$stim_orientation <- scale(auditory$stim_orientation)
  confidencedata <- rbind(visual, auditory)
  
  boundary_function <- function(k, m, sigma, exp) {
    b = k + (m*(sigma^exp))
    return(b)
  }
  #=================== Set Inital Values =================
  parameter_names <- c("k1", "k2", "k3","k4", 
                       "m1", "m2", "m3", "m4",
                       "sigma1", "sigma2", "sigma3","sigma4", "exp", "scaler")
  #=================== Flexible Coding  =================
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
    parameters$subject <- sub  
    #filter for that subject's data
    confidencedata <- confidencedata[confidencedata$subject_name == sub,] 
  }
  
  #create dataframe to store likelihoods in each rep
  likelihoods <- matrix(0, nrow(confidencedata), 1)
  
  #need to set an initial value for this to work
  log.likelihoods <- 0
  boundaries_all = matrix(0, 9, 8)
  
  for (rep in 1:reps) {   
    
    #likelihood function 
    model <- function(par){
      
      if (any(par[9:14] <= 0) | #none of the sigma parameters can be less than 0 
          par[13] < 1|
          #any(par[5:8] > 20) | #upper limit on sigma parameters
          (sum(diff(par[9:12]) <= 0) != 3) |
          par[9:12] > 20) { #make sure sigma is monotnically increasing
        log.likelihoods <- Inf
      }
      
      if (log.likelihoods < Inf) {
      
      sigma = par[(8 + confidencedata$stim_reliability_level)]
      sigma[confidencedata$stim_type == 'Auditory'] = par[14]*sigma[confidencedata$stim_type == 'Auditory']
      exp = par[13]
      
      k_params = c(Inf, rev(par[2:4]), par[1:4], Inf)
      m_params = c(Inf, rev(par[6:8]), par[5:8], Inf)
      
      sigmas_all = c(par[9:12], par[14]*par[9:12])     
      
      for (s in 1:length(sigmas_all)) {
        boundaries_all[,s] <- boundary_function(k_params, m_params, sigmas_all[s], exp)
      }
      boundaries_all[1:4,] = -boundaries_all[1:4,]
      
      if (any(apply(boundaries_all, 2, sort) != (boundaries_all))) {
        log.likelihoods <- Inf
      }
    }
      
      if (log.likelihoods < Inf) {
        b_low = boundary_function(k_params[confidencedata$r], m_params[confidencedata$r], sigma, exp)
        b_high = boundary_function(k_params[confidencedata$r + 1], m_params[confidencedata$r + 1], sigma, exp)
        
        b_low[confidencedata$r <= 4] = -b_low[confidencedata$r <= 4]
        b_high[confidencedata$r <= 3] = -b_high[confidencedata$r <= 3]
        
        #get likelihood for upper boundary
        likelihoods = (pnorm(b_high, mean = confidencedata$stim_orientation, sd = sigma)) - 
          (pnorm(b_low, mean = confidencedata$stim_orientation, sd = sigma))
      }
      
      #where likelihoods have been rounded to 0, set to minimum value
      likelihoods[likelihoods == 0] <- min(likelihoods[likelihoods>0])
      
      #calculate 
      log.likelihoods <- -sum(log(likelihoods))
      
      return(log.likelihoods)
    }
    
    
    
    # =================== Optimisation =================
    
    #if first rep - use starting_parameters
    if (rep == 1) { 
      parameter.fits <- optim(par = starting_parameters,
                              fn = model)
    } 
    
    # if rep is greater than one - use parameter values from previous iteration
    else {
      parameter.fits <- optim(par = as.double(parameters[(rep - 1), 1:n_params]),
                              fn = model) 
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
