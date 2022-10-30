## FUNCTIONS FOR BAYESIAN MODEL 
exp_model_function <- function(sub, reps = 50, starting_parameters) {
  #============= Converting to Z Scores ======================
  #setwd("/home/s4323621/confidence")
  setwd("/Users/rebeccawest/Dropbox/Documents/Cluster_Files")
  visual <- read.csv(file = 'visual_taskA_data.csv')
  auditory <- read.csv(file = 'visual_taskA_auditory_data.csv')
  
  visual$stim_orientation <- scale(visual$stim_orientation)
  auditory$stim_orientation <- scale(auditory$stim_orientation)
  confidencedata <- rbind(visual, auditory)
  #=================== Set Inital Values =================
  parameter_names <- c("k1", "k2", "k3","k4", 
                       "m1", "m2", "m3", "m4",
                       "sigma1", "sigma2", "sigma3","sigma4", "exp")
  
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

    confidencedata <- confidencedata[confidencedata$subject_name == sub,] 
  }
  
  #create dataframe to store likelihoods in each rep
  likelihoods <- matrix(0, nrow(confidencedata), 1)
  
  #need to set an initial value for this to work
  log.likelihoods <- 0
  
  for (rep in 1:reps) {   
    
    #likelihood function 
    linear_model <- function(par){
    
      for (rel in 1:4) {
        confidencedata_rel = confidencedata[confidencedata$stim_reliability_level == rel,]
        sigma = par[(8 + rel)]
        b_params = c(-Inf, -(par[1:3] + (par[5:7]*(sigma^par[13]))), (par[4] + (par[8]*(sigma^par[13]))), rev(par[1:3] + (par[5:7]*(sigma^par[13]))), Inf)
        
        if (any(par[9:13] <= 0) | #none of the sigma parameters can be less than 0 
            any(par[9:13] > 20) |
            par[13] < 1 |
            any(sort(b_params) != b_params) | #make sure boundaries are monotonically increasing
            (sum(diff(par[9:12]) <= 0) != 3)) { #make sure sigma is monotnically increasing
          log.likelihoods <- Inf
        }
        
        if (log.likelihoods < Inf) {
          
          #get likelihood for upper boundary
          likelihoods[((nrow(confidencedata_rel)*(rel-1))+1):(nrow(confidencedata_rel)*rel)] = ((pnorm((b_params[(confidencedata_rel$r + 1)]), mean = confidencedata_rel$stim_orientation, 
                                sd = sigma)) - 
                           #get likelihood for lower boundary 
                           (pnorm((b_params[confidencedata_rel$r]), mean = confidencedata_rel$stim_orientation, 
                                  sd = sigma))) 
          
        
        
        }
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
                              fn = linear_model)
    } 
    
    # if rep is greater than one - use parameter values from previous iteration
    else {
      parameter.fits <- optim(par = as.double(parameters[(rep - 1), 1:n_params]),
                              fn = linear_model) 
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
