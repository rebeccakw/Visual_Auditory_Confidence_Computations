bayesian_model_function <- function(sub, reps = 50, starting_parameters, modality = 'visual') {
  #============= Converting to Z Scores ======================
  #setwd("/home/s4323621/confidence")
  setwd("/Users/rebeccawest/Dropbox/Documents/Cluster_Files")
  
  if (modality == 'visual') {
    confidencedata <- read.csv(file = 'visual_taskA_data.csv')
  } 
  if (modality == 'auditory') {
    confidencedata <- read.csv(file = 'visual_taskA_auditory_data.csv')
  }
  
  #============= Bayesian Specific Stuff=========================
  if (modality == 'visual') {
    cat_mean = -4
    cat_sigma = 5
  } 
  if (modality == 'auditory') {
    cat_mean = 2300
    cat_sigma = 475
  }
  
  
  #write function to convert from d space to measurement space for each boundary 
  measurement <- function(d, sigma)  {
    measurement = (d*((sigma^2)+(cat_sigma^2)))/(2*cat_mean)
    return(measurement)  
  }

  #confidencedata$stim_orientation <- scale(confidencedata$stim_orientation)

  #=================== Set Inital Values =================
  parameter_names <- c("b1", "b2", "b3","b4", 
                       "sigma1", "sigma2", "sigma3","sigma4", "psi")
  
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
    #filter for that subject's data
    confidencedata <- confidencedata[confidencedata$subject_name == sub,] 
  }
  
  #create dataframe to store likelihoods in each rep
  likelihoods <- matrix(0, nrow(confidencedata), 1)
  
  #need to set an initial value for this to work
  log.likelihoods <- 0
  
  for (rep in 1:reps) {   
    
    #likelihood function 
    model <- function(par){
      
      #boundaries are negative here because measurement() function turns them positive
      #and I change the sign manually later
      boundaries = c(-Inf, par[1:4], rev(par[1:3]), -Inf) 
      
      if (any(par[5:9] <= 0) | #none of the sigma parameters can be less than 0 
          any(par[5:8] > 20) | #upper limit on sigma parameters
          any(sort(boundaries[2:5]) != boundaries[2:5]) | #make sure boundaries are monotonically increasing
          (sum(diff(par[5:8]) <= 0) != 3)) { #make sure sigma is monotnically increasing
        log.likelihoods <- Inf
      }
        
      #confidencedata_rel = confidencedata[confidencedata$stim_reliability_level == rel,]
      baseline_sigma = par[(4 + confidencedata$stim_reliability_level)]
      odn = par[9]*(abs(sin((confidencedata$stim_orientation*pi)/90)))
      sigma = baseline_sigma + odn
      
      lower_boundaries = boundaries[confidencedata$r]
      upper_boundaries = boundaries[(confidencedata$r + 1)]
      b_low = apply(cbind(lower_boundaries, sigma), 1, function(x) measurement(x[1], x[2]))
      b_high = apply(cbind(upper_boundaries, sigma), 1, function(x) measurement(x[1], x[2]))
      b_low[confidencedata$r <= 4] = -b_low[confidencedata$r <= 4]
      b_high[confidencedata$r <= 3] = -b_high[confidencedata$r <= 3]
      
      
      if (log.likelihoods < Inf) {
        
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
