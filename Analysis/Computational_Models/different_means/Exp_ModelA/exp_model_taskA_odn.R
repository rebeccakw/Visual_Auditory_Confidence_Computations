exp_model_function_odn <- function(confidencedata, sub, reps = 50, starting_parameters, modality = 'visual') {
  #============= Converting to Z Scores ======================
  sample_mean = mean(confidencedata$stim_orientation)
  sample_sd = sd(confidencedata$stim_orientation)
  
  zscore <- function(value) {
    z = (value - sample_mean)/sample_sd
    return(z)
  }
  
  confidencedata$unscaled_stim_orientation <- confidencedata$stim_orientation
  confidencedata$stim_orientation <- scale(confidencedata$stim_orientation)

  boundary_function <- function(k, m, sigma, exp) {
    b = k + (m*(sigma^exp))
    return(b)
  }
  #=================== Set Inital Values =================
  parameter_names <- c("k1", "k2", "k3","k4", 
                       "m1", "m2", "m3", "m4",
                       "sigma1", "sigma2", "sigma3","sigma4", "exp", "psi")
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
  
  for (rep in 1:reps) {   
    
    #likelihood function 
    model <- function(par){

      if (any(par[9:14] <= 0) | #none of the sigma parameters can be less than 0 
          par[13] < 1|
          any(par[9:12] > 20) | #upper limit on sigma parameters
          (sum(diff(par[9:12]) <= 0) != 3)) { #make sure sigma is monotnically increasing
        log.likelihoods <- Inf
      }
      if (log.likelihoods < Inf) {
      #rename parameters 
      kbounds = par[1:4]
      mbounds = par[5:8]
      exp = par[13]
      
      #get sigma for all trials 
      baseline_sigma = par[(8 + confidencedata$stim_reliability_level)]
      odn = zscore(par[14]*(abs(sin((confidencedata$unscaled_stim_orientation*pi)/90))))
      all_sigma = baseline_sigma + odn
      
      #calculate boundaries for all trials 
      boundaries = t(sapply(all_sigma, function(x) boundary_function(kbounds, mbounds, x, exp)))
      #get upper and lower boundaries 
      all_boundaries = cbind(-Inf, #lower limit
                             -boundaries[,1:3], #make first 3 boundaries negative
                             boundaries[,4], #category boundary
                             t(apply(boundaries[,1:3], 1, function(x) rev(x))), #positive boundaries
                             Inf) #upper limit 
     
      #make sure boundaries are monotonically increasing 
       if (any(t(apply(all_boundaries, 1, function(x) sort(x))) != all_boundaries)) {
         log.likelihoods <- Inf
       }
      
      #get upper boundary for trial responses
      b_high = all_boundaries[cbind(1:nrow(confidencedata), (confidencedata$r + 1))]
      #get lower boundary for trial responses
      b_low = all_boundaries[cbind(1:nrow(confidencedata), (confidencedata$r))]
      
      if (log.likelihoods < Inf) {
        
        #get likelihood for upper boundary
        likelihoods = (pnorm(b_high, mean = confidencedata$stim_orientation, sd = all_sigma)) - 
          (pnorm(b_low, mean = confidencedata$stim_orientation, sd = all_sigma))
  
        
        #where likelihoods have been rounded to 0, set to minimum value
        likelihoods[likelihoods <= 0] <- min(likelihoods[likelihoods>0])
        
        #calculate 
        log.likelihoods <- -sum(log(likelihoods))
      }
      }
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
