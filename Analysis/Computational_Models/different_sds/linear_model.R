linear_model_function <- function(confidencedata, sub, reps = 100, starting_parameters, model) {
#Rebecca West 2022
#This function estimates parameters for the scaled evidence strength model class
#Models include linear, quadratic and exponent 
#============= Function Arguments ==============
  #1. confidencedata = data that we are fitting
  #2. sub = subject number
  #3. Number of times to run optimisation 
  #each repetition uses the values returned by optim in the previous iteration as starting values 
  #4. Starting values for first repetition
  #5. Model - either 'lin', 'quad' or 'exp'
#============= Standardise Data (if needed) ======================
#confidencedata$stim_orientation <- scale(confidencedata$stim_orientation)
# if fitting data from auditory modality 
  if (all(confidencedata$stim_orientation == 0)) {
    confidencedata$stim_orientation = confidencedata$stim_frequency
  }
#=================== Set Inital Values =================
n_params = length(starting_parameters)
n_intensities = length(unique(confidencedata$stim_reliability_level))
  
#define function to calculate boundary parameters   
boundary <- function(k, m, s, exp) {
    b = k+(m*(s^exp))
    return(b)
  }
  #=================== Create Storage =================
  #create dataframe to store parameter estimates
  parameters <- data.frame(matrix(0, reps, (n_params)))
  
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
  
  #=================== Start Optimisation Loop =================
  for (rep in 1:reps) {   
    
    #likelihood function 
    linear_model <- function(par) {
      
      k_params = c(0, par[1:7], Inf)
      m_params = c(0, par[8:14], Inf)
      sigs = par[15:18]
      
      if (model == 'lin') {
        exp = 1
      } else if (model == 'quad') {
        exp = 2
      } else if (model == 'exp') {
        exp = par[19]
      }
        bs <- sapply(1:length(k_params), function(x) 
          boundary(k_params[x], m_params[x], sigs, exp))

      if (any(bs[,2:8] <= 0) | #none of the boundaries can be less than 0 
          #any(bs[,2:8] > 15) |
          any(sigs > 20) |
          exp < 1 |
          any(sigs <= 0) | #none of the sigma parameters can be less than 0 
          any(t(apply(bs, 1, function(x) sort(x))) != bs)) { #make sure boundaries are increasing
        log.likelihoods <- Inf
      }
      
      if (log.likelihoods < Inf) {
       
        sigma = sigs[confidencedata$stim_reliability_level] #get baseline sigma for each reliability level
        
        #get likelihood for upper postive boundary
        likelihoods = ((pnorm(boundary(k_params[(confidencedata$r + 1)], m_params[(confidencedata$r + 1)], sigma, exp), mean = confidencedata$stim_orientation, 
                              sd = sigma)) - 
                         #get likelihood for lower positive boundary 
                         (pnorm(boundary(k_params[(confidencedata$r)], m_params[(confidencedata$r)], sigma, exp), mean = confidencedata$stim_orientation, 
                                sd = sigma))) + 
          #get likelihood for upper negative boundary
          ((pnorm(-(boundary(k_params[(confidencedata$r)], m_params[(confidencedata$r)], sigma, exp)), mean = confidencedata$stim_orientation, 
                  sd = sigma)) -
             #get likelihood for lower positive boundary 
             (pnorm(-(boundary(k_params[(confidencedata$r + 1)], m_params[(confidencedata$r + 1)], sigma, exp)), mean = confidencedata$stim_orientation, 
                    sd = sigma)))
        }
        
        likelihoods[likelihoods <= 0] <- min(likelihoods[likelihoods > 0])
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
