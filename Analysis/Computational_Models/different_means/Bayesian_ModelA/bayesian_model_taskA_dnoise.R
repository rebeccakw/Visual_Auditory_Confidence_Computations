bayesian_model_function_dnoise <- function(sub, reps = 50, starting_parameters, modality = 'visual') {
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
    cat_mean = -0.6247561
    cat_sigma = 0.7809452
  } 
  if (modality == 'auditory') {
    cat_mean = -0.6449836
    cat_sigma = 0.7608733
  }
  
  
  #write function to convert from d space to measurement space for each boundary 
  measurement <- function(d, sigma)  {
    measurement = (d*((sigma^2)+(cat_sigma^2)))/(2*cat_mean)
    return(measurement)  
  }
  
  get_draws <- function(d_bound, d_sd) {
   if (abs(d_bound) == Inf) { #don't get draws for upper and lower limits
      spaced_draws = rep(d_bound, 101)
    } else {
      spaced_draws <- seq((qnorm(0.01, mean = d_bound, sd = d_sd)), (qnorm(0.99, mean = d_bound, sd = d_sd)), length.out = 101)
    }
    return(spaced_draws)
  }

  confidencedata$stim_orientation <- scale(confidencedata$stim_orientation)

  #=================== Set Inital Values =================
  parameter_names <- c("b1", "b2", "b3","b4", 
                       "sigma1", "sigma2", "sigma3","sigma4", "d_noise")
  
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
      
      if (any(par[5:9] <= 0) | #none of the sigma/dnoise parameters can be less than 0 
          any(par[5:8] > 20) | #upper limit on sigma parameters
          any(sort(cumsum(par[1:4]), decreasing = TRUE) != cumsum(par[1:4])) | #make sure boundaries are monotonically increasing
          (sum(diff(par[5:8]) <= 0) != 3)) { #make sure sigma is monotnically increasing
        log.likelihoods <- Inf
      }
      
      if (log.likelihoods < Inf) {
        
      b = cumsum(par[1:4])
      boundaries = c(-Inf, rev(b[2:4]), b, -Inf)  
      
      #confidencedata_rel = confidencedata[confidencedata$stim_reliability_level == rel,]
      sigma = par[(4 + confidencedata$stim_reliability_level)]
      dnoise = par[9]
      
      #get upper and lower boundary draws for each trial, based on response
      lower_boundaries <- sapply(boundaries[confidencedata$r], function(x) get_draws(x, dnoise))
      upper_boundaries <- sapply(boundaries[(confidencedata$r + 1)], function(x) get_draws(x, dnoise))
      
      b_low = apply(rbind(lower_boundaries, sigma), 2, function(x) measurement(x[1:101], x[102]))
      b_high = apply(rbind(upper_boundaries, sigma), 2, function(x) measurement(x[1:101], x[102]))
      
      #make relevant boundaries negative 
      b_low[,confidencedata$r <= 4] = apply(b_low[,confidencedata$r <= 4], 2, function(x) rev(-x))
      b_high[,confidencedata$r <= 3] = apply(b_high[,confidencedata$r <= 3], 2, function(x) rev(-x))
      
      density = matrix(dnorm(get_draws(boundaries[2], dnoise), boundaries[2], dnoise)/
                         sum(dnorm(get_draws(boundaries[2], dnoise), boundaries[2], dnoise)), 101, ncol(b_low))
      
      #need to transform these for assess likelihoods
      sigmas = matrix(rep(sigma, each = 101), nrow = 101, ncol = length(sigma))
      oris = matrix(rep(confidencedata$stim_orientation, each = 101), nrow = 101, ncol = length(sigma))
        
        #matrix multiplication 
        likelihoods = density*((pnorm(b_high, mean = oris, sd = sigmas)) - 
          (pnorm(b_low, mean = oris, sd = sigmas)))
        
        #get a likelihood for each trial by summing across columns 
        likelihoods = apply(likelihoods, 2, function(x) sum(x)) 
        
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
                              fn = model)
    } 
    
    # if rep is greater than one - use parameter values from previous iteration
    else {
      parameter.fits <- optim(par = as.numeric(parameters[(rep - 1), 1:n_params]),
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
    
    print(rep)
    rep = rep + 1 
  }
  
  #save parameters from the best fit
  best_estimate = which.min(parameters$log_likelihood)
  parameters <- data.frame(parameters[best_estimate,])
  
  return(parameters)
}
