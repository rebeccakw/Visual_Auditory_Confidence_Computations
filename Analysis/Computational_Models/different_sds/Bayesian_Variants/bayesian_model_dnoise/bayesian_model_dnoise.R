bayesian_model <- function(confidencedata, sub, reps = 100, starting_parameters, modality = 'visual') {
  
  #confidencedata$stim_orientation = scale(confidencedata$stim_orientation)
  
  #=================== Required Functions =================
  if (modality == 'visual') {
    cat_one_sigma = 0.3427148
    cat_two_sigma = 1.370859
  } else if (modality == 'auditory') {
    cat_one_sigma = 0.3429777
    cat_two_sigma = 1.371911
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
  
  get_draws <- function(d_bound, d_sd) {
    if (any(d_bound == Inf | d_bound == 0)) { #don't get draws for upper and lower limits
      spaced_draws = rep(d_bound, 101)
    } else {
      spaced_draws <- seq((qnorm(0.01, mean = d_bound, sd = d_sd)), (qnorm(0.99, mean = d_bound, sd = d_sd)), length.out = 101)
    }
    return(spaced_draws)
  }
  
  #====================================  
  n_params = length(starting_parameters)
  n_intensities = length(unique(confidencedata$stim_reliability_level))
  
  #create dataframe to store parameter estimates
  parameters <- data.frame(matrix(0, reps, n_params))
  parameter_names <- c("bound1", "bound2", "bound3","bound4", "bound5", "bound6", "bound7", 
                       "sigma1", "sigma2", "sigma3","sigma4", "d_noise")
  colnames(parameters) <- parameter_names
  
  if (length(sub) > 1) {
    parameters$subject <- 0
  } else {
    parameters$subject <- sub  
    confidencedata = confidencedata[confidencedata$subject_name == sub,]
  }
  
  #create dataframe to store likelihoods in each rep
  draws = 101
  n_trials = nrow(confidencedata)
  
  #need to set an initial value for this to work
  log.likelihoods <- 0
  
  for (rep in 1:reps) {   
    
    bayesian_model <- function(par){
      
      # || is true if at least one of the conditions is True
      if (any(par[8:12] <= 0) | #none of the sigma parameters can be 0 
          any(par[8:12] > 20) |
          any(diff(round(par[8:11], 1)) >= 0) | #make sure sigma parameters are at least one decimal point different
          any(par[2:7] > 0)) { #make sure boundaries are increasing in orientation space
        log.likelihoods <- Inf
      }
      
      if (log.likelihoods < Inf) {
        
        boundaries = c(0, cumsum(par[1:7]), Inf)
        dnoise = par[12]
        
        #get upper and lower boundary draws for each trial, based on response
        b_low <- lower_boundaries <- sapply(boundaries[confidencedata$r], function(x) get_draws(x, dnoise))
        b_high <- upper_boundaries <- sapply(boundaries[(confidencedata$r + 1)], function(x) get_draws(x, dnoise))
        
        sigma = par[(7+confidencedata$stim_reliability_level)] #get baseline sigma for each reliability level
        
        #transform boundary parameters to orientation space using sigma for that trial
        b_low[,confidencedata$r != 1] = apply(rbind(lower_boundaries[,confidencedata$r != 1], sigma[confidencedata$r != 1]), 2, function(x) d_to_ori(x[1:(length(x)-1)], x[length(x)]))
        b_high[,confidencedata$r != 8] = apply(rbind(upper_boundaries[,confidencedata$r != 8], sigma[confidencedata$r != 8]), 2, function(x) d_to_ori(x[1:(length(x)-1)], x[length(x)]))
        
        if (any(sapply(par[8:11], function(x) any(d_to_ori(par[1:7], x) < 0)))) {
          log.likelihoods <- Inf
        }
        #make sure lower boundary positions are not less than 0 
        #if (any(b_low[51 , (confidencedata$r != 1)] < 0) | #filter for lower trials indexing lower boundary limit (0)
        #any(b_high[51,] < 0) | (any(is.na(b_low))) | (any(is.na(b_high)))) {
        #  log.likelihoods <- Inf
        #}
        
        if (log.likelihoods < Inf) {
          #get density of boundary draws -- using boundary 1 
          #density are normalised (divided by total density) so use density values from boundary 1 for all boundary parameters.
          #different boundaries are just shifts in the mean and d_noise sd stays the same
          density <- dnorm(get_draws(boundaries[2], dnoise), boundaries[2], dnoise)/sum(dnorm(get_draws(boundaries[2], dnoise), boundaries[2], dnoise))
          
          orientations <- matrix(rep(confidencedata$stim_orientation, each = draws), draws, n_trials) #rep orientations so there is run for each draw of the boundary position
          sigmas <- matrix(rep(sigma, each = draws), draws, n_trials) #rep orientations so there is run for each draw of the boundary position
          
          #get likelihood for upper positive boundary 
          likelihoods = density*(((pnorm(b_high, mean = orientations, 
                                         sd = sigmas)) - 
                                    #get likelihood for lower positive boundary 
                                    (pnorm(b_low, mean = orientations, 
                                           sd = sigmas))) + 
                                   #get likelihood for upper negative boundary
                                   #could use sort() here too 
                                   ((pnorm(apply(-b_low, 2, function(x) rev(x)), mean = orientations, 
                                           sd = sigmas)) -
                                      #get likelihood for lower positive boundary 
                                      (pnorm(apply(-b_high, 2, function(x) rev(x)), mean = orientations, 
                                             sd = sigmas))))
          
          #am I allowed to do this?
          #I probably need to change this
          likelihoods[likelihoods == 0] <- min(likelihoods[likelihoods>0])
          
          log.likelihoods = -sum(log(colSums(likelihoods)))
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


