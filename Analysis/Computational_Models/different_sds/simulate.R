simulate <- function(sub, confidencedata, parameters, modality = 'visual', 
                     means = TRUE, model, class)  {
  #Rebecca West 2022
  #This function simulates data for all models 

  ##============Determine what data to simulate ===========
  #confidencedata$stim_orientation = scale(confidencedata$stim_orientation)
  if (all(confidencedata$stim_orientation == 0)) {
    confidencedata$stim_orientation = confidencedata$stim_frequency #hacky way of dealing with auditory data 
  }
  
  #if fitting group data 
  if (sub != 0) {
    confidencedata = confidencedata[confidencedata$subject_name == sub,]
  }
  
##===========LOAD LIBRARIES=============  
  library(dplyr) 

##============ FUNCTIONS================
  #bayesian function
  d_to_ori <- function(d, s, cat_one_sigma, cat_two_sigma) {
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
  
  #scaled evidence stregnth 
  boundary <- function(k, m, s, exp) {
    b = k + (m*(s^exp))
    return(b)
  }
  
##=========== GET PARAMETERS =================  
  #create dataframe with 100 repeats of every trial 
  confidencedata$trial = 1:nrow(confidencedata)
  simulated_data <- confidencedata[rep(seq_len(nrow(confidencedata)), each = 100), ]
  
  if (class == "bayesian") {
  b_list <- cumsum(as.numeric(parameters[1:7])) #list of boundary parameters
  sigmas <- as.numeric(parameters[8:11]) #sigma parameters
  sigma_list = sigmas[simulated_data$stim_reliability_level]
    if (model == 'bayes_prior') {
      cat_one_sigma <- parameters[12]
      cat_two_sigma <- parameters[13]
    } else {
      if (modality == 'visual') {
        cat_one_sigma = 0.3427148
        cat_two_sigma = 1.370859
      } else if (modality == 'auditory') {
        cat_one_sigma = 0.3429777
        cat_two_sigma = 1.371911
      }
    }
    d_boundaries <- t(matrix(b_list, length(b_list), nrow(simulated_data)))
    
    #transform boundary positions into orientation space 
    positive_ori_boundaries <- t(apply(cbind(d_boundaries, sigma_list), 1, function(x) d_to_ori(x[1:7], x[8], cat_one_sigma, cat_two_sigma)))
    negative_ori_boundaries <- t(apply(cbind(d_boundaries, sigma_list), 1, function(x) rev(-d_to_ori(x[1:7], x[8], cat_one_sigma, cat_two_sigma))))
    all_boundaries = cbind(negative_ori_boundaries, positive_ori_boundaries)
} else if (class == "linear") {
  ks = parameters[1:7]
  ms = parameters[8:14]
  sigmas <- parameters[15:18] 
  sigma_list = sigmas[simulated_data$stim_reliability_level]
  
  if (model == 'lin') {
    exponent = 1
  } else if (model == 'quad') {
    exponent = 2
  } else if (model == 'exp') {
    exponent = parameters[19]
  }
    #get boundary paramters 
    all_boundaries = t(sapply(sigma_list, function(x) c(-rev(boundary(ks, ms, x, exponent)), 
                                                        boundary(ks, ms, x, exponent))))  
} else if (class == "fixed") {
  b_list <- as.numeric(parameters[1:7]) #list of boundary parameters
  sigmas <- as.numeric(parameters[8:11]) #sigma parameters
  sigma_list = sigmas[simulated_data$stim_reliability_level]
    #get boundary paramters 
    boundaries = c(-rev(b_list), b_list)
    all_boundaries = t(matrix(boundaries, length(boundaries), nrow(simulated_data)))
}
  
  #get draws of psychological orientations 
  simulated_data$psychological_orientation <- rnorm(nrow(simulated_data), mean = simulated_data$stim_orientation, sd = sigma_list)
  line_up = cbind(simulated_data$psychological_orientation, all_boundaries)
  #sort boundaries and psychological orientation 
  sorted_line_up = t(apply(line_up, 1, sort)) 
  
  #find position of stimulus orientation and save it as the response 
  indexes = apply(cbind(simulated_data$psychological_orientation, sorted_line_up), 
                    1, function(x) which(x[1] == x[2:16]))
  #responses = c("8", "7", "6", "5", "4", "3", "2", "1", "2", "3", "4", "5", "6", "7", "8")
  responses = c("4", "3", "2", "1", "-1", "-2", "-3", "-4", "-3", "-2", "-1", "1", "2", "3", "4")
  simulated_data$simulated_r = as.numeric(responses[indexes])
  
  #get mean response for each trial 
  if (means == TRUE) {
    summarised_simulated_data <- simulated_data %>%
      group_by(trial) %>%
      summarise(subject_name = unique(subject_name),
                stim_type = unique(stim_type),
                task_type = unique(task_type),
                resp_category = mean(resp_category),
                resp_confidence = mean(resp_confidence),
                resp_correct = mean(resp_correct),
                stim_orientation = mean(stim_orientation),
                stim_frequency = mean(stim_frequency),
                stim_reliability = mean(stim_reliability),
                stim_category = mean(stim_category),
                stim_reliability_level = mean(stim_reliability_level),
                #ideal_accuracy = mean(ideal_accuracy),
                r = mean(r),
                psychological_r = mean(simulated_r),
                psychological_orienation = mean(psychological_orientation),
                psychological_confidence = mean(abs(simulated_r)),
                psychological_category = mean((simulated_r > 0) + 1))
    return(summarised_simulated_data)
  } else {
    return(simulated_data)    
  }
}