simulate <- function(sub, confidencedata, parameters, modality = 'visual', means = TRUE) {
 
  library(dplyr) 
  #confidencedata$stim_orientation = scale(confidencedata$stim_orientation)

  if (sub != 0) {
    confidencedata <- confidencedata[confidencedata$subject_name == sub,]
  } 
  
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
  
  b_list <- cumsum(as.numeric(parameters[1:7])) #list of boundary parameters
  sigmas <- as.numeric(parameters[8:11]) #sigma parameters
  d_noise = as.numeric(parameters[12])
  
  #create dataframe with 100 repeats of every trial 
  confidencedata$trial = 1:nrow(confidencedata)
  simulated_data <- confidencedata[rep(seq_len(nrow(confidencedata)), each = 1), ]
  
  sigma_list = sigmas[simulated_data$stim_reliability_level]

  #get draws of psychological orientations 
  simulated_data$psychological_orientation <- rnorm(nrow(simulated_data), mean = simulated_data$stim_orientation, sd = sigma_list)
  
  #get boundary draws
  simulated_data$boundary_draw <- rnorm(nrow(simulated_data), mean = 0, sd = d_noise)
  
  #create matrix with all boundaries and upper/lower limits 
  d_boundaries <- t(matrix(b_list, length(b_list), nrow(simulated_data)))
  
  #get the draw position for each boundary 
  d_boundaries_draws <- d_boundaries + simulated_data$boundary_draw
  
  #transform boundary positions into orientation space 
  pos_ori_boundaries <- t(apply(cbind(d_boundaries_draws, sigma_list), 1, function(x) d_to_ori(x[1:7], x[8])))
  neg_ori_boundaries <- t(apply(cbind(d_boundaries_draws, sigma_list), 1, function(x) -rev(d_to_ori(x[1:7], x[8]))))
  
  all_boundaries = cbind(neg_ori_boundaries, pos_ori_boundaries)
  
  #find where orientation lies between boundaries positions
  line_up = cbind(simulated_data$psychological_orientation, all_boundaries)
  
  #it actually would make a lot of sense to use order() here but I can't figure out how to get it to work
  sorted_line_up = t(apply(line_up, 1, sort)) 
  
  #find position of stimulus orientation and save it as the response 
  indexes = apply(cbind(simulated_data$psychological_orientation, sorted_line_up), 
                                     1, function(x) which(x[1] == x[2:16]))
  responses = c("8", "7", "6", "5", "4", "3", "2", "1", "2", "3", "4", "5", "6", "7", "8")
  simulated_data$simulated_r = as.numeric(responses[indexes])
  
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
                stim_reliability = mean(stim_reliability),
                stim_category = mean(stim_category),
                stim_reliability_level = mean(stim_reliability_level),
                ideal_accuracy = mean(ideal_accuracy),
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