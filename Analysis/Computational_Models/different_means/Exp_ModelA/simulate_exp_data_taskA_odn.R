simulate <- function (confidencedata, sub, parameters) {
  
  library(dplyr)
  boundary_function <- function(k, m, s, exp) {
    b = k + (m*(s^exp))
    return(b)
  }
  
  confidencedata$unscaled_stim_orientation = confidencedata$stim_orientation
  
  sample_mean = mean(confidencedata$stim_orientation)
  sample_sd = sd(confidencedata$stim_orientation)
  #use this for odn noise 
  z_score <- function(value) {
    z = (value - sample_mean)/sample_sd
    return(z)
  }
  
  confidencedata$stim_orientation = scale(confidencedata$stim_orientation)
  if (sub != 0) {
    confidencedata <- confidencedata[confidencedata$subject_name == sub,]
  } 
  
  #create dataframe with 100 repeats of every trial 
  confidencedata$trial = 1:nrow(confidencedata)
  simulated_data <- confidencedata[rep(seq_len(nrow(confidencedata)), each = 100), ]
  
  #Rename parameters so it's easier to follow  
  k_params = parameters[1:4]
  m_params = parameters[5:8]
  sigmas <- as.numeric(parameters[9:12]) #sigma parameters
  exp = as.numeric(parameters[13]) 
  
  baseline_sigma = sigmas[simulated_data$stim_reliability_level]
  odn = z_score(as.numeric(parameters[14])*
                  (abs(sin((simulated_data$unscaled_stim_orientation*pi)/90))))
  sigma_list = baseline_sigma + odn

  #calculate boundaries for all trials 
  boundaries = t(sapply(sigma_list, function(x) boundary_function(k_params, m_params, x, exp)))
  #get upper and lower boundaries 
  all_boundaries = cbind(-boundaries[,1:3], #make first 3 boundaries negative
                         boundaries[,4], #category boundary
                         t(apply(boundaries[,1:3], 1, function(x) rev(x)))) #positive boundaries
  
  #get draws of psychological orientations 
  simulated_data$psychological_orientation <- rnorm(nrow(simulated_data), mean = simulated_data$stim_orientation, sd = sigma_list)
  
  #find where orientation lies between boundaries positions
  line_up = cbind(simulated_data$psychological_orientation, all_boundaries)
  
  #it actually would make a lot of sense to use order() here but I can't figure out how to get it to work
  sorted_line_up = t(apply(line_up, 1, sort)) 
  
  #find position of stimulus orientation and save it as the response 
  indexes = apply(cbind(simulated_data$psychological_orientation, sorted_line_up), 
                  1, function(x) which(x[1] == x[2:9]))
  #responses = c("1", "2", "3", "4", "5", "6", "7", "8")
  responses = c(-4, -3, -2, -1, 1, 2, 3, 4) 
  simulated_data$simulated_r = as.numeric(responses[indexes])
  simulated_data$psychological_category = ((simulated_data$simulated_r > 0) + 1)
  simulated_data$psychological_confidence = abs(simulated_data$simulated_r)

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
}
