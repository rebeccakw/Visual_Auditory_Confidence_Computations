simulate <- function (sub, parameters) {

  library(dplyr)
  setwd("/Users/s4323621/Dropbox/Documents/Cluster_Files")
  visual <- read.csv(file = 'visual_taskA_data.csv')
  auditory <- read.csv(file = 'visual_taskA_auditory_data.csv')
  
  visual$stim_orientation <- scale(visual$stim_orientation)
  auditory$stim_orientation <- scale(auditory$stim_orientation)
  confidencedata <- rbind(visual, auditory)
  
  confidencedata <- confidencedata[confidencedata$subject_name == sub,]
  
  boundary_function <- function(k, m, sigma, exp) {
    b = k + (m*(sigma^exp))
    return(b)
  }
 
  #Rename parameters so it's easier to follow  
  k_params = parameters[1:4]
  m_params = parameters[5:8]
  sigmas <- as.numeric(parameters[9:12]) #sigma parameters
  exp = as.numeric(parameters[13]) #decision noise parameter
  scaler = as.numeric(parameters[14])
  
  confidencedata$trial = 1:nrow(confidencedata)
  #=======================================================
  #Let's get down to business  
  #create dataframe with 100 repeats of every trial 
  simulated_data <- confidencedata[rep(seq_len(nrow(confidencedata)), each = 100), ]
  
  sigma_list = sigmas[simulated_data$stim_reliability_level]
  sigma_list[simulated_data$stim_type == 'Auditory'] = scaler*sigma_list[simulated_data$stim_type == 'Auditory']
  
  #get draws of psychological orientations 
  simulated_data$psychological_orientation <- rnorm(nrow(simulated_data), mean = simulated_data$stim_orientation, sd = sigma_list)
  
  boundaries <- matrix(rep(c(k_params, m_params), each = nrow(simulated_data)), nrow = nrow(simulated_data), ncol = length(c(k_params, m_params)))
  
  boundaries_pos = t(apply(cbind(boundaries, sigma_list), 1, function(x) boundary_function(x[1:4], x[5:8], x[9], exp)))
  boundaries_neg = t(apply(cbind(boundaries, sigma_list), 1, function(x) -rev(boundary_function(x[2:4], x[6:8], x[9], exp))))
  
  boundaries_all = cbind(boundaries_neg, boundaries_pos)
  #find where orientation lies between boundaries positions
  line_up = cbind(simulated_data$psychological_orientation, boundaries_all)
  
  #it actually would make a lot of sense to use order() here but I can't figure out how to get it to work
  sorted_line_up = t(apply(line_up, 1, sort)) 
  
  #find position of stimulus orientation and save it as the response 
  indexes = apply(cbind(simulated_data$psychological_orientation, sorted_line_up), 
                  1, function(x) which(x[1] == x[2:9]))
  responses = c(-4, -3, -2, -1, 1, 2, 3, 4) 
  simulated_data$simulated_r = responses[indexes]
  simulated_data$psychological_category = ((simulated_data$simulated_r > 0) + 1)
  simulated_data$psychological_confidence = abs(simulated_data$simulated_r)
  
  simulated_data = simulated_data %>%
    group_by(trial) %>%
    summarise(stim_type = unique(stim_type),
      psychological_category = mean(psychological_category),
              psychological_confidence = mean(psychological_confidence),
              psychological_orientation = mean(psychological_orientation),
              psychological_r = mean(simulated_r),
              stim_orientation = mean(stim_orientation),
              stim_reliability_level = mean(stim_reliability_level),
              resp_category = mean(resp_category),
              resp_confidence = mean(resp_confidence),
              subject_name = mean(subject_name),
              stim_category = mean(stim_category))
  
  return(simulated_data)
}

  