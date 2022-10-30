simulate <- function(sub, parameters, modality = 'visual') {
  
  library(dplyr)
  #=======================================================
  setwd("/Users/s4323621/Dropbox/Documents/Cluster_Files")
  #Get the data 
  if (modality == 'visual') {
    confidencedata <- read.csv(file = 'visual_taskA_data.csv')
    cat_mean = -0.6247561
    cat_sigma = 0.7809452
  } 
  if (modality == 'auditory') {
    confidencedata <- read.csv(file = 'visual_taskA_auditory_data.csv')
    cat_mean = -0.6449836
    cat_sigma = 0.7608733
  }
 
  #filter for subject data
  confidencedata$stim_orientation = scale(confidencedata$stim_orientation)
  confidencedata <- confidencedata[confidencedata$subject_name == sub,]
  confidencedata$trial = 1:nrow(confidencedata)
  
  #=======================================================
  #Relevant functions 
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
  #=======================================================
  #Rename parameters so it's easier to follow  
  b_list = cumsum(parameters[1:4])
  sigmas <- as.numeric(parameters[5:8]) #sigma parameters
  d_noise = as.numeric(parameters[9]) #decision noise parameter
  
  #=======================================================
  #Let's get down to business  
  #create dataframe with 100 repeats of every trial 
  simulated_data <- confidencedata[rep(seq_len(nrow(confidencedata)), each = 100), ]
  
  sigma_list = sigmas[simulated_data$stim_reliability_level]
  
  #get draws of psychological orientations 
  simulated_data$psychological_orientation <- rnorm(nrow(simulated_data), mean = simulated_data$stim_orientation, sd = sigma_list)
  
  #get boundary draws
  simulated_data$boundary_draw <- rnorm(nrow(simulated_data), mean = 0, sd = d_noise)
  
  #create matrix with all boundaries and upper/lower limits 
  d_boundaries <- t(matrix(b_list, length(b_list), nrow(simulated_data)))
  
  #get the draw position for each boundary 
  d_boundaries_draws <- d_boundaries + simulated_data$boundary_draw
  d_boundaries_ori_pos = t(apply(cbind(d_boundaries_draws, sigma_list), 1, function(x) measurement(x[1:4], x[5])))
  d_boundaries_ori_neg = t(apply(-d_boundaries_ori_pos[,2:4], 1, function(x) rev(x)))
  
  d_boundaries_ori = cbind(d_boundaries_ori_neg, d_boundaries_ori_pos)
  #find where orientation lies between boundaries positions
  line_up = cbind(simulated_data$psychological_orientation, d_boundaries_ori)
  
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
    summarise(psychological_category = mean(psychological_category),
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