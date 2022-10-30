simulate <- function(sub, parameters, modality = 'visual') {
  
  library(dplyr)
  #=======================================================
  #Get the data 
  setwd("/Users/s4323621/Dropbox/Documents/Cluster_Files")
  if (modality == 'visual') {
    confidencedata <- read.csv(file = 'visual_taskA_data.csv')
    cat_mean = -4
    cat_sigma = 5
  }
  if (modality == 'auditory') {
    confidencedata <- read.csv(file = 'visual_taskA_auditory_data.csv')
    cat_mean = 2300
    cat_sigma = 475
  }
  #filter for subject data
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
  b = cumsum(parameters[1:4])
  #THIS IS SO MESSY AND DIFFERENT TO THE MODEL SO DON'T DO IT 
  b_list = rev(c(rev(b[2:4]), b[1], -b[2:4])) #list of boundary parameters
  sigmas <- as.numeric(parameters[5:8]) #sigma parameters
  psi <- as.numeric(parameters[9])
  d_noise = as.numeric(parameters[10]) #decision noise parameter
  
  #=======================================================
  #Let's get down to business  
  #create dataframe with 100 repeats of every trial 
  simulated_data <- confidencedata[rep(seq_len(nrow(confidencedata)), each = 100), ]
  
  baseline_sigmas = sigmas[simulated_data$stim_reliability_level]
  odn = psi*(abs(sin((simulated_data$stim_orientation*pi)/90)))
  sigma_list = baseline_sigmas + odn 
  
  #get draws of psychological orientations 
  simulated_data$psychological_orientation <- rnorm(nrow(simulated_data), mean = simulated_data$stim_orientation, sd = sigma_list)
  
  #get boundary draws
  simulated_data$boundary_draw <- rnorm(nrow(simulated_data), mean = 0, sd = d_noise)
  
  #create matrix with all boundaries and upper/lower limits 
  d_boundaries <- t(matrix(b_list, length(b_list), nrow(simulated_data)))
  
  #get the draw position for each boundary 
  d_boundaries_draws <- d_boundaries + simulated_data$boundary_draw
  d_boundaries_ori = t(apply(cbind(d_boundaries_draws, sigma_list), 1, function(x) measurement(x[1:7], x[8])))
  
  #find where orientation lies between boundaries positions
  line_up = cbind(simulated_data$psychological_orientation, d_boundaries_ori)
  
  #it actually would make a lot of sense to use order() here but I can't figure out how to get it to work
  sorted_line_up = t(apply(line_up, 1, sort)) 
  
  #find position of stimulus orientation and save it as the response 
  indexes = apply(cbind(simulated_data$psychological_orientation, sorted_line_up), 
                  1, function(x) which(x[1] == x[2:9]))
  responses = c(-4, -3, -2, -1, 1, 2, 3, 4) 
  simulated_data$simulated_r = responses[indexes]
  simulated_data$psychological_category = ((simulated_data$simulated_r > 4) + 1)
  simulated_data$psychological_confidence = 0
  simulated_data$psychological_confidence[simulated_data$psychological_category == 1] = 5 - simulated_data$simulated_r[simulated_data$psychological_category == 1]
  simulated_data$psychological_confidence[simulated_data$psychological_category == 2] = simulated_data$simulated_r[simulated_data$psychological_category == 2] - 4
  
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