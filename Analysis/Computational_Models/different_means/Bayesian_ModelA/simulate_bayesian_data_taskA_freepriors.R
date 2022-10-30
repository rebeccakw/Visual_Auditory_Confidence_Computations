simulate <- function (sub, parameters, modality = 'visual') {
  
  setwd("/Users/rebeccawest/Dropbox/Documents/Cluster_Files")
  
  if (modality == 'visual') {
    confidencedata <- read.csv(file = 'visual_taskA_data.csv')
  }
  
  if (modality == 'auditory') {
    confidencedata <- read.csv(file = 'visual_taskA_auditory_data.csv')
  }
  
  #write function to convert from d space to measurement space for each boundary 
  measurement <- function(d, sigma, cat_mean, cat_sigma)  {
    measurement = (d*((sigma^2)+(cat_sigma^2)))/(2*cat_mean)
    return(measurement)  
  }
  
  confidencedata$stim_orientation <- scale(confidencedata$stim_orientation)
  confidencedata <- confidencedata[confidencedata$subject_name == sub,]
  
  n_samples <- 101
  
  # Create data storage
  cat_mean = parameters[9]
  cat_sigma = parameters[10]
  responses <- data.frame(matrix(NA, 1, n_samples))
  
  rel1 = c(-measurement(parameters[1:3], parameters[5], cat_mean, cat_sigma), measurement(parameters[4], parameters[5], cat_mean, cat_sigma), rev(measurement(parameters[1:3], parameters[5], cat_mean, cat_sigma)))
  rel2 = c(-measurement(parameters[1:3], parameters[6], cat_mean, cat_sigma), measurement(parameters[4], parameters[6], cat_mean, cat_sigma), rev(measurement(parameters[1:3], parameters[6], cat_mean, cat_sigma)))
  rel3 = c(-measurement(parameters[1:3], parameters[7], cat_mean, cat_sigma), measurement(parameters[4], parameters[7], cat_mean, cat_sigma), rev(measurement(parameters[1:3], parameters[7], cat_mean, cat_sigma)))
  rel4 = c(-measurement(parameters[1:3], parameters[8], cat_mean, cat_sigma), measurement(parameters[4], parameters[8], cat_mean, cat_sigma), rev(measurement(parameters[1:3], parameters[8], cat_mean, cat_sigma)))
  boundaries_all = rbind(rel1, rel2, rel3, rel4) 
  
  for (trial in 1:nrow(confidencedata)) {
    rel <- confidencedata$stim_reliability_level[trial]
    ori <-  confidencedata$stim_orientation[trial]
    baseline_sigma <- parameters[(4+rel)] 
    #odn <- par[9]*(abs(sin((ori*pi)/90)))
    sigma <-  baseline_sigma #+ odn
    measurements <- rnorm(n_samples, mean = ori, 
                          sd = sigma)
    
    boundaries <- boundaries_all[rel,]
    
    for (measurement in 1:n_samples) {
      responses[measurements[measurement] < boundaries[1], measurement] <- -4
      responses[(measurements[measurement] > boundaries[1] &  measurements[measurement] < boundaries[2]), measurement] <- -3
      responses[(measurements[measurement] > boundaries[2] &  measurements[measurement] < boundaries[3]), measurement] <- -2
      responses[(measurements[measurement] > boundaries[3] &  measurements[measurement] < boundaries[4]), measurement] <- -1
      responses[(measurements[measurement] > boundaries[4] &  measurements[measurement] < boundaries[5]), measurement] <- 1
      responses[(measurements[measurement] > boundaries[5] &  measurements[measurement] < boundaries[6]), measurement] <- 2
      responses[(measurements[measurement] > boundaries[6] &  measurements[measurement] < boundaries[7]), measurement] <- 3
      responses[(measurements[measurement] > boundaries[7]), measurement] <- 4
    }
    
    confidencedata$simulated_resp_buttonid[trial] = mean(as.matrix(responses))
    confidencedata$psychological_confidence[trial] = mean(as.matrix(abs(responses)))
    confidencedata$psychological_category[trial] = mean(1 + (responses > 0))
    confidencedata$psychological_orientation[trial] = mean(measurements)
  }
  
  return(confidencedata)
}
