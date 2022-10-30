simulate <- function (sub, parameters, modality = 'visual') {

  setwd("/Users/rebeccawest/Dropbox/Documents/Cluster_Files")
  if (modality == 'visual') {
    confidencedata <- read.csv(file = 'visual_taskA_data.csv')
  } 
  if (modality == 'auditory') {
    confidencedata <- read.csv(file = 'visual_taskA_auditory_data.csv')
  }
  
  confidencedata$stim_orientation <- scale(confidencedata$stim_orientation)
  
  confidencedata <- confidencedata[confidencedata$subject_name == sub,]
  
  n_samples <- 101
  
  # Create data storage
  responses <- data.frame(matrix(NA, 1, n_samples))
  
  for (trial in 1:nrow(confidencedata)) {
    rel = confidencedata$stim_reliability_level[trial]
    ori = confidencedata$stim_orientation[trial]
    sigma <- parameters[(4+rel)] 
    measurements <- rnorm(n_samples, mean = ori, 
                          sd = sigma)
    
    boundaries <- c(-parameters[1:3], parameters[4], rev(parameters[1:3]))
    
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
