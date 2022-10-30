simulate_exp <- function (sub, parameters, modality = 'visual') {
  
  setwd("/Users/s4323621/Dropbox/Documents/Cluster_Files")
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
  
  #write function to convert from d space to measurement space for each boundary 
  measurement <- function(b, sigma) {
    d = -log((1/b) - 1)
    ori = (d*((sigma^2)+(cat_sigma^2)))/(2*cat_mean)
    return(ori)  
  }
  
  confidencedata$stim_orientation <- scale(confidencedata$stim_orientation)
  confidencedata <- confidencedata[confidencedata$subject_name == sub,]
  
  n_samples <- 101
  
  # Create data storage
  responses <- data.frame(matrix(NA, 1, n_samples))
  
  rel1 = c(-measurement(parameters[1:3], parameters[5]), measurement(parameters[4], parameters[5]), rev(measurement(parameters[1:3], parameters[5])))
  rel2 = c(-measurement(parameters[1:3], parameters[6]), measurement(parameters[4], parameters[6]), rev(measurement(parameters[1:3], parameters[6])))
  rel3 = c(-measurement(parameters[1:3], parameters[7]), measurement(parameters[4], parameters[7]), rev(measurement(parameters[1:3], parameters[7])))
  rel4 = c(-measurement(parameters[1:3], parameters[8]), measurement(parameters[4], parameters[8]), rev(measurement(parameters[1:3], parameters[8])))
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
