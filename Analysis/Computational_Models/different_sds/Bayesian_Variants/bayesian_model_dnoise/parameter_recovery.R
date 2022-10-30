#Parameter Recovery
library(doParallel)
library(foreach)
library(dplyr)
library(beepr)

cores = 8
cl <- makeCluster(cores)
registerDoParallel(cl)

##=========== Define Parameters ====================
source("/Users/s4323621/Dropbox/Documents/multimodal_confidence/bayesian_model_dnoise/generate_bayesian_parameters_dnoiseonly.R")

confidencedata <- read.csv("/Users/s4323621/Dropbox/Documents/multimodal_confidence/data/visual_taskB.csv")

#set number of sets of starting_parameters for each subject
n_subjects = 8
parameters_per_subject = 3
n_starting_parameters = n_subjects*parameters_per_subject #3
inital_starting_parameters = c(0.4,	-0.07,	-0.222,	-0.609,	-1.008,	-1.556,	-1.732,	
                               1.5,	1.2,	1,	0.7,	0.5)

#from standard bayesian model
inital_starting_parameters = c( 0.358,    -0.012,    -0.064,    -0.266,    -0.323,    -0.398,    -0.459,     
                                1.2,     1.1,     1,
                                0.9, 0.5)


tic = Sys.time()
starting_parameters <- foreach(sub = 1:n_starting_parameters,
                               .combine = "rbind")  %dopar% {
                                 output = generate_bayesian_parameters(inital_starting_parameters, 
                                                                           rep(1, 12), 'visual')
                               }
toc = Sys.time()
toc - tic
beep()


starting_parameters = starting_parameters[which(!is.na(starting_parameters[,1])),]
n_starting_parameters = nrow(starting_parameters) #3
combinations = cbind(rep(1:n_subjects, each = parameters_per_subject), starting_parameters)

##====================== Load Experimental Data ==================
confidencedata <- read.csv("/Users/rebeccawest/Dropbox/Documents/multimodal_confidence/data/visual_taskB.csv")
confidencedata$stim_orientation <- scale(confidencedata$stim_orientation)


data_size = 720*0.5
confidencedata <- data.frame(matrix(vector(), data_size*nrow(starting_parameters), 4, 
                                    dimnames = list(c(), c("stim_orientation", "stim_reliability_level",
                                                           "r", "simulation_set"))))
confidencedata$simulation_set = rep(1:nrow(starting_parameters), each = data_size)
for (set in 1:nrow(starting_parameters)) {
  confidencedata$stim_reliability_level[confidencedata$simulation_set == set] = rep(1:4, data_size/4)
  confidencedata$stim_orientation[confidencedata$simulation_set == set] = c(rnorm(data_size/2, 0, 3), 
                                                                            rnorm(data_size/2, 0, 12))
}

#sequenced data 
#confidencedata$simulation_set = rep(1:nrow(starting_parameters), each = data_size)
#for (set in 1:nrow(starting_parameters)) {
#  confidencedata$stim_reliability_level[confidencedata$simulation_set == set] = rep(1:4, data_size/4)
#  confidencedata$stim_orientation[confidencedata$simulation_set == set] = c(seq(-12, 12, length.out = data_size/2), 
#                                                                            seq(-12, 12, length.out = data_size/2))
#}
##============= Simulate Data =====================
source("/Users/s4323621/Dropbox/Documents/multimodal_confidence/bayesian_model_dnoise/simulate_bayesian_data_dnoise.R")
#use this for simulated data 
tic = Sys.time()
simulated_data <- foreach(x = 1:nrow(starting_parameters),
                          .combine = "rbind")  %dopar% {
                            output = simulate(
                              0, #subject number
                              confidencedata[confidencedata$simulation_set == x,],
                              as.numeric(starting_parameters[x,]),
                              means = FALSE) #starting_parameters
                          }
toc=Sys.time()
print(toc-tic)
beep()

colnames(simulated_data)[7] <- "simulation_set"
colnames(simulated_data)[17] <- "simulation_set"
simulated_data$r <- simulated_data$simulated_r

#you need this if the means argument is true 
#simulated_data$resp_confidence <- round(simulated_data$psychological_confidence)
#simulated_data$resp_category <- round(simulated_data$psychological_category)
#simulated_data$r[simulated_data$resp_category == 1] <- 5 - simulated_data$resp_confidence[simulated_data$resp_category == 1]
#simulated_data$r[simulated_data$resp_category == 2] <- simulated_data$resp_confidence[simulated_data$resp_category == 2] + 4

##================ Fit Model to Data ==================
source("/Users/s4323621/Dropbox/Documents/multimodal_confidence/bayesian_model_dnoise/bayesian_model_dnoise.R")
n_parameters = ncol(starting_parameters)
print("running loops")
output <- foreach(x = 1:nrow(combinations),
                  .combine = "rbind")  %dopar% {
                    output = bayesian_model(simulated_data[simulated_data$simulation_set == x,], 
                                                     #as.numeric(combinations[x, 1]), #subject_number
                                                     1:8,
                                                     100,  #number of repetitions
                                                     as.numeric(combinations[x, 2:(n_parameters+1)]),
                                                     'visual') #starting_parameters 
                  }
toc=Sys.time()
print(toc-tic)
beep()

parameter_recovery = cbind(stack(output[,1:n_parameters]), stack(data.frame(starting_parameters)))
parameter_recovery$set = rep(1:nrow(starting_parameters), each = n_parameters)
colnames(parameter_recovery) <- c("fitted", "parameter", "generating", "NA", "set")

write.csv(parameter_recovery, "/Users/s4323621/Dropbox/Documents/multimodal_confidence/recovery_data/bayesian_dnoise_720*0.5.csv")

ggplot(parameter_recovery) + geom_point(aes(generating, fitted)) + 
  facet_wrap_equal(~parameter, scales = "free", ncol = 3) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2) + xlab("Simulated Parameter") + 
  ylab("Recovered Parameter") + theme(strip.text.x = element_text(size = 20), axis.text=element_text(size=14),
                                      axis.title=element_text(size=25))
