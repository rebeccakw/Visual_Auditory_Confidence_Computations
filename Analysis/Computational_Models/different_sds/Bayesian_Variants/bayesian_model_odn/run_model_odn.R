library(doParallel)
library(foreach)
library(beepr)

wd <- "/Users/rebeccawest/Dropbox/Documents/EndofDay_Confidence_Models/Bayesian_Model/Standard_odn"
setwd(wd)
#source(paste0(wd, "/generate_bayesian_parameters.R"))
source(paste0(wd, "/generate_bayesian_parameters_odn.R"))

#read in dataset 
confidencedata <- read.csv("/Users/rebeccawest/Dropbox/Documents/EndofDay_Confidence_Models/control_exp_data.csv")
confidencedata <- read.csv("/Users/rebeccawest/Dropbox/Documents/Cluster_Files/visual_taskB_data.csv")
#confidencedata$stim_orientation = scale(confidencedata$stim_orientation)

cores = 8
cl <- makeCluster(cores)
registerDoParallel(cl)

#set number of sets of starting_parameters for each subject
n_starting_parameters = 8 #3
n_subjects = 8 #number of subjects
combinations <- expand.grid(1:n_starting_parameters, 1:n_subjects)

#inital_starting_parameters = c(1,	-0.07,	-0.222,	-0.609,	-1.008,	-1.556,	-1.732,	
#                               4,	2,	1.2,	0.7)
#inital_starting_parameters = c(0.7,	-0.07,	-0.222,	-0.609,	-1.008,	-1.556,	-1.732,	
#                               2.5,	2,	1.2,	0.7)
inital_starting_parameters = c(0.7,	-0.07,	-0.222,	-0.609,	-1.008,	-1.556,	-1.732,	
                               2.4,	1.9,	1.1,	0.6, 0.5)

tic = Sys.time()
starting_parameters <- foreach(sub = 1:n_starting_parameters,
                               .combine = "rbind")  %dopar% {
                                 output = generate_bayesian_parameters_odn(inital_starting_parameters, 
                                                                       rep(2, 12), confidencedata,
                                                                       'visual')
                               }
toc = Sys.time()
toc - tic
beep()

starting_parameters = starting_parameters[which(!is.na(starting_parameters[,1])),]
combinations <- expand.grid(1:nrow(starting_parameters), 1:n_subjects)


source(paste0(wd, "/bayesian_model_odn.R"))
tic = Sys.time()
output <- foreach(x = 1:nrow(combinations),
                  .combine = "rbind")  %dopar% {
                    output = bayesian_model_odn(confidencedata, #get data for subject
                                                starting_parameters[combinations[x, 1],], #get starting parameters 
                                                combinations[x, 2], #give function subject number
                                                100, #number of times to run the optimisation
                                                wd, 'visual') 
                  }

toc = Sys.time()
toc - tic
beep()

#fit the group data 
tic = Sys.time()
output <- foreach(x = 1:nrow(starting_parameters),
                  .combine = "rbind")  %dopar% {
                    output = bayesian_model(confidencedata, #get data for subject
                                            starting_parameters[x,], #get starting parameters 
                                           1:8, #give function subject number
                                            100, #number of times to run the optimisation
                                            wd, 'visual') 
                  }

toc = Sys.time()
toc - tic
beep()


library(dplyr)
best_estimate = output %>%
  group_by(subject) %>%
  filter(log_likelihood != 0) %>%
  summarise(index = which.min(log_likelihood), 
            lowest_likelihood = log_likelihood[index])


best_fits = matrix(0, max(output$subject), ncol(output))
for (sub in unique(output$subject)) {
  best_fits[sub,] = as.numeric(output[(output$subject == sub & output$log_likelihood == best_estimate$lowest_likelihood[sub]),])
}
colnames(best_fits) = colnames(output)

#save model output 
write.csv(output, paste0(wd, "/bayesian_model_output_with_psi", toc, ".csv"))

source(paste0(wd,"/bayesian_model_simulations_odn.R"))

simulated_data <- foreach(sub = 1:n_subjects,
                          .combine = "rbind")  %dopar% {
                            output = simulate_data(confidencedata, 
                                                   as.numeric(best_fits[sub, 1:12]), 
                                                   sub, 'visual')
                          }

#for group data 
#simulated_data = simulate_data(confidencedata, best_fits[,1:11], 1:8)

source("/Users/s4323621/Dropbox/Documents/EndofDay_Confidence_Models/Bayesian_Model/simulate_bayesian_data.R")

