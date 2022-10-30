library(doParallel)
library(foreach)
library(dplyr)

#source model function 
source("/Users/s4323621/Dropbox/Documents/EndofDay_Confidence_Models/Exp_ModelA/generate_exp_parametersA.R")
source("/Users/rebeccawest/Dropbox/Documents/EndofDay_Confidence_Models/Exp_ModelA/generate_exp_parametersA.R")

cores = 8
cl <- makeCluster(cores)
registerDoParallel(cl)


#set number of sets of starting_parameters for each subject
n_starting_parameters = 20 #3

tic = Sys.time()
starting_parameters <- foreach(x = 1:n_starting_parameters,
                               .combine = "rbind")  %dopar% {
                                 output = generate_exp_parameters()
                               }
toc = Sys.time()
toc - tic

successes = length(starting_parameters[starting_parameters != 0])/13
starting_parameters = matrix(starting_parameters[starting_parameters != 0], successes, 13)

n_starting_parameters = nrow(starting_parameters) #3
n_subjects = 8 #number of subjects
combinations <- expand.grid(1:n_starting_parameters, 1:n_subjects)


source("/Users/rebeccawest/Dropbox/Documents/EndofDay_Confidence_Models/Exp_ModelA/exp_model_taskA.R")

tic=Sys.time()
output <- foreach(x = 1:nrow(combinations),
                  .combine = "rbind")  %dopar% {
                    output = exp_model_function(combinations[x, 2], #subject_number
                                                 100,  #number of repetitions
                                                 starting_parameters[combinations[x, 1],],
                                                 'auditory') #starting_parameters 
                  }
toc=Sys.time()
print(toc-tic)

#fitting the group data
#tic=Sys.time()
#output <- foreach(x = 1:nrow(starting_parameters),
#                  .combine = "rbind")  %dopar% {
#                    output = quadratic_model_function(1:7, #subject_number
#                                                  100,  #number of repetitions
#                                                 starting_parameters[x,]) #starting_parameters 
#                  }
#toc=Sys.time()
#print(toc-tic)

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


source("/Users/s4323621/Dropbox/Documents/EndofDay_Confidence_Models/Exp_ModelA/simulate_exp_data_taskA.R")
best_fits <- read.csv('/Users/s4323621/Dropbox/Documents/Best_Fitting_Parameters/TaskA/Exp/auditory.csv')

tic=Sys.time()
print("running loops")
simulated_data <- foreach(sub = 1:8,
                          .combine = "rbind")  %dopar% {
                            output = simulate(sub, as.numeric(best_fits[sub, 1:13]),
                                              'auditory')
                          }
toc=Sys.time()
print(toc-tic)


write.csv(output, "/Users/s4323621/Dropbox/Documents/EndofDay_Confidence_Models/Linear_ModelA/standardised_linear_model_visual_output_.csv")
