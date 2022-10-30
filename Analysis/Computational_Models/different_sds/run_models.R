rm(list = ls())
datas = c("visual", "auditory") #c("visual", "auditory")
models = c("fixed", "lin", "quad", "exp", "bayes", "bayes_prior") #c("fixed", "lin", "quad", "exp", "bayes", "bayes_prior")

setwd("/Users/s4323621/Dropbox/Documents/multimodal_confidence/updated")
n_starting_parameters = 20 #how many sets of starting parameters to use for each model 

#load libraries
library(doParallel)
library(foreach)
library(dplyr)
library(beepr)

#set up parallel processing 
cores = 8
cl <- makeCluster(cores)
registerDoParallel(cl)

for (data in datas) {
##========= LOAD DATA =================
#load data - this won't work for multiple data sets 
if (data == "visual") {
    confidencedata <- read.csv("/Users/s4323621/Dropbox/Documents/multimodal_confidence/data/bimodal/visual_data.csv")
} else if (data == "auditory") {
  confidencedata <- read.csv("/Users/s4323621/Dropbox/Documents/multimodal_confidence/data/bimodal/auditory_data.csv")
}

for (model in models) {
#select model class 
if (model == 'lin' | model == 'quad' | model == 'exp') {
class = "linear"
} else if (model == 'bayes' | model == 'bayes_prior') {
  class = "bayesian"
} else if (model == 'fixed') {
  class = "fixed"
}
modality = data
##========= GENERATE STARTING PARAMETERS =================
#generate starting parameters  - only works with a single model right now 
source(paste0(getwd(), "/generate_", class, "_parameters.R"))
tic = Sys.time()
if (class == "fixed") {
  starting_parameters <- foreach(sub = 1:n_starting_parameters,
                                 .combine = "rbind")  %dopar% {
                                   output = generate_fixed_parameters()
                                 }
} else if (class == "linear") {
  starting_parameters <- foreach(sub = 1:n_starting_parameters,
                                         .combine = "rbind")  %dopar% {
                                           output = generate_linear_parameters(model = model)
                                         }
} else if (class == "bayesian") {
  inital_starting_parameters = c(1,	-0.07,	-0.222,	-0.609,	-1.008,	-1.556,	-1.732,	
                                 4,	2,	1.2,	0.7)
  if (model == "bayes_prior") {
    fit_prior = TRUE
  } else {fit_prior = FALSE}
  starting_parameters <-foreach(sub = 1:n_starting_parameters,
                                         .combine = "rbind")  %dopar% {
                                           output = generate_bayesian_parameters(inital_starting_parameters, #means 
                                                                                 rep(3, 12), #sds
                                                                                 data, 
                                                                                 fit_prior = fit_prior)
                                         }
}
  toc = Sys.time()
  toc - tic
  beep()
  
starting_parameters = starting_parameters[which(!is.na(starting_parameters[,1])),]
n_subjects = length(unique(confidencedata$subject_name))
combinations <- expand.grid(1:nrow(starting_parameters), 1:n_subjects)

##========= RUN MODEL =================
source(paste0(getwd(), "/", class, "_model.R"))

tic = Sys.time()
if (class == "fixed") {
  output <- foreach(x = 1:nrow(combinations),
                    .combine = "rbind")  %dopar% {
                      output = fixed_model_function(confidencedata, combinations[x, 2], #subject_number
                                                    100,  #number of repetitions
                                                    starting_parameters[combinations[x, 1],]) #starting_parameters 
                    }
} else if (class == "linear") {
  output <- foreach(x = 1:nrow(combinations),
                    .combine = "rbind")  %dopar% {
                      output = linear_model_function(confidencedata, sub = combinations[x, 2], #subject_number
                                                    reps = 100,  #number of repetitions
                                                    starting_parameters[combinations[x, 1],], #starting_parameters 
                                                    model = model) 
                    }
} else if (class == "bayesian") {
  output <- foreach(x = 1:nrow(combinations),
                    .combine = "rbind")  %dopar% {
                      output = bayesian_model_function(confidencedata, 
                                                       combinations[x, 2], #give function subject number  
                                                       100,
                                                       starting_parameters[combinations[x, 1],], #starting parameters 
                                                       data, #modality 
                                                       fit_prior = fit_prior) 
                    }
}
toc = Sys.time()
toc - tic
beep()

#get best fitting parameters for each subject 
best_estimate = output %>%
  group_by(subject) %>%
  summarise(index = which.min(log_likelihood), 
            lowest_likelihood = log_likelihood[index])

best_fits = matrix(0, max(output$subject), ncol(output))
for (sub in unique(output$subject)) {
  best_fits[sub,] = as.numeric(output[(output$subject == sub & output$log_likelihood == best_estimate$lowest_likelihood[sub]),])
}
colnames(best_fits) = colnames(output)

write.csv(best_fits, paste0("/Users/s4323621/Dropbox/Documents/multimodal_confidence/parameters/", 
                              data, "_", model, ".csv"))

##=========== SIMULATE DATA==============
source(paste0(getwd(), "/simulate.R"))
simulated_data <- foreach(sub = 1:n_subjects,
                          .combine = "rbind")  %dopar% {
                            output = simulate(sub, 
                                              confidencedata, 
                                              as.numeric(best_fits[sub, 1:ncol(starting_parameters)]), #best fitting parameters
                                              means = TRUE, #do we want the mean response for each trial
                                              model = model,
                                              class = class)
                          }

write.csv(simulated_data, paste0("/Users/s4323621/Dropbox/Documents/multimodal_confidence/simulated_data/", 
                            data, "_", model, ".csv"))
rm(output)
}
}

