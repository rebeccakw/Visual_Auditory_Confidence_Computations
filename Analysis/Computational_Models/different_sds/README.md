# **Different SDs Models**

* run_models.R - file used to run all models. Working directory will need to be updated.  
 
## Models 
Log-likelihood functions for all models 
* bayesian_model.R - Log-likelihood function for standard Bayesian model and Bayesian model with free category distrubition parameters. 
* fixed_model.R - Log-likelihood function for distance model (unscaled evidence strength model)
* linear_model.R - Log-likelihood function for linear, quadratic and free-exponent models

## Starting Values 
These functions randomly sample parameters values from a defined normal distribution and ensure that the sampled parameters meet the constraints of each model. Used as starting values for log-likelihood functions. See paper for model constraints. 
* generate_bayesian_parameters.R - for standard bayesian model and bayesian model with free category distribution parameters 
* generate_fixed_parameters.R - for distance model 
* generate_linear_parameters.R - for linear, quadratic and free-exponent model 

## Simulate Model Responses 
Simulates model responses for all models (takes a set of parameters as input)
* simulate.R 

## Variants of Bayesian Model 
Have their own folders because they take longer to run. Each folder has it's own run, sample parameters, calculate log-likelihood & simulate responses file. 
* Bayesian model with decision noise 
* Bayesian model with orientation dependent noise 