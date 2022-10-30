generate_linear_parameters <- function(model) {
#============================
#Rebecca West 2022
#This function randomly generates parameters for the scaled evidence strength model class.
#Models in this class: Linear, Quadratic & Free-Exponent
#Parameters include: ks + ms (boundaries), sigmas and exponent.
#Parameters are sampled from a normal distribution with a defined mean (referred to with mu) 
#and standard deviation (referred to as sd). All parameters also have a lower and upper bound. 
  
#Function continues to resample parameters from the define distribution until all parameters meet 
#the requirements of that model
#============================
#Function used to calculate boundary position 
  boundary <- function(k, m, s, exp) {
    b = k+(m*(s^exp))
    return(b)
  }

##============= Sample Sigma Parameters ===========
#I do this before the boundary parameters for simplicity
#Upper and lower bounds defined
sigma_bounds = c(0, 40)

#Define sigma parameters
#define means
sigma_mu = c(1.23,	1.00,	0.81,	0.79) 

#define SDs 
#has to be smaller for the linear models otherwise it takes too long to sample 
if (model == 'lin') {
  sigma_sd = rep(0.2, length(sigma_mu))
} else {
  sigma_sd = rep(0.5, length(sigma_mu))
}

#sample the sigma parameters
repeat {
  # sampling sigma first
  sigmas <- (rnorm(length(sigma_mu), 0, 1)*sigma_sd)+sigma_mu
  
  #if sampling sigma - need them to be ordered consecutively
  if (all(sort(sigmas,  decreasing = TRUE) == sigmas) & #order check
      all(sigmas > sigma_bounds[1]) & #lower bound check
      all(sigmas < sigma_bounds[2])) { #upper bound check
    break
  }
}

#============= Define Boundary Parameters ============
#Once we have the sigma parameters, sample the boundary parameters
#Linear Model 
if (model == 'lin') {
  kboundaries_mu = c(0.26,	0.50,	0.92,	0.28,	-0.42,	-0.86,	-0.18)
  mboundaries_mu = c(-0.19,	-0.36, -0.53,	0.73,	2.14,	3.17,	2.78)
  exponent = 1 #fixed to 1 for linear model 
  kboundaries_sd = rep(1, 7)
  mboundaries_sd = rep(1, 7)
#Quadratic Model 
} else if (model == 'quad') {
  kboundaries_mu = c(0.17, 0.32, 0.66, 0.64, 0.62, 0.69, 1.18)
  mboundaries_mu = c(-0.09, -0.18, -0.26, 0.36, 1.06, 1.57, 1.38)
  kboundaries_sd = rep(0.5, 7)
  mboundaries_sd = rep(0.5, 7)
  exponent = 2 #fit to 2 for quadratic model 
#Free-Exponent Model 
} else if (model == 'exp') {
  kboundaries_mu = c(0.17, 0.32, 0.66, 0.64, 0.62, 0.69, 1.18)
  mboundaries_mu = c(-0.09, -0.18, -0.26, 0.36, 1.06, 1.57, 1.38)
  kboundaries_sd = rep(0.5, 7)
  mboundaries_sd = rep(0.5, 7)
  
  #sample exponent parameter
  exponent_mu = 2
  exponent_sd = 1
  exponent_bounds = c(1, 10)
  repeat {
    exponent <- rnorm(1, exponent_mu, exponent_sd)
    if ((exponent > exponent_bounds[1]) &
        (exponent < exponent_bounds[2])) { 
      break
    }
  }
}

boundaries_bounds = c(0, 70) #limits on boundary position (once calculated using k + m*sigma)

#============= Sample Boundary Parameters ============
#set a timer because sometimes this will run on an infinite loop 
#if the sampled sigma parameters suck 
start_time = Sys.time()
time_out = FALSE

repeat {    
  #sample k and m parameters
  kbounds <- (rnorm(length(kboundaries_mu), 0, 1)*kboundaries_sd)+kboundaries_mu
  mbounds <- (rnorm(length(mboundaries_mu), 0, 1)*mboundaries_sd)+mboundaries_mu
  
  #calculate boundaries using sampled k & m + chosen sigmas & exponent
  bs <- sapply(1:length(kbounds), function(x) 
    boundary(kbounds[x], mbounds[x], sigmas, exponent))
  
  #check if boundaries meet conditions 
  if (all(bs > boundaries_bounds[1]) & #lower bound check
      all(bs < boundaries_bounds[2]) & # upper bound check
      all(t(apply(bs, 1, function(x) sort(x))) == bs)) { #make sure boundaries are monotonically increasing
    break
  }
  
  #check timer 
  run_time = Sys.time()
  if (difftime(run_time, start_time, units = "mins") > 10) {
    time_out = TRUE
    break
  }
}
#============= Return Parameters ============
if (model == 'lin' | model == 'quad') {
generating_parameters <- c(kbounds, mbounds, sigmas)
} else if (model == 'exp') {
generating_parameters <- c(kbounds, mbounds, sigmas, exponent)
}

#if we couldn't find a suitable set of parameters 
if (time_out == TRUE) {
  generating_parameters = rep(NA, length(generating_parameters))
}

return(generating_parameters)
}