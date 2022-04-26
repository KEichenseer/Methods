gradient <- function(x, coeff, sdy) { # sigma is labelled "sdy"
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  return(A + max(c(K-A,0))/((1+(exp(Q*(x-M))))) + rnorm(length(x),0,sdy))
}

loglik <- function(x, yobs, yest, sdyest, coeff, sdy) {
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  
  ll1 <- sum(sapply(1:nbin, function(b) sum(dnorm(yobs[[b]], yest[b], sdyest[b],log=TRUE))))
  pred = A + max(c(K-A,0))/((1+(exp(Q*(x-M)))))
  ll2 <- sum(dnorm(yest, mean = pred, sd = sdy, log = TRUE))
  return(ll1+ll2)
}

logprior <- function(coeff, yest) {
  return(sum(c(
    dunif(coeff[1], -4, 40, log = TRUE),
    dunif(coeff[2], -4, 40, log = TRUE),
    dnorm(coeff[3], 45, 10, log = TRUE),
    dlnorm(coeff[4], -2, 1, log = TRUE),
    sum(sapply(1:nbin, function(b) dunif(yest[[b]], -4, 40, log = TRUE))))))
}

logposterior <- function(x, y, coeff, sdy){
  return (loglik(x, y, coeff, sdy) + logprior(coeff))
}

MH_propose_coeff <- function(coeff, prop_sd_coeff){
  return(rnorm(4,mean = coeff, sd= prop_sd_coeff))
}
MH_propose_yest <- function(coeff, prop_sd_yest){
  return(rnorm(nbin,mean = yest, sd = prop_sd_yest))
}


# Main MCMCM function
run_MCMC <- function(nIter, x, y, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                     prop_sd_coeff, prop_sd_yest){
  ### Initialisation
  coefficients = array(dim = c(nIter,4)) # set up array to store coefficients
  coefficients[1,] = coeff_inits # initialise coefficients
  sdy = rep(NA_real_,nIter) # set up vector to store sdy
  sdy[1] = sdy_init # intialise sdy
  A_sdy = 3 # parameter for the prior on the inverse gamma distribution of sdy
  B_sdy = 0.1 # parameter for the prior on the inverse gamma distribution of sdy
  n <- length(nbin)
  shape_sdy <- A_sdy+n/2 # shape parameter for the inverse gamma
  
  yestimate = array(dim = c(nIter,nbin)) # set up array to store coefficients
  yestimate[1,] = yest_inits # initialise coefficients
  sdyest = array(dim = c(nIter,nbin)) # set up vector to store sdy
  sdyest[1] = sdyest_inits # intialise sdy
  A_sdyest = 3 # parameter for the prior on the inverse gamma distribution of sdy
  B_sdyest = 0.1 # parameter for the prior on the inverse gamma distribution of sdy
  nest <- sapply(y,length)
  shape_sdyest <- A_sdy+nest/2 # shape parameter for the inverse gamma
  
  ### The MCMC loop
  for (i in 2:nIter){ 
    
    ## 1. Gibbs step to estimate sdy
    sdy[i] = sqrt(1/rgamma(
      1,shape_sdy,B_sdy+0.5*sum((y-gradient(x,coefficients[i-1,],0))^2)))
    
    ## 2. Gibbs step to estimate yest and sdyest
    mu[i]  <- rnorm(n = 1, mean = ybar, sd = sqrt(1 / (n * tau[i - 1])))    
    tau[i] <- rgamma(n = 1, shape = n / 2, scale = 2 / ((n - 1) * s2 + n * (mu[i] - ybar)^2))
    # https://stats.stackexchange.com/questions/266665/gibbs-sampler-examples-in-r
    sdyest[i] = sqrt(1/rgamma(
      1,shape_sdy,B_sdy+0.5*sum((y-gradient(x,coefficients[i-1,],0))^2)))
    
    ## 2. Metropolis-Hastings step to estimate the regression coefficients
    proposal = MH_propose(coefficients[i-1,],proposal_sd =  c(.5,.5,.5,0.01)) # new proposed values
    if(any(proposal[4] <= 0)) HR = 0 else # Q needs to be >0
      # Hastings ratio of the proposal
      HR = exp(logposterior(x = x, y = y, coeff = proposal, sdy = sdy[i]) -
                 logposterior(x = x, y = y, coeff = coefficients[i-1,], sdy = sdy[i]))
    # accept proposal with probability = min(HR,1)
    if (runif(1) < HR){ 
      coefficients[i,] = proposal
      # if proposal is rejected, keep the values from the previous iteration
    }else{
      coefficients[i,] = coefficients[i-1,]
    }
  } # end of the MCMC loop
  
  ###  Function output
  output = data.frame(A = coefficients[,1],
                      K = coefficients[,2],
                      M = coefficients[,3],
                      Q = coefficients[,4],
                      sdy = sdy)
  return(output)
}





x1 <- seq(-1,4,0.01)
plot(x1,dnorm(x1,0,1))#*dnorm(x1,1,1))
points(x1,dnorm(x1,2,1), col = "dodgerblue")
points(x1,dnorm(x1,0,1)*dnorm(x1,2,1)*4, col = "red")
