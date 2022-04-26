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
    sum(sapply(1:nbin, function(b) dunif(yest[b], -4, 40, log = TRUE))))))
}

logposterior <- function(x, yobs, yest, sdyest, coeff, sdy){
  return (loglik(x, yobs, yest, sdyest, coeff, sdy) + logprior(coeff,yest))
}

MH_propose_coeff <- function(coeff, prop_sd_coeff){
  return(rnorm(4,mean = coeff, sd= prop_sd_coeff))
}
MH_propose_yest <- function(yest, prop_sd_yest){
  return(mvnfast::rmvn(1,mu = yest, sigma = prop_sd_yest))
}


# Main MCMCM function
run_MCMC <- function(nIter, x, yobs, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                     prop_sd_coeff, prop_sd_yest, nbin){
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
  sdyest[1,] = sdyest_inits # intialise sdy
  A_sdyest = 3 # parameter for the prior on the inverse gamma distribution of sdy
  B_sdyest = 0.1 # parameter for the prior on the inverse gamma distribution of sdy
  nest = sapply(yobs,length)
  shape_sdyest =  A_sdyest+nest/2 # shape parameter for the inverse gamma

  ymean = sapply(yobs,mean)
  yvar = sapply(yobs,var)

  ### The MCMC loop
  for (i in 2:nIter){

    ## 1. Gibbs step to estimate sdy
    sdy[i] = sqrt(1/rgamma(
      1,shape_sdy,B_sdy+0.5*sum((yestimate[i-1]-gradient(x,coefficients[i-1,],0))^2)))

    ## 2. Gibbs step to estimate sdyest
    for(j in 1:nbin) sdyest[i,j] = sqrt(1/rgamma(
      1,shape_sdyest[j],B_sdyest+0.5*sum((yobs[[j]]-yestimate[i-1,j])^2)))

    #https://stats.stackexchange.com/questions/525674/gibbs-sampler-for-normal-and-inverse-gamma-distribution-in-r
    # https://stats.stackexchange.com/questions/266665/gibbs-sampler-examples-in-r

    ## 3. Metropolis-Hastings step to estimate yest
    proposal_yest = MH_propose_yest(yestimate[i-1,],prop_sd_yest =  prop_sd_yest) # new proposed values

    ## 4. Metropolis-Hastings step to estimate the regression coefficients
    proposal_coeff = MH_propose_coeff(coefficients[i-1,],prop_sd =  prop_sd_coeff) # new proposed values
    if(any(proposal_coeff[4] <= 0)) HR = 0 else # Q needs to be >0
      # Hastings ratio of the proposal
      HR = exp(logposterior(x = x, yobs = yobs, yest = proposal_yest, sdyest = sdyest[i,], coeff = proposal_coeff, sdy = sdy[i]) -
                 logposterior(x = x, yobs = yobs, yest = yestimate[i-1,],  sdyest = sdyest[i,], coeff = coefficients[i-1,], sdy = sdy[i]))
    # accept proposal with probability = min(HR,1)
    if (runif(1) < HR){
      coefficients[i,] = proposal_coeff
      yestimate[i,] = proposal_yest

      # if proposal is rejected, keep the values from the previous iteration
    }else{
      coefficients[i,] = coefficients[i-1,]
      yestimate[i,] = yestimate[i-1,]

    }
  } # end of the MCMC loop

  ###  Function output
  output = list(data.frame(A = coefficients[,1],
                      K = coefficients[,2],
                      M = coefficients[,3],
                      Q = coefficients[,4],
                      sdy = sdy),
                yestimate = yestimate,
                sdyest = sdyest)
  return(output)
}


nbin = 7
prop_sd_yest <- matrix(0.99,nrow = nbin, ncol = nbin)
diag(prop_sd_yest) <- 1
prop_sd_coeff <- c(1,1,1,0.01)

coeff_inits = c(10,30,45,0.2)
yest_inits = c(30,25,20,15,10,0,0)
sdyest_inits = rep(2,7)

x <- seq(10,70,10)
y <- lapply(c(30,27,20,13,8,5,3), function(x) rnorm(7,x,2))
plot(0,0,xlim = c(10,70), ylim = c(-1,35))
for(i in 1:nbin) points(rep(x[i],7),y[[i]])

m <- run_MCMC(nIter = 50000, x, y, coeff_inits = coeff_inits,
                     sdy_init = 1, yest_inits = yest_inits,
                     sdyest_inits = sdyest_inits,
                     prop_sd_coeff, prop_sd_yest,
                      nbin = 7)

ind <- seq(1,50000,10)
plot(m[[1]][ind,2])

plot(m[[1]][ind,1],m[[1]][ind,2])

plot(m[[2]][ind,1], type = "l")

plot(m[[3]][ind,6])

ind <- 5
hist(rnorm(10000,median(m[[2]][,ind]),median(m[[3]][,ind])),100)
abline(v = y[[ind]], lwd = 2, col = "turquoise")


plot(seq(10,70,0.1), gradient(seq(10,70,0.1), apply(m[[1]][20000:50000,1:4],2,median), 0), ylim = c(0,30))
points(x,apply(m[[2]],2,median), col = "red", pch = 19)

x1 <- seq(-1,4,0.01)
plot(x1,dnorm(x1,0,1))#*dnorm(x1,1,1))
points(x1,dnorm(x1,2,1), col = "dodgerblue")
points(x1,dnorm(x1,0,1)*dnorm(x1,2,1)*4, col = "red")
