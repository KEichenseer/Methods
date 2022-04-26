gradient <- function(x, coeff, sdy) { # sigma is labelled "sdy"
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  return(A + max(c(K-A,0))/((1+(exp(Q*(x-M))))) + rnorm(length(x),0,sdy))
}

loglik <- function(x, ymean, yest, sdyest, coeff, sdy) {
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]

  ll1 <- sum(dnorm(ymean, yest, sdyest,log=TRUE))
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

logposterior <- function(x, ymean, yest, sdyest, coeff, sdy){
  return (loglik(x, ymean, yest, sdyest, coeff, sdy) + logprior(coeff,yest))
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
  B_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  n <- length(nbin)
  shape_sdy <- A_sdy+n/2 # shape parameter for the inverse gamma

  yestimate = array(dim = c(nIter,nbin)) # set up array to store coefficients
  yestimate[1,] = yest_inits # initialise coefficients
  sdyest = array(dim = c(nIter,nbin)) # set up vector to store sdy
  sdyest[1,] = sdyest_inits # intialise sdy
  A_sdyest = 3 # parameter for the prior on the inverse gamma distribution of sdy
  B_sdyest = 1 # parameter for the prior on the inverse gamma distribution of sdy
  nest = sapply(yobs,length)
  shape_sdyest =  A_sdyest+nest/2 # shape parameter for the inverse gamma

  ymean = sapply(yobs,mean)
  yvar = sapply(yobs,var)

  logpost = rep(NA,nIter)

  ### The MCMC loop
  for (i in 2:nIter){

    ## 1. Gibbs step to estimate sdy
    sdy[i] = sqrt(1/rgamma(
      1,shape_sdy,B_sdy+0.5*sum((yestimate[i-1,]-gradient(x,coefficients[i-1,],0))^2)))

    ## 2. Gibbs step to estimate sdyest
   # for(j in 1:nbin) sdyest[i,j] = sqrt(1/rgamma(
   #    1,shape_sdyest[j],B_sdyest+0.5*sum((yobs[[j]]-yestimate[i-1,j])^2)))
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
      logpostold = logposterior(x = x, ymean = ymean, yest = yestimate[i-1,],  sdyest = sdyest[i,], coeff = coefficients[i-1,], sdy = sdy[i])
      logpostnew = logposterior(x = x, ymean = ymean, yest = proposal_yest, sdyest = sdyest[i,], coeff = proposal_coeff, sdy = sdy[i])
      HR = exp(logpostnew -
                 logpostold)
    # accept proposal with probability = min(HR,1)
    if (runif(1) < HR){
      coefficients[i,] = proposal_coeff
      yestimate[i,] = proposal_yest
      logpost[i] = logpostnew
      # if proposal is rejected, keep the values from the previous iteration
    }else{
      coefficients[i,] = coefficients[i-1,]
      yestimate[i,] = yestimate[i-1,]
      logpost[i] = logpostold

    }
  } # end of the MCMC loop

  ###  Function output
  output = list(data.frame(A = coefficients[,1],
                      K = coefficients[,2],
                      M = coefficients[,3],
                      Q = coefficients[,4],
                      sdy = sdy,
                      logpost = logpost),
                yestimate = yestimate,
                sdyest = sdyest)
  return(output)
}


nbin = 7
npb = 3
prop_sd_yest <- matrix(0.9,nrow = nbin, ncol = nbin)
diag(prop_sd_yest) <- 1
prop_sd_coeff <- c(.5,.5,.5,0.01)

coeff_inits = c(0,30,45,0.1)
yest_inits = c(30,29,28,10,9,8,7) #c(30,25,20,15,10,0,0)
sdyest_inits = rep(2,7)

x <- seq(10,70,10)
npb <- c(3,3,3,3,3,3,12)
ym <- c(30,27,20,12,10,14,4)
y <- lapply(1:nbin, function(x) rnorm(npb[x],ym[x],2))
plot(0,0,xlim = c(10,70), ylim = c(-1,35))
for(i in 1:nbin) points(rep(x[i],npb[i]),y[[i]])

m2 <- run_MCMC(nIter = 50000, x = x, yobs = y, coeff_inits = coeff_inits,
                     sdy_init = 1, yest_inits = yest_inits,
                     sdyest_inits = sdyest_inits,
                     prop_sd_coeff, prop_sd_yest,
                      nbin = 7)
mn <- m1
ind <- seq(1,2000,1)
plot(ind,mn[[1]][ind,5], type = "l")
plot(ind,mn[[2]][ind,3], type = "l")

plot(ind,m[[3]][ind,7], type = "l")

plot(m[[1]][ind,1],m[[1]][ind,2])
abline(a=0,b=1)
plot(m[[1]][ind,2], type = "l")

plot(m[[3]][ind,6])

ind <- 5
hist(rnorm(10000,median(m[[2]][,ind]),median(m[[3]][,ind])),100)
abline(v = y[[ind]], lwd = 2, col = "turquoise")

burnin = 20000+1
nIter = 100000
plot(seq(10,70,0.1), gradient(seq(10,70,0.1), apply(m[[1]][burnin:nIter,1:4],2,median), 0), ylim = c(-1,35))
points(x,apply(m[[2]][burnin:nIter,],2,median), col = "red", pch = 19)
sapply(1:nbin, function(a) points(c(x[a],x[a]),c(median(m[[2]][burnin:nIter,a])+c(-2,2)*median(m[[3]][burnin:nIter,a])), type = "l", col= "red"))
for(i in 1:nbin) points(rep(x[i],7),y[[i]], col = rgb(0,0.8,0.6,0.7))

points(x,apply(m[[2]][burnin:nIter,],2,median), col = rgb(0.5,1,0,0.5), pch = 19)
points(x,apply(m2[[2]][burnin:nIter,],2,median), col = rgb(0.5,0,1,0.5), pch = 19)

burnin = 20000+1
nIter = 100000
plot(seq(10,70,0.1), gradient(seq(10,70,0.1), apply(m2[[1]][burnin:nIter,1:4],2,median), 0), ylim = c(-1,35))
points(x,apply(m2[[2]][burnin:nIter,],2,median), col = "red", pch = 19)
sapply(1:nbin, function(a) points(c(x[a],x[a]),c(median(m2[[2]][burnin:nIter,a])+c(-2,2)*median(m2[[3]][burnin:nIter,a])), type = "l", col= "red"))
for(i in 1:nbin) points(rep(x[i],7),y[[i]], col = rgb(0,0.8,0.6,0.7))


plot(seq(10,70,0.1), gradient(seq(10,70,0.1), apply(m3[[1]][burnin:nIter,1:4],2,median), 0), ylim = c(-1,35))
points(x,apply(m3[[2]][burnin:nIter,],2,median), col = "red", pch = 19)
sapply(1:nbin, function(a) points(c(x[a],x[a]),c(median(m3[[2]][burnin:nIter,a])+c(-2,2)*median(m3[[3]][burnin:nIter,a])), type = "l", col= "red"))
for(i in 1:nbin) points(rep(x[i],7),y[[i]], col = rgb(0,0.8,0.6,0.7))

burnin = 20000+1
nIter = 50000
mn <- m2
plot(seq(10,70,0.1), gradient(seq(10,70,0.1), apply(mn[[1]][burnin:nIter,1:4],2,median), 0), ylim = c(-1,35))
points(x,apply(mn[[2]][burnin:nIter,],2,median), col = "orange", pch = 19)
sapply(1:nbin, function(a) points(c(x[a],x[a]),quantile(mn[[2]][burnin:nIter,a], probs = c(0.025,0.975)), type = "l", col= rgb(.8,.8,0,0.5), lwd =3))
sapply(1:nbin, function(a) points(c(x[a],x[a]),c(median(mn[[2]][burnin:nIter,a])+c(-2,2)*median(mn[[3]][burnin:nIter,a])), type = "l", col= "red"))

for(i in 1:nbin) points(rep(x[i],npb[i]),y[[i]], col = rgb(0,0.8,0.6,0.7))


x1 <- seq(-1,4,0.01)
plot(x1,dnorm(x1,0,1))#*dnorm(x1,1,1))
points(x1,dnorm(x1,2,1), col = "dodgerblue")
points(x1,dnorm(x1,0,1)*dnorm(x1,2,1)*4, col = "red")
