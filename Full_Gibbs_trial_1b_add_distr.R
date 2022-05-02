gradient <- function(x, coeff, sdy) { # sigma is labelled "sdy"
  coeff = unlist(coeff)
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  return(A + max(c(K-A,0))/((1+(exp(Q*(x-M))))) + rnorm(length(x),0,sdy))
}

loglik <- function(x, ymean, yest, sdyest, coeff, sdy) {
  # extract regression coefficients
  coeff = unlist(coeff)
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  # likelihood for yest - oxygen isotope data normal distributed
  ll1a <- sum(dnorm(yest, ymean, sdyest,log=TRUE),na.rm=T)
  # likelihood for yest - other PDFs
  # likelihood for
  #ll1b <- sum(dunif(yest[1],20,40, log = TRUE),
  #            dnorm(yest[6],6,5, log = TRUE))
  pred = A + max(c(K-A,0))/((1+(exp(Q*(x-M)))))
  ll2 <- sum(dnorm(yest, mean = pred, sd = sdy, log = TRUE))
  return(c(ll2+ll1a))
}

#logprior <- function(coeff, yest) {
# coeff = unlist(coeff)
#  return(sum(c(
#    dunif(coeff[1], -4, 40, log = TRUE),
#    dunif(coeff[2], -4, 40, log = TRUE),
#    dnorm(coeff[3], 45, 10, log = TRUE),
#    dlnorm(coeff[4], -2.2, 0.8, log = TRUE))))
#}

logprior <- function(coeff) {
  coeff = unlist(coeff)
  return(sum(c(
    dunif(coeff[1], -4, coeff[2], log = TRUE),
    dunif(coeff[2], coeff[1], 40, log = TRUE),
    dnorm(coeff[3], 45, 10, log = TRUE),
    dlnorm(coeff[4], -2.2, 0.8, log = TRUE))))
}


logposterior <- function(x, ymean, yest, sdyest, coeff, sdy){
  return (loglik(x, ymean, yest, sdyest, coeff, sdy) + logprior(coeff))
}

MH_propose_coeff <- function(coeff, prop_sd_coeff){
  return(rnorm(4,mean = coeff, sd= prop_sd_coeff))
}
MH_propose_yest <- function(yest, prop_sd_yest){
  return(mvnfast::rmvn(1,mu = yest, sigma = prop_sd_yest))
}

# function for plotting the 95 % CI shading
error_polygon <- function(x,en,ep,color) {
  polygon( c(x[1], x, x[length(x)], x[length(x)], rev(x), x[1]),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
}
#
# Main MCMCM function
run_MCMC <- function(nIter, x, yobs, yd_mu, yd_sd, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                     prop_sd_coeff, prop_sd_yest, nbin){
  ### Initialisation
  coefficients = array(dim = c(nIter,4)) # set up array to store coefficients
  coefficients[1,] = coeff_inits # initialise coefficients
  sdy = rep(NA_real_,nIter) # set up vector to store sdy
  sdy[1] = sdy_init # intialise sdy


  n_p <- length(yobs)

  n_d <- length(yd_mu)

  n_d_ind <- (n_p+1):(n_p+n_d)

  nbin <- n_p + length(yd_mu)

  yestimate = array(dim = c(nIter,nbin)) # set up array to store coefficients
  yestimate[1,] = yest_inits # initialise coefficients

  sdyest = array(dim = c(nIter,n_p)) # set up vector to store sdy
  sdyest[1,] = sdyest_inits # intialise sdy

  A_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  B_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  shape_sdy <- A_sdy+nbin/2 # shape parameter for the inverse gamma


  #### Investigate these: Need to be broad for single obser
  A_sdyest = 1 # parameter for the prior on the inverse gamma distribution of sdyest
  B_sdyest = 1 # parameter for the prior on the inverse gamma distribution of sdyest
  ####
  yn = sapply(yobs,length)
  #yn[which(is.na(yobs))] = NA
  shape_sdyest =  A_sdyest+yn/2 # shape parameter for the inverse gamma
  ### n-1?!
  ymean = c(sapply(yobs,mean),yd_mu) # mean of observations and means of given distributions
  yvar = sapply(yobs,var)
  #yvar[which(is.na(yvar))] <- max(yvar,na.rm=T)
  sumobs <- sapply(yobs,sum)

  logpost = rep(NA,nIter)


  ### The MCMC loop
  for (i in 2:nIter){
    pred = gradient(x,coefficients[i-1,],0)

    ### 1. Gibbs steps for data that have multiple points (estimate global mean and sd)

    ## 1.1.a Gibbs step to estimate yestimate
    yestimate[i,1:n_p] = rnorm(n_p,
                          #(sumobs/sdyest[i-1,]^2 + pred/sdy[i-1]) / (yn/sdyest[i-1,]^2) + (1/sdy[i-1]^2),
                         1/(1/sdy[i-1]^2 + yn/sdyest[i-1,1:n_p]^2)*(pred[1:n_p]/sdy[i-1]^2 + sumobs/sdyest[i-1,1:n_p]^2),
                         #1/sqrt(yn/sdyest[i-1,]^2 + 1/sdy[i-1]^2))
                         sqrt((1/sdy[i-1]^2 + yn/sdyest[i-1,1:n_p]^2)^(-1)))

    ## 1.1.b Gibbs step to estimate sdyest
    for(j in 1:n_p) sdyest[i,j] = sqrt(1/rgamma(1,
                                                 shape_sdyest[j],
                                                 (B_sdyest+0.5*sum((yobs[[j]]-yestimate[i,j])^2))))

    ### 1.2. Gibbs steps for data that have sample mean and sd given (estimate global mean only)
    ## 1.2.a Gibbs step to estimate yestimate
    ### distribution from here: https://people.eecs.berkeley.edu/~jordan/courses/260-spring10/lectures/lecture5.pdf
    yestimate[i,n_d_ind] = rnorm(n_d,
                                 sdy[i-1]^2/(yd_sd^2+sdy[i-1]^2)*yd_mu + yd_sd^2/(yd_sd^2+sdy[i-1]^2)*pred[n_d_ind],
                                 sqrt(1/sdy[i-1]^2 + 1/yd_sd^2))


     ### 3. Gibbs step to estimate sdy
    sdy[i] = sqrt(1/rgamma(1,
                           shape_sdy,
                          (B_sdy+0.5*sum((yestimate[i,]-pred)^2))))


    #https://stats.stackexchange.com/questions/525674/gibbs-sampler-for-normal-and-inverse-gamma-distribution-in-r
    # https://stats.stackexchange.com/questions/266665/gibbs-sampler-examples-in-r


    ## 4. Metropolis-Hastings step to estimate the regression coefficients
    proposal_coeff = MH_propose_coeff(coefficients[i-1,],prop_sd =  prop_sd_coeff) # new proposed values



    if(any(proposal_coeff[4] <= 0)) HR = 0 else {# Q needs to be >0
      # Hastings ratio of the proposal
      logpostold = logposterior(x = x, ymean = ymean, yest = yestimate[i,],
                                sdyest = c(sdyest[i,],yd_sd), coeff = coefficients[i-1,],
                                sdy = sdy[i])
      logpostnew = logposterior(x = x, ymean = ymean, yest = yestimate[i,],
                                sdyest = c(sdyest[i,],yd_sd), coeff = proposal_coeff,
                                sdy = sdy[i])
        HR = exp(logpostnew -
                 logpostold)
    }
    # accept proposal with probability = min(HR,1)
    if (runif(1) < HR){
      accept = 1
      coefficients[i,] = proposal_coeff
      #yestimate[i,] = proposal_yest
      logpost[i] = logpostnew
      # if proposal is rejected, keep the values from the previous iteration
    }else{
      accept = 0
      coefficients[i,] = coefficients[i-1,]
      #yestimate[i,] = yestimate[i-1,]
      logpost[i] = logpostold

    }
    #if(accept == 0) main = "rejected" else main = "accepted"
    #  plot(seq(0,90,0.1), gradient(seq(0,90,0.1), coefficients[i-1,],0), main = main,
    #       ylim = c(-1,40), type = "l", lwd = 3, ylab = "Temperature", xlab = "Latitude")
    #  points(seq(0,90,0.1), gradient(seq(0,90,0.1), proposal_coeff,0),
     #        type = "l", lwd = 3, lty = 3, col = "red")
#
    #  points(x,  yestimate[i,], pch = 19)
    #  points(x,  yestimate[i-1,], pch = 1)
     # points(x,ymean,col = "dodgerblue", pch = 17)
# print(i)
# print(coefficients[i,])
# print(yestimate[i,])
# print(sdyest[i,])
# print(sdy[i])
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


nbin = 4
ndbin = 3
#prop_sd_yest <- matrix(0.01,nrow = nbin, ncol = nbin)
#diag(prop_sd_yest) <- 0.5
prop_sd_coeff <- c(1,1,1,0.1)

coeff_inits = c(0,30,45,0.1)
yest_inits = rep(25,nbin+ndbin) #c(30,25,20,15,10,0,0)
sdyest_inits = rep(2,nbin)

x <- c(seq(10,60,10))
npb <- c(7,7,7,3,3,3)
ym <- c(32,30,24,17,24,10)
ysd <- c(1,1,3,2,4,2)


nbin = 4
x <- c(30,40,50,60,10,20,70)
npb <- c(5,5,5,3)
ym <- (c(25,18,15,25))
ysd <- c(1,1,1,5)
y <- lapply(1:nbin, function(x) rnorm(npb[x],ym[x],ysd[x]))
y[which(npb==0)] <- NA
plot(0,0,xlim = c(0,90), ylim = c(-1,40), type = "n")
for(i in (1:nbin)) if(npb[i] !=0)points(rep(x[i],npb[i]),y[[i]])

yd_mu <- c(33,30,10)
yd_sd <- c(2,1,2)

nIter = 50000
sdy_init = 1
yobs = y
system.time({ m2 <-  run_MCMC(nIter = nIter, x = x, yobs = y, yd_mu = yd_mu, yd_sd = yd_sd,
                 coeff_inits = coeff_inits,
                 sdy_init = sdy_init, yest_inits = yest_inits,
                 sdyest_inits = sdyest_inits,
                 prop_sd_coeff=prop_sd_coeff, prop_sd_yest=prop_sd_yest)
})
burnin = 25000+1
#nIter = 30000
mn <- m2

latitude <- 0:90
sample_it <- sample((burnin+1):nIter,2000)
grad_025 <- sapply(1:length(latitude), function(f) quantile(apply(mn[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.025))
grad_975 <- sapply(1:length(latitude), function(f) quantile(apply(mn[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.975))

plot(seq(0,90,0.1), gradient(seq(0,90,0.1), apply(mn[[1]][burnin:nIter,1:4],2,median), 0),
     ylim = c(-1,40), type = "l", lwd = 3, ylab = "Temperature", xlab = "Latitude")
### Calculate confidence intervals
error_polygon(latitude,grad_025,grad_975,rgb(0,0,0,0.15))

for(i in 1:5) points(latitude,gradient(latitude,unlist(mn[[1]][sample(burnin:nIter,1),1:4]),0), type = "l", lwd = 2, col = rgb(0,0,0,0.3))

#for(i in 1:nbin) if(!(is.na(prior_df[i,1])) & prior_df[i,1]=="dunif")  lines(x[c(i,i)],c(prior_df[i,2],prior_df[i,3]), lwd = 7, col = rgb(0,0.9,0,0.5))
#for(i in 1:nbin) {if(!(is.na(prior_df[i,1])) & prior_df[i,1]=="dnorm")  {beanplot::beanplot(rnorm(1000,prior_df[i,2],prior_df[i,3]), add = T,
#                                                                                            at = x[i], maxwidth = 5, side = "second",
#                                                                what = c(0,1,0,0), col = rgb(0,0.7,0.7,0.5), border = NA)}}

for(d in 1:3) beanplot::beanplot(rnorm(2000,yd_mu[d],yd_sd[d]), add = T,
                   at = x[c(5:7)[d]], maxwidth = 5, side = "second",
                   what = c(0,1,0,0), col = rgb(0,0.7,0.7,0.5), border = NA)


offset <- rnorm(7,0,0.05)
points(x[1:7]+offset,apply(mn[[2]][burnin:nIter,1:7],2,mean), col = "red", pch = 19, cex = 1.25)
sapply(1:7, function(a) points(c(x[a]+offset[a],x[a]+offset[a]),quantile(mn[[2]][burnin:nIter,a], probs = c(0.25,0.75)), type = "l", col= rgb(1,0,0,0.5), lwd =5))
sapply(1:7, function(a) points(c(x[a]+offset[a],x[a]+offset[a]),quantile(mn[[2]][burnin:nIter,a], probs = c(0.125,0.875)), type = "l", col= rgb(1,0,0,0.5), lwd =3))
sapply(1:7, function(a) points(c(x[a]+offset[a],x[a]+offset[a]),quantile(mn[[2]][burnin:nIter,a], probs = c(0.025,0.975)), type = "l", col= rgb(1,0,0,0.5), lwd =1))

#sapply(1:nbin, function(a) points(c(x[a],x[a]),c(median(mn[[2]][burnin:nIter,a])+c(-2,2)*median(mn[[3]][burnin:nIter,a])), type = "l", col= "red"))

for(i in (1:nbin)) if(npb[i] !=0) points(rep(x[i],npb[i])+rnorm(npb[i],0,0.25),y[[i]], col = rgb(0,0.3,1,0.55), pch = 17)

points(x[1:4],sapply(y,mean), col = rgb(0,.85,.425,.75), pch = 4, lwd = 2)





legend("topright",legend=c("regression line", "temperature estimate","dO18 data", "coral reef range (uniform prior)", "additional info. (normal prior)"),
       cex = 0.8, col = c("black", "red", rgb(0,0.3,1,0.55),  rgb(0,0.9,0,0.5), rgb(0,0.7,0.7,0.5)), lwd = c(2,NA,NA,4,4),
       pch = c(NA,19,17,NA,NA), pt.cex = c(NA,1.25,1,NA,NA))



plot(mn[[1]][,4])#,mn[[2]][,6])

loglik(x,sapply(y,mean), mn[[2]][nIter,1:nbin], mn[[3]][nIter,1:nbin], mn[[1]][nIter,1:4], mn[[1]][nIter,5])
logprior(mn[[1]][nIter,1:4])


# Write Jags model for comparison

dunif(coeff[1], -4, 40, log = TRUE),
dunif(coeff[2], -4, 40, log = TRUE),
dnorm(coeff[3], 45, 10, log = TRUE),
dlnorm(coeff[4], -2, .1, log = TRUE)
