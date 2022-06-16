# gradient <- function(x, coeff, sdy) { # sigma is labelled "sdy"
#   coeff = unlist(coeff)
#   A = coeff[1]
#   K = coeff[2]
#   M = coeff[3]
#   Q = coeff[4]
#   return(A + max(c(K-A,0))/((1+(exp(Q*(x-M))))) + rnorm(length(x),0,sdy))
# }
#
#
# loglik_norm <- function(x, yest, ymean, sdyest, coeff, sdy) {
#   # extract regression coefficients
#   coeff = unlist(coeff)
#   A = coeff[1]
#   K = coeff[2]
#   M = coeff[3]
#   Q = coeff[4]
#   # likelihood for yest - oxygen isotope data normal distributed
#   ll1 <- sum(dnorm(yest, ymean, sdyest,log=TRUE),na.rm=T)
#
#   # likelihood for yest - other PDFs
#   # likelihood for
#   #ll1b <- sum(dunif(yest[1],20,40, log = TRUE),
#   #            dnorm(yest[6],6,5, log = TRUE))
#   pred = A + max(c(K-A,0))/((1+(exp(Q*(x-M)))))
#   ll2 <- sum(dnorm(yest, mean = pred, sd = sdy, log = TRUE))
#   return(c(ll2+ll1))
# }
#
#
# loglik_skew <- function(x, yest, mu, sigma, lambda, coeff, sdy) {
#   # extract regression coefficients
#   coeff = unlist(coeff)
#   A = coeff[1]
#   K = coeff[2]
#   M = coeff[3]
#   Q = coeff[4]
#   # likelihood for yest - oxygen isotope data normal distributed
#   ll1 <- sum(log(2/sigma)+dnorm((yest-mu)/sigma,log=T)+pnorm(lambda*(yest-mu)/sigma,log=T))
#
#
#   # likelihood for yest - other PDFs
#   # likelihood for
#   #ll1b <- sum(dunif(yest[1],20,40, log = TRUE),
#   #            dnorm(yest[6],6,5, log = TRUE))
#   pred = A + max(c(K-A,0))/((1+(exp(Q*(x-M)))))
#   ll2 <- sum(dnorm(yest, mean = pred, sd = sdy, log = TRUE))
#   return(c(ll2+ll1))
# }
#
#
# logprior <- function(coeff) {
#   coeff = unlist(coeff)
#   return(sum(c(
#     dunif(coeff[1], -4, coeff[2], log = TRUE),
#     dunif(coeff[2], coeff[1], 40, log = TRUE),
#     dnorm(coeff[3], 45, 10, log = TRUE),
#     dlnorm(coeff[4], -2.2, 0.8, log = TRUE))))
# }


##############
##############  Parametrise with difference between hot and cold end instead

gradient <- function(x, coeff, sdy) { # parametrise with difference between cold and hot end instead
  coeff = unlist(coeff)
  A = coeff[1]
  DKA = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  return(A + DKA/((1+(exp(Q*(x-M))))) + rnorm(length(x),0,sdy))
}

loglik_norm <- function(x, yest, ymean, sdyest, coeff, sdy) {
  # extract regression coefficients
  coeff = unlist(coeff)
  A = coeff[1]
  DKA = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  # likelihood for yest - oxygen isotope data normal distributed
  ll1 <- sum(dnorm(yest, ymean, sdyest,log=TRUE),na.rm=T)

  # likelihood for yest - other PDFs
  # likelihood for
  #ll1b <- sum(dunif(yest[1],20,40, log = TRUE),
  #            dnorm(yest[6],6,5, log = TRUE))
  pred = A + DKA/((1+(exp(Q*(x-M)))))
  ll2 <- sum(dnorm(yest, mean = pred, sd = sdy, log = TRUE))
  return(c(ll2+ll1))
}


loglik_skew <- function(x, yest, mu, sigma, lambda, coeff, sdy) {
  # extract regression coefficients
  coeff = unlist(coeff)
  A = coeff[1]
  DKA = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  # likelihood for yest - oxygen isotope data normal distributed
  ll1 <- sum(log(2/sigma)+dnorm((yest-mu)/sigma,log=T)+pnorm(lambda*(yest-mu)/sigma,log=T))


  # likelihood for yest - other PDFs
  # likelihood for
  #ll1b <- sum(dunif(yest[1],20,40, log = TRUE),
  #            dnorm(yest[6],6,5, log = TRUE))
  pred = A + DKA/((1+(exp(Q*(x-M)))))
  ll2 <- sum(dnorm(yest, mean = pred, sd = sdy, log = TRUE))
  return(c(ll2+ll1))
}


logprior <- function(coeff) {
  coeff = unlist(coeff)
  return(sum(c(
    log(2/12)+dnorm((coeff[1] - (-2.85))/12,log=T)+pnorm(10*(coeff[1] - (-2.85))/12,log=T),  #dunif(coeff[1], -4, 40, log = TRUE),
    dtnorm(coeff[2], 0, Inf,25,12, log = TRUE),
    dnorm(coeff[3], 45, 10, log = TRUE),
    dlnorm(coeff[4], -2.2, 0.8, log = TRUE))))
}

# function to generate truncated normal
dtnorm <- function(x,lower,upper,mean,sd, log = FALSE) {
  ret <- numeric(length(x))
  ret[x < lower | x > upper] <- if (log)
    -Inf
  else 0
  ret[upper < lower] <- NaN
  ind <- x >= lower & x <= upper
  if (any(ind)) {
    denom <- pnorm(upper, mean, sd) - pnorm(lower, mean,
                                            sd)
    xtmp <- dnorm(x, mean, sd, log)
    if (log)
      xtmp <- xtmp - log(denom)
    else xtmp <- xtmp/denom
    ret[x >= lower & x <= upper] <- xtmp[ind]
  }
  ret
} # from msm



logposterior_norm <- function(x, yest, ymean, sdyest, coeff, sdy){
  return (loglik_norm(x, yest, ymean, sdyest, coeff, sdy))
}

logposterior_skew <- function(x, yest, mu, sigma, lambda, coeff, sdy){
  return (loglik_skew(x, yest, mu, sigma, lambda, coeff, sdy))
}
#
# logpostold_lik_wrapper <- function(n_p,n_norm,n_skew) {
#   out <- 0
#   if(n_p != 0) out <- out + logposterior_norm(x = x[n_p_ind], yest = yestimate[i,n_p_ind], ymean = yobs_mean,
#                                               sdyest = sdyest[i,], coeff = coefficients[i-1,],
#                                               sdy = sdy[i])
#   if(n_norm != 0) out <- out + logposterior_norm(x = x[n_norm_ind], yest = yestimate[i,n_norm_ind], ymean = ynorm_mu,
#                                                  sdyest = ynorm_sd, coeff = coefficients[i-1,],
#                                                  sdy = sdy[i])
#   if(n_skew != 0) out <- out +  logposterior_skew(x = x[n_skew_ind], yest = yestimate[i,n_skew_ind], mu = yskew_mu, yskew_sigma, yskew_lambda,
#                                                   coeff = coefficients[i-1,], sdy[i])
#   return(out)
# }
#
# logpostnew_lik_wrapper <- function(n_p,n_norm,n_skew) {
#   out <- 0
#   if(n_p != 0) out <- out + logposterior_norm(x = x[n_p_ind], yest = yestimate[i,n_p_ind], ymean = yobs_mean,
#                                               sdyest = sdyest[i,], coeff = proposal_coeff,
#                                               sdy = sdy[i])
#   if(n_norm != 0) out <- out + logposterior_norm(x = x[n_norm_ind], yest = yestimate[i,n_norm_ind], ymean = ynorm_mu,
#                                                  sdyest = ynorm_sd, coeff = proposal_coeff,
#                                                  sdy = sdy[i])
#   if(n_skew != 0) out <- out +  logposterior_skew(x = x[n_skew_ind], yest = yestimate[i,n_skew_ind], mu = yskew_mu, yskew_sigma, yskew_lambda,
#                                                   coeff = proposal_coeff, sdy[i])
#   return(out)
#
# }


MH_propose_coeff <- function(coeff, prop_sd_coeff){
  return(rnorm(4,mean = coeff, sd= prop_sd_coeff))
}



#MH_propose_yest <- function(yest, prop_sd_yest){
#  return(mvnfast::rmvn(1,mu = yest, sigma = prop_sd_yest))
#}

# Gibbs sampling of mu with skew normal likelihood and normal prior
### this is for single observations, so no mean(x) or mean(y)
skew_mu <- function(x, y, sigma, rho, mu_prior, sigma_prior) {
  n1 = 1
  rnorm(length(x),(n1*sigma_prior^2*(x-sigma*rho*y)+sigma^2*(1-rho^2)*mu_prior)/
          (n1*sigma_prior^2+sigma^2*(1-rho^2)),
        sqrt((sigma_prior^2*sigma^2*(1-rho^2))/(n1*sigma_prior^2+sigma^2*(1-rho^2))) )
  #rnorm(length(x),(n1*sigma_prior^2*((x)-sigma*rho*(y))+sigma^2*(1-rho^2)*mu_prior)/
  #        (n1*sigma_prior^2+sigma^2*(1-rho^2)),
  #      sqrt((sigma_prior^2*sigma^2*(1-rho^2))/(n1*sigma_prior^2+sigma^2*(1-rho^2))) )
}





# function for plotting the 95 % CI shading
error_polygon <- function(x,en,ep,color) {
  polygon( c(x[1], x, x[length(x)], x[length(x)], rev(x), x[1]),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
}
#
# Main MCMCM function
run_MCMC <- function(nIter = 1000, obsmat = NULL, distrmat = NULL, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                     prop_sd_coeff, quiet = FALSE){
  ### Initialisation
  coefficients = array(dim = c(nIter,4)) # set up array to store coefficients
  coefficients[1,] = coeff_inits # initialise coefficients
  sdy = rep(NA_real_,nIter) # set up vector to store sdy
  sdy[1] = sdy_init # intialise sdy



  x = NULL
  if(!is.null(obsmat)) {
  ylist <- lapply(unique(obsmat$sample), function(x) obsmat$temperature[which(obsmat$sample == x)])
  x <- as.numeric(sapply(unique(obsmat$sample), function(x) unique(obsmat$latitude[which(obsmat$sample == x)])))
  n_p <- length(ylist)
  n_p_ind <- 1:n_p
  } else n_p <- 0

  if(!is.null(distrmat)) {

    x <- c(x,as.numeric(distrmat$latitude))

    ynorm_mu <- as.numeric(distrmat$location[which(distrmat$distribution == "normal")])
    ynorm_sd <- as.numeric(distrmat$scale[which(distrmat$distribution == "normal")])

    yskew_mu <- as.numeric(distrmat$location[which(distrmat$distribution == "skew-normal")])
    yskew_sigma <- as.numeric(distrmat$scale[which(distrmat$distribution == "skew-normal")])
    yskew_lambda <- as.numeric(distrmat$shape[which(distrmat$distribution == "skew-normal")])

    n_norm <- length(ynorm_mu)
    n_norm_ind <- (n_p+1):(n_p+n_norm)

    n_skew <- length(yskew_mu)

    n_skew_ind <- (n_p+n_norm+1):(n_p+n_norm+n_skew)


  } else {
  n_norm <- 0
  n_skew <- 0
  }

  nbin <- n_p + n_norm + n_skew

  yestimate = array(dim = c(nIter,nbin)) # set up array to store coefficients
  yestimate[1,] = yest_inits # initialise coefficients

  A_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  B_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  shape_sdy <- A_sdy+nbin/2 # shape parameter for the inverse gamma


  sdyest = array(dim = c(nIter,n_p)) # set up vector to store sdy
  sdyest[1,] = sdyest_inits # intialise sdy


  ymean = NULL
  if(n_p != 0) {
    #### Investigate these: Need to be broad for single obser
    A_sdyest = 1 # parameter for the prior on the inverse gamma distribution of sdyest
    B_sdyest = 1 # parameter for the prior on the inverse gamma distribution of sdyest
    ####
    yn = sapply(ylist,length)
    #yn[which(is.na(yobs))] = NA
    shape_sdyest =  A_sdyest+yn/2 # shape parameter for the inverse gamma
    ### n-1?!
    yobs_mean= c(sapply(ylist,mean))
    yobs_var = sapply(ylist,var)
    sumobs <- sapply(ylist,sum)

  }

    if(n_skew != 0) yskew_rho <- -yskew_lambda/sqrt(1+yskew_lambda^2) # not sure why but this needs to be negative

  logpost = rep(NA,nIter)
  # start progress bar
  #if (!quiet) cli::cli_progress_bar('Sampling', total = nIter)
  ### The MCMC loop
  for (i in 2:nIter){
    # update progress bar
    #if (!quiet) cli::cli_progress_update(set = i, status = paste0('iteration ', i))

    pred = gradient(x,coefficients[i-1,],0)

    ### 1. Gibbs steps for data that have multiple points (estimate global mean and sd)

    ## 1.1.a Gibbs step to estimate yestimate
    if(n_p != 0) {
    yestimate[i,1:n_p] = rnorm(n_p,
                               #(sumobs/sdyest[i-1,]^2 + pred/sdy[i-1]) / (yn/sdyest[i-1,]^2) + (1/sdy[i-1]^2),
                               1/(1/sdy[i-1]^2 + yn/sdyest[i-1,1:n_p]^2)*(pred[1:n_p]/sdy[i-1]^2 + sumobs/sdyest[i-1,1:n_p]^2),
                               #1/sqrt(yn/sdyest[i-1,]^2 + 1/sdy[i-1]^2))
                               sqrt((1/sdy[i-1]^2 + yn/sdyest[i-1,1:n_p]^2)^(-1)))

    ## 1.1.b Gibbs step to estimate sdyest
    for(j in 1:n_p) sdyest[i,j] = sqrt(1/rgamma(1,
                                                shape_sdyest[j],
                                                (B_sdyest+0.5*sum((ylist[[j]]-yestimate[i,j])^2))))
    }

    ### 1.2. Gibbs steps for data that have sample mean and sd given (estimate global mean only)
    ## 1.2.a Gibbs step to estimate yestimate
    ### distribution from here: https://people.eecs.berkeley.edu/~jordan/courses/260-spring10/lectures/lecture5.pdf
    if(n_norm != 0) {
    yestimate[i,n_norm_ind] = rnorm(n_norm,
                                 sdy[i-1]^2/(ynorm_sd^2+sdy[i-1]^2)*ynorm_mu + ynorm_sd^2/(ynorm_sd^2+sdy[i-1]^2)*pred[n_norm_ind],
                                 sqrt(1/sdy[i-1]^2 + 1/ynorm_sd^2))
    }
    ### 1.3. Gibbs steps for data that have location, scale and shape parameter given (skew-normal). Estimate global mean only)
    ## 1.3.a Gibbs step to estimate yestimate
    ### distribution from here: http://koreascience.or.kr/article/JAKO200504840590864.pdf
    if(n_skew != 0) {
    z <- (yskew_mu - yestimate[i-1,n_skew_ind])/yskew_sigma
    y <- truncnorm::rtruncnorm(n_skew,0,Inf,yskew_rho*z,sqrt(1-yskew_rho^2))
    yestimate[i,n_skew_ind] = skew_mu(x=yskew_mu, y=y, sigma=yskew_sigma, rho=yskew_rho,
                                      mu_prior=pred[n_skew_ind], sigma_prior=sdy[i-1])
    }



    ### 3. Gibbs step to estimate sdy
    sdy[i] = sqrt(1/rgamma(1,
                           shape_sdy,
                           (B_sdy+0.5*sum((yestimate[i,]-pred)^2))))


    #https://stats.stackexchange.com/questions/525674/gibbs-sampler-for-normal-and-inverse-gamma-distribution-in-r
    # https://stats.stackexchange.com/questions/266665/gibbs-sampler-examples-in-r

    # https://statswithr.github.io/book/inference-and-decision-making-with-multiple-parameters.html

    # skewed!!!
    # https://arxiv.org/pdf/1305.3080.pdf
    # p.7 explains the above https://www.researchgate.net/publication/340400073_Bayesian_Inference_for_Skew-Symmetric_Distributions

    # look into this for skew normal Gibbs https://academic.oup.com/biostatistics/article/11/2/317/268224
    # and in the Suppl. M.:
    # implementation: https://scristia.github.io/RcppComputingClub/sn_mix.html
    #                 https://github.com/scristia/SkewNormalMix/blob/master/R/sn.gibbs.mix.R

    # Outgrowing the Procrustean Bed of Normality: The Utility of Bayesian Modeling for Asymmetrical Data Analysis

    ## 4. Metropolis-Hastings step to estimate the regression coefficients
    proposal_coeff = MH_propose_coeff(coefficients[i-1,],prop_sd =  prop_sd_coeff) # new proposed values


    if(any(proposal_coeff[4] <= 0) | i == 2) {
      HR = 0
      } else {# Q needs to be >0
      # Hastings ratio of the proposal
      logpostold = 0

      if(n_p != 0) logpostold <- logpostold + logposterior_norm(x = x[n_p_ind], yest = yestimate[i,n_p_ind], ymean = yobs_mean,
                                                  sdyest = sdyest[i,], coeff = coefficients[i-1,],
                                                  sdy = sdy[i])
      if(n_norm != 0) logpostold <- logpostold + logposterior_norm(x = x[n_norm_ind], yest = yestimate[i,n_norm_ind], ymean = ynorm_mu,
                                                     sdyest = ynorm_sd, coeff = coefficients[i-1,],
                                                     sdy = sdy[i])
      if(n_skew != 0) logpostold <- logpostold +  logposterior_skew(x = x[n_skew_ind], yest = yestimate[i,n_skew_ind], mu = yskew_mu, yskew_sigma, yskew_lambda,
                                                      coeff = coefficients[i-1,], sdy[i])

      logpostold = logpostold + logprior(coefficients[i-1,])

      logpostnew = 0

      if(n_p != 0) logpostnew <- logpostnew + logposterior_norm(x = x[n_p_ind], yest = yestimate[i,n_p_ind], ymean = yobs_mean,
                                                  sdyest = sdyest[i,], coeff = proposal_coeff,
                                                  sdy = sdy[i])
      if(n_norm != 0) logpostnew <- logpostnew + logposterior_norm(x = x[n_norm_ind], yest = yestimate[i,n_norm_ind], ymean = ynorm_mu,
                                                     sdyest = ynorm_sd, coeff = proposal_coeff,
                                                     sdy = sdy[i])
      if(n_skew != 0) logpostnew <- logpostnew +  logposterior_skew(x = x[n_skew_ind], yest = yestimate[i,n_skew_ind], mu = yskew_mu, yskew_sigma, yskew_lambda,
                                                      coeff = proposal_coeff, sdy[i])
      logpostnew = logpostnew + logprior(proposal_coeff)

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
                           DKA = coefficients[,2],
                           M = coefficients[,3],
                           Q = coefficients[,4],
                           sdy = sdy,
                           logpost = logpost),
                yestimate = yestimate,
                sdyest = sdyest,
                lat = x)
  return(output)
}

### Next steps: implement uniform distributions, maybe with this accept-reject step:
### https://arxiv.org/pdf/0907.4010.pdf

nbin = 1
ndbin = 1
nsbin = 1
#prop_sd_yest <- matrix(0.01,nrow = nbin, ncol = nbin)
#diag(prop_sd_yest) <- 0.5
prop_sd_coeff <- c(3,3,3,0.1)

coeff_inits = c(0,30,45,0.1)
yest_inits = rep(25,nbin+ndbin+nsbin) #c(30,25,20,15,10,0,0)
sdyest_inits = rep(2,nbin)

x <- c(60,50,10)

yobs <- rnorm(3,10,2)
ylist <- list(yobs)
ynorm_mu <- 20
ynorm_sd <- 2

yskew_mu <- 32 # c(35,32,20)
yskew_sigma <- 9.3 # c(1,2,3)
yskew_lambda <- 2.00000000000000 #c(4,3,5)


xseq <- seq(15,45,0.1)
s = 1
plot(xseq,2/sd1*dnorm((xseq-yskew_mu[s])/yskew_sigma[s])*pnorm(yskew_lambda[s]*(xseq-yskew_mu[s])/yskew_sigma[s]),
       type = "l", lty = 1, lwd = 2, col = rgb(0,0,1,0.5))
s = 2
points(xseq,2/sd1*dnorm((xseq-yskew_mu[s])/yskew_sigma[s])*pnorm(yskew_lambda[s]*(xseq-yskew_mu[s])/yskew_sigma[s]),
     type = "l", lty = 1, lwd = 2, col = rgb(0,.8,.8,0.5))

s = 3
points(xseq,2/sd1*dnorm((xseq-yskew_mu[s])/yskew_sigma[s])*pnorm(yskew_lambda[s]*(xseq-yskew_mu[s])/yskew_sigma[s]),
       type = "l", lty = 1, lwd = 2, col = rgb(0,.9,0,0.5))


#nbin = 4
#x <- c(30,40,50,60,10,20,70)
#npb <- c(5,5,5,3)
#ym <- (c(25,18,15,25))
#ysd <- c(1,1,1,5)
#y <- lapply(1:nbin, function(x) rnorm(npb[x],ym[x],ysd[x]))
#y[which(npb==0)] <- NA
plot(0,0,xlim = c(0,90), ylim = c(-1,40), type = "n",xlab = "latitude", ylab = "temperature")
points(rep(x[1],length(yobs)),yobs)
points(x[2],ynorm_mu, col = "blue")
points(rep(x[2],2),ynorm_mu+c(-2,2)*ynorm_sd, col = "blue", type = "l")
points(x[3],yskew_mu, col = "red")

####################################################################
####
####  TEST 1
####
####################################################################


obsmat <- data.frame(cbind(c(15,15,15,50,50,50),c(1,1,1,2,2,2),c(28,33,34,18,14,22)))# column latitude, sample, temperature
colnames(obsmat) <- c("latitude", "sample", "temperature")
distrmat <- data.frame(cbind(c(5,30),c("normal","skew-normal"),c(35,30), c(4,4), c(NA,5)))# column latitude, distribution, location, scale, shape
colnames(distrmat) <- c("latitude", "distribution", "location", "scale", "shape")

prop_sd_coeff <- c(3,3,3,0.1)
coeff_inits = c(0,30,45,0.1)
yest_inits = rep(25,length(unique(obsmat$sample))+nrow(distrmat)) #c(30,25,20,15,10,0,0)
sdyest_inits = rep(2,length(unique(obsmat$sample)))


nIter = 50000
sdy_init = 1

system.time({ m1 <-  run_MCMC(nIter = nIter, obsmat = obsmat, distrmat = distrmat,
                              coeff_inits = coeff_inits,
                              sdy_init = sdy_init, yest_inits = yest_inits,
                              sdyest_inits = sdyest_inits,
                              prop_sd_coeff=prop_sd_coeff)
})

burnin = 25000+1

nIter = nIter
x = x
ylist = ylist
ynorm_mu = ynorm_mu
ynorm_sd = ynorm_sd
yskew_mu  = yskew_mu
yskew_sigma = yskew_sigma
yskew_lambda = yskew_lambda
coeff_inits = coeff_inits
sdy_init = sdy_init
yest_inits = yest_inits
sdyest_inits = sdyest_inits
prop_sd_coeff=prop_sd_coeff


latitude <- 0:90
sample_it <- sample((burnin+1):nIter,2000)
grad_025 <- sapply(1:length(latitude), function(f) quantile(apply(m1[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.025))
grad_975 <- sapply(1:length(latitude), function(f) quantile(apply(m1[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.975))

plot(seq(0,90,0.1), gradient(seq(0,90,0.1), apply(m1[[1]][burnin:nIter,1:4],2,median), 0),
     ylim = c(-1,60), type = "l", lwd = 3, ylab = "Temperature", xlab = "Latitude")
### Calculate confidence intervals
error_polygon(latitude,grad_025,grad_975,rgb(0,0,0,0.15))

for(i in 1:5) points(latitude,gradient(latitude,unlist(m1[[1]][sample(burnin:nIter,1),1:4]),0), type = "l", lwd = 2, col = rgb(0,0,0,0.3))


index <- 1

for(d in index) beanplot::beanplot(rnorm(2000,as.numeric(distrmat[d,3]),as.numeric(distrmat[d,4])), add = T,
                                 at = as.numeric(distrmat[d,1]), maxwidth = 5, side = "second",
                                 what = c(0,1,0,0), col = rgb(0,0.7,0.7,0.5), border = NA)

z1 <- list(NULL)
y1 <- list(NULL)
yskew_rho <-  as.numeric(distrmat[,5])/sqrt(1+as.numeric(distrmat[,5])^2)
N = 2000
i1 <- 0
for(d in 2) {
i1=i1+1
z1[[i1]] <- truncnorm::rtruncnorm(n = N, a = 0, b = Inf, mean = 0, sd = 1)
#y1 <- epsilon1 + a1*z1 + rnorm(N,0,omega1)
y1[[i1]] <- as.numeric(distrmat[d,3]) + as.numeric(distrmat[d,4])*yskew_rho[d]*z1[[i1]] +
  as.numeric(distrmat[d,4])*sqrt(1-(yskew_rho[d]^2))*rnorm(N,0,1)

beanplot::beanplot(y1[[i1]], add = T,
                   at = as.numeric(distrmat[d,1]), maxwidth = 5, side = "second",
                   what = c(0,1,0,0), col = rgb(0,0.7,0.7,0.5), border = NA)
}

for(i in unique(obsmat$sample)) {
  subs <- subset(obsmat, sample == i)
  points(subs$latitude,subs$temperature, pch = 4, bg = NA, col = rgb(0,0.5,0.85,0.67), cex = 2, lwd = 2)}

np <- ncol(m1[[2]])
offset <- rep(0,np)
points(m1$lat+offset,apply(m1[[2]][burnin:nIter,1:np],2,mean), col = rgb(1,0,0,0.5), pch = 19, cex = 1.25)
sapply(1:np, function(a) points(c(m1$lat[a]+offset[a],m1$lat[a]+offset[a]),quantile(m1[[2]][burnin:nIter,a], probs = c(0.25,0.75)), type = "l", col= rgb(1,0,0,0.5), lwd =5))
sapply(1:np, function(a) points(c(m1$lat[a]+offset[a],m1$lat[a]+offset[a]),quantile(m1[[2]][burnin:nIter,a], probs = c(0.125,0.875)), type = "l", col= rgb(1,0,0,0.5), lwd =3))
sapply(1:np, function(a) points(c(m1$lat[a]+offset[a],m1$lat[a]+offset[a]),quantile(m1[[2]][burnin:nIter,a], probs = c(0.025,0.975)), type = "l", col= rgb(1,0,0,0.5), lwd =1))



####################################################################
####
####  TEST 2
####
####################################################################


obsmat <- data.frame(cbind(c(15,15,15,50,50,50,25,25,25,25),c(1,1,1,2,2,2,3,3,3,3),c(45,38,34,18,17,26,30,38,41,44)))# column latitude, sample, temperature
colnames(obsmat) <- c("latitude", "sample", "temperature")
distrmat <- data.frame(cbind(c(5,30,35),c("normal","skew-normal","skew-normal"),c(35,30,25), c(4,4,3), c(NA,5,-7)))# column latitude, distribution, location, scale, shape
colnames(distrmat) <- c("latitude", "distribution", "location", "scale", "shape")

prop_sd_coeff <- c(3,3,3,0.1)
coeff_inits = c(0,30,45,0.1)
yest_inits = rep(25,length(unique(obsmat$sample))+nrow(distrmat)) #c(30,25,20,15,10,0,0)
sdyest_inits = rep(2,length(unique(obsmat$sample)))


nIter = 50000
sdy_init = 1

system.time({ m1 <-  run_MCMC(nIter = nIter, obsmat = obsmat, distrmat = distrmat,
                              coeff_inits = coeff_inits,
                              sdy_init = sdy_init, yest_inits = yest_inits,
                              sdyest_inits = sdyest_inits,
                              prop_sd_coeff=prop_sd_coeff)
})

burnin = 25000+1


latitude <- 0:90
sample_it <- sample((burnin):nIter,2000)
grad_025 <- sapply(1:length(latitude), function(f) quantile(apply(m1[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.025))
grad_975 <- sapply(1:length(latitude), function(f) quantile(apply(m1[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.975))

plot(seq(0,90,0.1), gradient(seq(0,90,0.1), apply(m1[[1]][burnin:nIter,1:4],2,median), 0),
     ylim = c(-1,60), type = "l", lwd = 3, ylab = "Temperature", xlab = "Latitude")
### Calculate confidence intervals
error_polygon(latitude,grad_025,grad_975,rgb(0,0,0,0.15))

for(i in 1:5) points(latitude,gradient(latitude,unlist(m1[[1]][sample(burnin:nIter,1),1:4]),0), type = "l", lwd = 2, col = rgb(0,0,0,0.3))


index <- 1

for(d in index) beanplot::beanplot(rnorm(2000,as.numeric(distrmat[d,3]),as.numeric(distrmat[d,4])), add = T,
                                   at = as.numeric(distrmat[d,1]), maxwidth = 5, side = "second",
                                   what = c(0,1,0,0), col = rgb(0,0.7,0.7,0.6), border = NA)

z1 <- list(NULL)
y1 <- list(NULL)
yskew_rho <-  as.numeric(distrmat[,5])/sqrt(1+as.numeric(distrmat[,5])^2)
N = 2000
i1 <- 0
for(d in 2:3) {
  i1=i1+1
  z1[[i1]] <- truncnorm::rtruncnorm(n = N, a = 0, b = Inf, mean = 0, sd = 1)
  #y1 <- epsilon1 + a1*z1 + rnorm(N,0,omega1)
  y1[[i1]] <- as.numeric(distrmat[d,3]) + as.numeric(distrmat[d,4])*yskew_rho[d]*z1[[i1]] +
    as.numeric(distrmat[d,4])*sqrt(1-(yskew_rho[d]^2))*rnorm(N,0,1)

  beanplot::beanplot(y1[[i1]], add = T,
                     at = as.numeric(distrmat[d,1]), maxwidth = 5, side = "second",
                     what = c(0,1,0,0), col = rgb(0,0.7,0.7,0.6), border = NA)
}

for(i in unique(obsmat$sample)) {
  subs <- subset(obsmat, sample == i)
  points(subs$latitude,subs$temperature, pch = 4, bg = NA, col = rgb(0,0.5,0.85,0.67), cex = 2, lwd = 2)}

np <- ncol(m1[[2]])
offset <- rep(0,np)
points(m1$lat+offset,apply(m1[[2]][burnin:nIter,1:np],2,mean), col = rgb(1,0,0,0.5), pch = 19, cex = 1.25)
sapply(1:np, function(a) points(c(m1$lat[a]+offset[a],m1$lat[a]+offset[a]),quantile(m1[[2]][burnin:nIter,a], probs = c(0.25,0.75)), type = "l", col= rgb(1,0,0,0.5), lwd =5))
sapply(1:np, function(a) points(c(m1$lat[a]+offset[a],m1$lat[a]+offset[a]),quantile(m1[[2]][burnin:nIter,a], probs = c(0.125,0.875)), type = "l", col= rgb(1,0,0,0.5), lwd =3))
sapply(1:np, function(a) points(c(m1$lat[a]+offset[a],m1$lat[a]+offset[a]),quantile(m1[[2]][burnin:nIter,a], probs = c(0.025,0.975)), type = "l", col= rgb(1,0,0,0.5), lwd =1))





####################################################################
####
####  TEST 3
####
####################################################################


obsmat <- data.frame(cbind(c(15,15,15,50,50,50,25,25,25),c(1,1,1,2,2,2,3,3,3),c(45,38,34,18,17,16,30,32,34)))# column latitude, sample, temperature
colnames(obsmat) <- c("latitude", "sample", "temperature")
distrmat <- data.frame(cbind(c(5,30,35),c("normal","skew-normal","skew-normal"),c(15,27,30), c(5,4,3), c(NA,5,-7)))# column latitude, distribution, location, scale, shape
colnames(distrmat) <- c("latitude", "distribution", "location", "scale", "shape")

prop_sd_coeff <- c(3,3,3,0.1)
coeff_inits = c(0,30,45,0.1)
yest_inits = rep(25,length(unique(obsmat$sample))+nrow(distrmat)) #c(30,25,20,15,10,0,0)
sdyest_inits = rep(2,length(unique(obsmat$sample)))


nIter = 50000
sdy_init = 1

system.time({ m1 <-  run_MCMC(nIter = nIter, obsmat = obsmat, distrmat = distrmat,
                              coeff_inits = coeff_inits,
                              sdy_init = sdy_init, yest_inits = yest_inits,
                              sdyest_inits = sdyest_inits,
                              prop_sd_coeff=prop_sd_coeff)
})

burnin = 25000+1


latitude <- 0:90
sample_it <- sample((burnin):nIter,2000)
grad_025 <- sapply(1:length(latitude), function(f) quantile(apply(m1[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.025))
grad_975 <- sapply(1:length(latitude), function(f) quantile(apply(m1[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.975))

plot(seq(0,90,0.1), gradient(seq(0,90,0.1), apply(m1[[1]][burnin:nIter,1:4],2,median), 0),
     ylim = c(-1,60), type = "l", lwd = 3, ylab = "Temperature", xlab = "Latitude")
### Calculate confidence intervals
error_polygon(latitude,grad_025,grad_975,rgb(0,0,0,0.15))

for(i in 1:5) points(latitude,gradient(latitude,unlist(m1[[1]][sample(burnin:nIter,1),1:4]),0), type = "l", lwd = 2, col = rgb(0,0,0,0.3))


index <- 1

for(d in index) beanplot::beanplot(rnorm(2000,as.numeric(distrmat[d,3]),as.numeric(distrmat[d,4])), add = T,
                                   at = as.numeric(distrmat[d,1]), maxwidth = 5, side = "second",
                                   what = c(0,1,0,0), col = rgb(0,0.7,0.7,0.6), border = NA)

z1 <- list(NULL)
y1 <- list(NULL)
yskew_rho <-  as.numeric(distrmat[,5])/sqrt(1+as.numeric(distrmat[,5])^2)
N = 2000
i1 <- 0
for(d in 2:3) {
  i1=i1+1
  z1[[i1]] <- truncnorm::rtruncnorm(n = N, a = 0, b = Inf, mean = 0, sd = 1)
  #y1 <- epsilon1 + a1*z1 + rnorm(N,0,omega1)
  y1[[i1]] <- as.numeric(distrmat[d,3]) + as.numeric(distrmat[d,4])*yskew_rho[d]*z1[[i1]] +
    as.numeric(distrmat[d,4])*sqrt(1-(yskew_rho[d]^2))*rnorm(N,0,1)

  beanplot::beanplot(y1[[i1]], add = T,
                     at = as.numeric(distrmat[d,1]), maxwidth = 5, side = "second",
                     what = c(0,1,0,0), col = rgb(0,0.7,0.7,0.6), border = NA)
}

for(i in unique(obsmat$sample)) {
  subs <- subset(obsmat, sample == i)
  points(subs$latitude,subs$temperature, pch = 4, bg = NA, col = rgb(0,0.5,0.85,0.67), cex = 2, lwd = 2)}

np <- ncol(m1[[2]])
offset <- rep(0,np)
points(m1$lat+offset,apply(m1[[2]][burnin:nIter,1:np],2,mean), col = rgb(1,0,0,0.5), pch = 19, cex = 1.25)
sapply(1:np, function(a) points(c(m1$lat[a]+offset[a],m1$lat[a]+offset[a]),quantile(m1[[2]][burnin:nIter,a], probs = c(0.25,0.75)), type = "l", col= rgb(1,0,0,0.5), lwd =5))
sapply(1:np, function(a) points(c(m1$lat[a]+offset[a],m1$lat[a]+offset[a]),quantile(m1[[2]][burnin:nIter,a], probs = c(0.125,0.875)), type = "l", col= rgb(1,0,0,0.5), lwd =3))
sapply(1:np, function(a) points(c(m1$lat[a]+offset[a],m1$lat[a]+offset[a]),quantile(m1[[2]][burnin:nIter,a], probs = c(0.025,0.975)), type = "l", col= rgb(1,0,0,0.5), lwd =1))


####################################################################
####
####  TEST with Cambrian (Stage 3) ooid data
####
####################################################################

ooid <- read.csv("D://OneDrive - Durham University/projects/ooids/ooid_data_set_samples_processed.csv")

datum1 <- subset(ooid, age_ori == "Cambrian\n(Series 2/Atdabanian)" & location == "Northeastern Armorican\nMassif, France")
datum2 <- subset(ooid, age_ori == "Cambrian\n(Series 2/Atdabanian)" & location == "Flinders Range\nSouth Australia")
coord1 <- c(3,48) # "Amorican massif
coord2 <- c(138.7, -31.5) # Flinders range

ooidsubs <- rbind(datum1,datum2)

pallat_ooid <- c(abs(c(-46.34,16.42)))

distrmat <- data.frame(cbind(pallat_ooid,c("skew-normal","skew-normal"),ooidsubs$mu, ooidsubs$sigma, ooidsubs$lambda))# column latitude, distribution, location, scale, shape
colnames(distrmat) <- c("latitude", "distribution", "location", "scale", "shape")

isot <- read.csv(
  "D://OneDrive - Durham University/projects/isotopes/GrossmanJoachimski2022SI/Isotope_Grossmann_Joachimski_ScienfiticReports2022.csv",
  fileEncoding = "latin1")
isots <- subset(isot,stage_2020 == "Age 3")

### Hays no trend, seawater d18 O of -1.08:
isots_t <- 15.7 - 4.36*(isots$sauerstof_isotop_permille - isots$d18Osw..iceV) + 0.12*(isots$sauerstof_isotop_permille-0.6-isots$d18Osw..iceV)^2
pallat_iso <- rep(abs(-0.92),length(isots_t))
obsmat <- data.frame(cbind(pallat_iso,rep(1,nrow(isots)), isots_t))# column latitude, sample, temperature
colnames(obsmat) <- c("latitude", "sample", "temperature")

prop_sd_coeff <- c(3,3,3,0.1)
coeff_inits = c(0,30,45,0.1)
yest_inits = rep(25,length(unique(obsmat$sample))+nrow(distrmat)) #c(30,25,20,15,10,0,0)
sdyest_inits = rep(2,length(unique(obsmat$sample)))


nIter = 50000
sdy_init = 1

system.time({ m1 <-  run_MCMC(nIter = nIter, obsmat = obsmat, distrmat = distrmat,
                              coeff_inits = coeff_inits,
                              sdy_init = sdy_init, yest_inits = yest_inits,
                              sdyest_inits = sdyest_inits,
                              prop_sd_coeff=prop_sd_coeff)
})

burnin = 25000+1


latitude <- 0:90
sample_it <- sample((burnin):nIter,2000)
grad_025 <- sapply(1:length(latitude), function(f) quantile(apply(m1[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.025))
grad_975 <- sapply(1:length(latitude), function(f) quantile(apply(m1[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.975))

plot(seq(0,90,0.1), gradient(seq(0,90,0.1), apply(m1[[1]][burnin:nIter,1:4],2,median),0), xlim = c(0,90),
     ylim = c(-1,65), type = "l", lwd = 3, ylab = "Temperature", xlab = "Latitude", "Cambrian - Stage 3")
### Calculate confidence intervals
error_polygon(latitude,grad_025,grad_975,rgb(0,0,0,0.15))

for(i in 1:5) points(latitude,gradient(latitude,unlist(m1[[1]][sample(burnin:nIter,1),1:4]),0), type = "l", lwd = 2, col = rgb(0,0,0,0.3))


#index <- 1

#for(d in index) beanplot::beanplot(rnorm(2000,as.numeric(distrmat[d,3]),as.numeric(distrmat[d,4])), add = T,
#                                   at = as.numeric(distrmat[d,1]), maxwidth = 5, side = "second",
#                                   what = c(0,1,0,0), col = rgb(0,0.7,0.7,0.6), border = NA)

z1 <- list(NULL)
y1 <- list(NULL)
yskew_rho <-  as.numeric(distrmat[,5])/sqrt(1+as.numeric(distrmat[,5])^2)
N = 4000
i1 <- 0
beancols <- c( rgb(0.40,0.8,0,0.6), rgb(0,0.7,0.7,0.6))
for(d in 1:2) {
  i1=i1+1
  z1[[i1]] <- truncnorm::rtruncnorm(n = N, a = 0, b = Inf, mean = 0, sd = 1)
  #y1 <- epsilon1 + a1*z1 + rnorm(N,0,omega1)
  y1[[i1]] <- as.numeric(distrmat[d,3]) + as.numeric(distrmat[d,4])*yskew_rho[d]*z1[[i1]] +
    as.numeric(distrmat[d,4])*sqrt(1-(yskew_rho[d]^2))*rnorm(N,0,1)

  beanplot::beanplot(y1[[i1]], add = T,
                     at = as.numeric(distrmat[d,1]), maxwidth = 5, side = "second",
                     what = c(0,1,0,0), col = beancols[d], border = NA)
}

for(i in unique(obsmat$sample)) {
  subs <- subset(obsmat, sample == i)
  points(subs$latitude,subs$temperature, pch = 4, bg = NA, col = rgb(0,0.5,0.85,0.67), cex = 2, lwd = 2)}

np <- ncol(m1[[2]])
offset <- rep(0,np)
points(m1$lat+offset,apply(m1[[2]][burnin:nIter,1:np],2,mean), col = rgb(1,0,0,0.5), pch = 19, cex = 1.25)
sapply(1:np, function(a) points(c(m1$lat[a]+offset[a],m1$lat[a]+offset[a]),quantile(m1[[2]][burnin:nIter,a], probs = c(0.25,0.75)), type = "l", col= rgb(1,0,0,0.5), lwd =5))
sapply(1:np, function(a) points(c(m1$lat[a]+offset[a],m1$lat[a]+offset[a]),quantile(m1[[2]][burnin:nIter,a], probs = c(0.125,0.875)), type = "l", col= rgb(1,0,0,0.5), lwd =3))
sapply(1:np, function(a) points(c(m1$lat[a]+offset[a],m1$lat[a]+offset[a]),quantile(m1[[2]][burnin:nIter,a], probs = c(0.025,0.975)), type = "l", col= rgb(1,0,0,0.5), lwd =1))

legend("topright",legend=c("regression line", "temperature estimate","dO18 data", "aragonitic giant ooids", "calcitic ooids"),
       cex = 0.8, col = c("black", "red", rgb(0,0.5,0.85,.67),  beancols[1], beancols[2]), lwd = c(2,NA,NA,6,6),
       pch = c(NA,19,4,NA,NA), pt.cex = c(NA,1.25,2,NA,NA))
