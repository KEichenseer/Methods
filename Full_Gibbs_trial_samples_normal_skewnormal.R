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
  ll1 <- sum(dnorm(yest, ymean, sdyest,log=TRUE),na.rm=T)
  # likelihood for yest - other PDFs
  # likelihood for
  #ll1b <- sum(dunif(yest[1],20,40, log = TRUE),
  #            dnorm(yest[6],6,5, log = TRUE))
  pred = A + max(c(K-A,0))/((1+(exp(Q*(x-M)))))
  ll2 <- sum(dnorm(yest, mean = pred, sd = sdy, log = TRUE))
  return(c(ll2+ll1))
}

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

# Gibbs sampling of mu with skew normal likelihood and normal prior
skew_mu <- function(x, y, sigma, rho, mu_prior, sigma_prior) {

  rnorm(1,(n*sigma_prior^2*(mean(x)-sigma*rho*mean(y))+sigma^2*(1-rho^2)*mu_prior)/
          (n*sigma_prior^2+sigma^2*(1-rho^2)),
        sqrt((sigma_prior^2*sigma^2*(1-rho^2))/(n*sigma_prior^2+sigma^2*(1-rho^2))) )
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
                     prop_sd_coeff, prop_sd_yest, quiet = FALSE){
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

  # start progress bar
  if (!quiet) cli::cli_progress_bar('Sampling', total = nIter)
  ### The MCMC loop
  for (i in 2:nIter){
    # update progress bar
    #if (!quiet) cli::cli_progress_update(set = i, status = paste0('iteration ', i))

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

    ### 1.3. Gibbs steps for data that have location, scale and shape parameter given (skew-normal). Estimate global mean only)
    ## 1.3.a Gibbs step to estimate yestimate
    ### distribution from here: http://koreascience.or.kr/article/JAKO200504840590864.pdf
    z <- (yskew_mu - yestimate[i-1,n_skew_ind])/yskew_sigma
    y <- truncnorm::rtruncnorm(n,0,Inf,rho*z,sqrt(1-rho^2))
    yestimate[i,n_skew_ind] = skew_mu(x=yskew_mu[i-1,], y=y, sigma=yskew_sigma, rho=yskew_rho,
                                      mu_prior=pred[n_d_ind], sigma_prior=sdy[i-1])


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
