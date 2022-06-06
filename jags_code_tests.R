regression_model_with_errors <- function() {

  ## Likelihood
  for (i in 1:n){
    y_measurement[i] ~ dnorm(y_estimate[i], 1/(y_se[i]*y_se[i]))
    y_estimate[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta * x[i]
  }


  ## Priors
  sigma ~ dnorm(0, 1/(100*100)); T(0,)
  tau <- 1 / (sigma * sigma) # precision = 1/standard deviation ^2
  alpha ~ dnorm(0, 1/(100*100))
  beta ~ dnorm(0, 1/(100*100))

}


n <- 100
x <- rnorm(n,0,1)
y_measurement <- rnorm(n,x,1)
y_se <- rnorm(n,0,0.2)
plot(x, y_measurement, ylim = c(-3.3,5))
for(i in 1:n) points(rep(x[i],2), y_measurement[i]+c(-2,2)*y_se[i], type = "l")
regression_data <- list("x", "y_measurement", "y_se","n")

fit_se  <- R2jags::jags(data = regression_data,
                parameters.to.save = c("alpha",
                                       "beta",
                                       "sigma",
                                       "y_estimate"
                ),
                n.iter = 2000,
                n.thin = 1,
                n.chains =  3, # Other values set at default (for simplicity)
                model.file = regression_model_with_errors)

#plot(x,y, pch = 21, bg = rgb(0,0,0,0.2))

abline(a = fit$BUGSoutput$mean$alpha, b = fit$BUGSoutput$mean$beta, lwd = 2)

regrange <- seq(min(x), max(x), 0.1)
regmat2 <- matrix(NA,nrow = 1000, ncol = length(regrange))

for(i in 1:1000) {
  regmat2[i,] <- fit_se$BUGSoutput$sims.list$alpha[i] + fit_se$BUGSoutput$sims.list$beta[i]*regrange
}

reg_025se <- apply(regmat2, 2, function(x) quantile(x, probs = 0.025))
reg_975se <- apply(regmat2, 2, function(x) quantile(x, probs = 0.975))

plot(x, y_measurement, ylim = c(-3.3,5))
for(i in 1:n) points(rep(x[i],2), y_measurement[i]+c(-2,2)*y_se[i], type = "l")
error_polygon(regrange,reg_975se,reg_025se, col = rgb(1,0,0,0.2))

points(regrange, (fit_se$BUGSoutput$mean$alpha + regrange*fit_se$BUGSoutput$mean$beta), type = "l", lwd = 2, col = "red")

points(x,fit_se$BUGSoutput$mean$y_estimate, bg = rgb(1,0,0,0.4), pch = 21, col = NA)



regression_model_with_points <- function() {

  ## Likelihood
  for (j in 1:n){
  for (i in 1:ni[j]){
    y[i,j] ~ dnorm(y_estimate[j], 1/(y_se[j]*y_se[j]))
  }
    y_estimate[j] ~ dnorm(mu[j], tau)
    mu[j] <- alpha + beta * x[j]

    y_se[j] ~ dgamma(1,1)

  }


  ## Priors
  sigma ~ dnorm(0, 1/(100*100)); T(0,)
  tau <- 1 / (sigma * sigma) # precision = 1/standard deviation ^2
  alpha ~ dnorm(0, 1/(100*100))
  beta ~ dnorm(0, 1/(100*100))

}


n <- 6
x <- rnorm(n,0,1)
ni <- c(1,5,5,5,5,5)+2
y <- matrix(NA,nrow = max(ni), ncol = n)
for(i in 1:n) {
  y[1:ni[i],i] <- rnorm(ni[i],x[i],0.2)
}

y[1,1] <- 1

plot(x,y[1,], ylim = range(y,na.rm = T), xlim = range(x), col = rgb(1:n/n,n:1/n,0,1), pch = 19)
points(x,y[2,], ylim = range(y,na.rm = T), xlim = range(x), col = rgb(1:n/n,n:1/n,0,1), pch = 19)
points(x,y[3,], ylim = range(y,na.rm = T), xlim = range(x), col = rgb(1:n/n,n:1/n,0,1), pch = 19)
points(x,y[4,], ylim = range(y,na.rm = T), xlim = range(x), col = rgb(1:n/n,n:1/n,0,1), pch = 19)
points(x,y[5,], ylim = range(y,na.rm = T), xlim = range(x), col = rgb(1:n/n,n:1/n,0,1), pch = 19)


regression_data <- list("x", "y","n", "ni")

fit_p  <- R2jags::jags(data = regression_data,
                        parameters.to.save = c("alpha",
                                               "beta",
                                               "sigma",
                                               "y_estimate",
                                               "y_se"
                        ),
                        n.iter = 50000,
                        n.thin = 1,
                        n.chains =  1, # Other values set at default (for simplicity)
                        model.file = regression_model_with_points)

abline(a = fit_p$BUGSoutput$mean$alpha, b = fit_p$BUGSoutput$mean$beta, lwd = 2)
library("ggmcmc")
ggs_traceplot(ggs(coda::as.mcmc.list(fit_p$BUGSoutput)),family = "y_se")
hist(fit_p$BUGSoutput$sims.list$y_se[,2],100)
hist(rgamma(10000,1,1), 100, add = T,col = rgb(1,0,0,0.2))



# Main MCMCM function
run_MCMC <- function(nIter, x, yobs, prior_df, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                     prop_sd_coeff, prop_sd_yest, nbin){
  ### Initialisation
  coefficients = array(dim = c(nIter,4)) # set up array to store coefficients
  coefficients[1,] = coeff_inits # initialise coefficients
  sdy = rep(NA_real_,nIter) # set up vector to store sdy
  sdy[1] = sdy_init # intialise sdy

  yestimate = array(dim = c(nIter,nbin)) # set up array to store coefficients
  yestimate[1,] = yest_inits # initialise coefficients

  sdyest = array(dim = c(nIter,nbin)) # set up vector to store sdy
  sdyest[1,] = sdyest_inits # intialise sdy

  A_sdy = 3 # parameter for the prior on the inverse gamma distribution of sdy
  B_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  shape_sdy <- A_sdy+nbin/2 # shape parameter for the inverse gamma


  #### Investigate these: Need to be broad for single obser
  A_sdyest = 3 # parameter for the prior on the inverse gamma distribution of sdyest
  B_sdyest = 1 # parameter for the prior on the inverse gamma distribution of sdyest
  ####
  nest = sapply(yobs,length)
  nest[which(is.na(yobs))] = NA
  shape_sdyest =  A_sdyest+nest/2 # shape parameter for the inverse gamma
  ### n-1?!
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
      1,max(c(shape_sdyest[j],A_sdyest),na.rm=T),
      max(c(B_sdyest+0.5*sum((yobs[[j]]-yestimate[i-1,j])^2,na.rm=T), B_sdyest),na.rm=T)))
    #https://stats.stackexchange.com/questions/525674/gibbs-sampler-for-normal-and-inverse-gamma-distribution-in-r
    # https://stats.stackexchange.com/questions/266665/gibbs-sampler-examples-in-r

    ## 3. Metropolis-Hastings step to estimate yest
    proposal_yest = MH_propose_yest(yestimate[i-1,],prop_sd_yest =  prop_sd_yest) # new proposed values

    ## 4. Metropolis-Hastings step to estimate the regression coefficients
    proposal_coeff = MH_propose_coeff(coefficients[i-1,],prop_sd =  prop_sd_coeff) # new proposed values
    if(any(proposal_coeff[4] <= 0)) HR = 0 else {# Q needs to be >0
      # Hastings ratio of the proposal
      logpostold = logposterior(x = x, ymean = ymean, yest = yestimate[i-1,],
                                sdyest = sdyest[i,], coeff = coefficients[i-1,],
                                sdy = sdy[i], prior_df = prior_df)
      logpostnew = logposterior(x = x, ymean = ymean, yest = proposal_yest,
                                sdyest = sdyest[i,], coeff = proposal_coeff,
                                sdy = sdy[i], prior_df = prior_df)
      HR = exp(logpostnew -
                 logpostold)
    }
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


