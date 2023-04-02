ciPoly <- function(x,en,ep,color=rgb(0,0,0,0.2)) {
  polygon( c(x[1], x, x[length(x)], x[length(x)], rev(x), x[1]),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)}


set.seed(777)

N <- 100

x <- runif(N,2,5)
x_sd <- 0.2*x+0.25 #runif(N,0.025,0.25)
x_obs <- x + 1*rnorm(N,0,x_sd)

y <- 7+2*sin(x)
y_sd <- runif(N,0.1,0.2) # rep(1,N)#
y_obs <- y + rnorm(N,0,y_sd)


x_prec <- 1/(x_sd^2)
y_prec <- 1/(y_sd^2)


plot(x_obs,y_obs,pch = 21, bg = rgb(0,0,0,0.2))


#### Load packages
library(R2jags)
library(ggmcmc)


regression_model_with_errors <- function() {

  ## Likelihood
  for (i in 1:N){
    x_estimate[i] ~ dnorm(x_measurement[i], 1/(x_sd[i]*x_sd[i]))
    y_measurement[i] ~ dnorm(y_estimate[i], 1/(y_sd[i]*y_sd[i]))
    y_estimate[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta * x_estimate[i]
  }


  ## Priors

  tau ~ dgamma(1, 1)  # inverse-gamma prior for precision
  sigma <- 1/sqrt(tau)

  #tau <- 1/(sigma*sigma) # precision = 1/standard deviation ^2
  #sigma ~ dunif(0,1)
  alpha ~ dnorm(0, 1/(100*100))
  beta ~ dnorm(0, 1/(100*100))

}


y_measurement <- y_obs
x_measurement <- x_obs

regression_data <- list("x_measurement", "y_measurement", "x_sd", "y_sd","N")

fit_sd  <- jags(data = regression_data,
                parameters.to.save = c("alpha",
                                       "beta",
                                       "sigma",
                                       "y_estimate",
                                       "x_estimate"
                ),
                n.iter = 10000,
                n.thin = 1,
                n.chains =  3, # Other values set at default (for simplicity)
                model.file = regression_model_with_errors)

# points(fit_sd$BUGSoutput$mean$x_estimate,fit_sd$BUGSoutput$mean$y_estimate,
#        pch = 21, bg = rgb(1,0,0,0.2), col = rgb(1,0,0,0.75))
#
# abline(lm(y_obs~x_obs),lty=1,lwd=2)
# abline(lm(y~x),lty=2,lwd=2)
#


### Gibbs sampler

X <- cbind(rep(1,length(x)),x_obs)

# pre-compute
XtX = t(X) %*% X
beta_hat = solve(XtX, t(X) %*% y_obs) # compute this ahead of time cause we'll need it a lot
XtXi = solve(XtX)

n_params = 2
n_obs = N

beta = c(0,0) # starting value
sigma2 = 1 #starting value

n_iterations = 10000

beta_out = matrix(data=NA, nrow=n_iterations, ncol=n_params)
sigma2_out = matrix(data = NA, nrow = n_iterations, ncol=1)

for (i in 1:n_iterations){
  beta = mvnfast::rmvn(n=1, beta_hat, sigma2 * XtXi)[1,] # Beta is N(XtXiXty, sigma^2XtXi)

  part = (y_obs - X %*% beta)

  sigma2 = 1/rgamma(1, 1+n_obs/2, 1+t(part) %*% part * .5 ) # sigma^2 is IG(n/2, ....)

  # save the results.
  beta_out[i,] = beta
  sigma2_out[i,] = sigma2
}

# hist(beta_out[,1])
# abline(v=lm1$coefficients[1], col='red')
# summary(lm1)
plot(x_obs,y_obs,pch = 21, bg = rgb(0,0,0,0.2))

regrange <- seq(min(x_obs), max(x_obs), 0.1)
regmat <- matrix(NA,nrow = 15000, ncol = length(regrange))

for(i in 1:15000) {
  regmat[i,] <- fit_sd$BUGSoutput$sims.list$alpha[i] + fit_sd$BUGSoutput$sims.list$beta[i]*regrange
}

reg_025 <- apply(regmat, 2, function(x) quantile(x, probs = 0.025))
reg_975 <- apply(regmat, 2, function(x) quantile(x, probs = 0.975))

ciPoly(regrange, reg_975,reg_025, col = rgb(1,.2,0,0.25))
abline(a = fit_sd$BUGSoutput$mean$alpha, b = fit_sd$BUGSoutput$mean$beta, col = "orange",lwd=3)



regrange <- seq(min(x_obs), max(x_obs), 0.1)
regmat <- matrix(NA,nrow = 8000, ncol = length(regrange))

for(i in 2001:10000) {
  regmat[i-2000,] <- beta_out[i,1] + beta_out[i,2]*regrange
}

reg_025 <- apply(regmat, 2, function(x) quantile(x, probs = 0.025))
reg_500 <- apply(regmat, 2, function(x) quantile(x, probs = 0.500))

reg_975 <- apply(regmat, 2, function(x) quantile(x, probs = 0.975))


lm1 <- lm(y_obs~x_obs)
pred <- predict(lm1, newdata = data.frame(x_obs=regrange),interval = "confidence")
ciPoly(regrange, pred[,2],pred[,3], col = rgb(0,0,0,0.2))
points(regrange, pred[,1], type = "l", lwd = 2)

ciPoly(regrange, reg_975,reg_025, col = rgb(0.75,0,0.65,0.2))
#sapply(1:n, function(i) points(rep(x[i], 2), c(y[i]+y_se[i],y[i]-y_se[i]), type = "l", col = rgb(0,0,0,0.4)))
points(regrange, reg_500, type = "l", lwd = 2, col = "red")

#ciPoly(regrange, reg_025,reg_975, col = rgb(1,0,0,0.2))




###############################################################################



1.28262    +c(-1,1)*0.28148
quantile(beta_out[,1],probs = c(0.025,0.975))

1.85284    +c(-1,1)*0.09275
quantile(beta_out[,2],probs = c(0.025,0.975))

formula <- bf(y_obs ~ x_obs)

### Gibbs sampler for errors-in-variables regression

X <- cbind(rep(1,length(x)),x_obs)

# pre-compute

n_params = 2
n_obs = N

beta = c(1,1) # starting value
sigma2 = 1 #starting value

X_est = X
y_est = y_obs
y_pred = (X_est %*% beta)[,1]

n_iterations = 20000

beta_out = matrix(data=NA, nrow=n_iterations, ncol=n_params)
sigma2_out = matrix(data = NA, nrow = n_iterations, ncol=1)

x_est_out = matrix(data=NA, nrow=n_iterations, ncol=N)
y_est_out = matrix(data=NA, nrow=n_iterations, ncol=N)

for (i in 1:n_iterations){



  # y_est = rnorm(n_obs,
  #               sigma2/(y_sd^2+sigma2)*y_obs + y_sd^2/(y_sd^2+sigma2)*y_pred,
  #               sqrt(1/(1/sigma2 + 1/y_sd^2)))

   y_est =      # rnorm(n_obs,
  #                    sigma2/(y_sd^2+sigma2)*y_pred + y_sd^2/(y_sd^2+sigma2)*y_obs,
  #                    sqrt(1/(1/sigma2 + 1/y_sd^2)))

  rnorm(n_obs,
        sigma2/(y_sd^2+sigma2)*y_obs + y_sd^2/(y_sd^2+sigma2)*y_pred,
        sqrt(1/(1/sigma2 + 1/y_sd^2))) # I think that one at least should be correct now... If x_sd is
        #  set to very small values, then the results mostly agree with those of the jags implementation

    #rnorm(n_obs,
  #               1/(1/sigma2+1/y_sd)*(y_pred/sigma2 + y_obs/y_sd),
  #               sqrt(1/(1/sigma2 + 1/y_sd^2)))


  if(beta[2] == 0) x_pred = rep(beta[1],length(x_obs)) else x_pred = (y_pred-beta[1])/beta[2]

  # term <- sigma2/beta[2]^2
  # if(term==Inf) term <- 0

   #X_est[,2] = rnorm(n_obs,
    #                 (beta[2] * (x_obs) / x_sd^2 + (y_est - beta[1] - y_pred) * beta[2] / sigma2) / (beta^2 * (1 / x_sd^2) + 1 / sigma2),
    #                 1 / (beta[2]^2 * sum(1 / x_sd^2) + 1 / sigma2)
     #                 )
   X_est[,2] = rnorm(n_obs,
                     ((beta[2] * x_obs / x_sd^2) + ((y_est - beta[1] - y_pred) * beta[2]/ sigma2)) / ((beta[2]^2 / x_sd^2) + (1 / sigma2)),
                     1 / ((beta[1]^2 / x_sd^2) + (1 / sigma2))
   )

     # rnorm(n_obs,
     #                 (beta * (y - alpha) + x_obs / x_sd^2) / (beta^2 + 1/x_sd^2), #(x_obs/x_sd^2 + beta[2]*y_est)/(1/x_sd^2 + beta[2]^2),
     #                 sqrt(1/(1/x_sd^2 + beta[2]^2)))

     #rnorm(n_obs,
  #                   mean = (y_obs - beta[1] - beta[2]*y_est) * (x_sd^2) / (y_sd^2 + x_sd^2) + x_obs * y_sd^2 / (y_sd^2 + x_sd^2),
  #                   sd = sqrt(x_sd^2 * y_sd^2 / (y_sd^2 + x_sd^2)))
  #
  # Save the current x_est sample
    #
    # rnorm(n_obs,
    #       sigma2/(x_sd^2+sigma2)*x_pred + x_sd^2/(x_sd^2+sigma2)*x_obs,
    #       sqrt(1/(1/sigma2 + 1/x_sd^2)))

  # rnorm(n_obs,
  #       sigma2/(x_sd^2+sigma2)*x_obs + x_sd^2/(x_sd^2+sigma2)*x_pred,
  #       sqrt(1/(1/sigma2 + 1/x_sd^2)))
  #
  #




    #
  XtX = t(X_est) %*% X_est
  beta_hat = solve(XtX, t(X_est) %*% y_est)
  XtXi = solve(XtX)

  beta = mvnfast::rmvn(n=1, beta_hat, sigma2 * XtXi)[1,] # Beta is N(XtXiXty, sigma^2XtXi)

  y_pred = (X_est %*% beta)[,1]

  resid = (y_pred-y_est)

  sigma2 = 1/rgamma(1, 1+n_obs/2, 1+sum(resid^2) * .5 ) # sigma^2 is IG(n/2, ....)

  # save the results.
  beta_out[i,] = beta
  sigma2_out[i,] = sigma2
  x_est_out[i,] = X_est[,2]
  y_est_out[i,] = y_est
}

#plot(x_obs,y_obs)
#abline(a=fit_sd$BUGSoutput$mean$alpha, b=fit_sd$BUGSoutput$mean$beta)
#plot(beta_out[1000:5000,1])
abline(a=median(beta_out[1000:20000,1]), b=median(beta_out[1000:20000,2]), col = "green",lwd=2, lty = 2)

#abline(a=mean(beta_samples[1000:10000,1]), b=mean(beta_samples[1000:10000,2]), col = "red")
regrange <- seq(min(x_obs), max(x_obs), 0.1)
regmat <- matrix(NA,nrow = 10000, ncol = length(regrange))

for(i in 10001:20000) {
  regmat[i-10000,] <- beta_out[i,1] + beta_out[i,2]*regrange
}

reg_025 <- apply(regmat, 2, function(x) quantile(x, probs = 0.025))
reg_975 <- apply(regmat, 2, function(x) quantile(x, probs = 0.975))
ciPoly(regrange, reg_025,reg_975, col = rgb(0,1,0,0.2))


points(apply(x_est_out,2,mean),apply(y_est_out,2,mean),pch = 21, col = NA, bg = rgb(0,0,1,0.332))

points(apply(fit_sd$BUGSoutput$sims.list$x_estimate,2,mean),
       apply(fit_sd$BUGSoutput$sims.list$y_estimate,2,mean),
       pch = 21, col = NA, bg = rgb(1,0,0.5,0.2))













lm1 <- lm(y_obs~x_obs)
pred <- predict(lm1, newdata = data.frame(x_obs=regrange),interval = "confidence")
ciPoly(regrange, pred[,2],pred[,3], col = rgb(0,.5,1,0.2))



quantile(beta_out[1000:20000,1],probs=c(0.025,0.5,0.975))
fit_sd$BUGSoutput$summary["alpha",]

quantile(beta_out[1000:20000,2],probs=c(0.025,0.5,0.975))
fit_sd$BUGSoutput$summary["beta",]


x_measurement <- x_obs
y_measurement <- y_obs

# Set the priors
alpha <- 0
beta <- 0
tau <- 1

# Set the number of iterations
num_iter <- 10000

# Initialize the posterior samples
alpha_samples <- rep(0, num_iter)
beta_samples <- rep(0, num_iter)
mu_samples <- matrix(0, num_iter, N)
x_estimate_samples <- matrix(0, num_iter, N)
y_estimate_samples <- matrix(0, num_iter, N)

# Run the Gibbs sampler
for (i in 1:num_iter) {
  # Sample x_estimate
  for (j in 1:N) {
    var_x_estimate <- 1 / ((1 / (x_sd[j]^2)) + (1 / (tau^2)))
    mean_x_estimate <- var_x_estimate * ((x_measurement[j] / (x_sd[j]^2)) + (mu_samples[i,j] / (tau^2)))
    x_estimate_samples[i,j] <- rnorm(1, mean_x_estimate, sqrt(var_x_estimate))
  }

  # Sample y_estimate
  for (j in 1:N) {
    var_y_estimate <- 1 / ((1 / (y_sd[j]^2)) + (1 / (tau^2)))
    mean_y_estimate <- var_y_estimate * ((y_measurement[j] / (y_sd[j]^2)) + (mu_samples[i,j] / (tau^2)))
    y_estimate_samples[i,j] <- rnorm(1, mean_y_estimate, sqrt(var_y_estimate))
  }

  # Sample mu
  for (j in 1:N) {
    var_mu <- 1 / ((1 / (tau^2)) + ((1 / (beta^2)) * (1 / (x_sd[j]^2))))
    mean_mu <- var_mu * (((y_estimate_samples[i,j] / (y_sd[j]^2)) * beta) + alpha + ((x_estimate_samples[i,j] / (x_sd[j]^2)) * beta))
    mu_samples[i,j] <- rnorm(1, mean_mu, sqrt(var_mu))
  }

  # Sample alpha
  var_alpha <- 1 / ((1 / (tau^2)) + ((1 / (beta^2)) * sum(1 / (x_sd^2))))
  mean_alpha <- var_alpha * sum((mu_samples[i,] - (beta * x_estimate_samples[i,])) / (tau^2))
  alpha_samples[i] <- rnorm(1, mean_alpha, sqrt(var_alpha))

  # Sample beta
  var_beta <- 1 / ((1 / (tau^2)) + ((1 / (beta^2)) * sum((x_estimate_samples[i,]^2) / (x_sd^2))))
  mean_beta <- var_beta * sum((x_estimate_samples[i,] * (mu_samples[i,] - alpha_samples[i])) / (tau^2))
  beta_samples[i] <- rnorm(1, mean_beta, sqrt(var_beta))
}

