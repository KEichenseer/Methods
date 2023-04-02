
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


  y_est =
    rnorm(n_obs,
          sigma2/(y_sd^2+sigma2)*y_obs + y_sd^2/(y_sd^2+sigma2)*y_pred,
          sqrt(1/(1/sigma2 + 1/y_sd^2))) # I think that one at least should be correct now... If x_sd is
  #  set to very small values, then the results mostly agree with those of the jags implementation


  if(beta[2] == 0) x_pred = rep(beta[1],length(x_obs)) else x_pred = (y_pred-beta[1])/beta[2]


  X_est[,2] =
  #   rnorm(n_obs,
  #                   ((beta[2] * x_obs / x_sd^2) + ((y_est - beta[1] - y_pred) * beta[2]/ sigma2)) / ((beta[2]^2 / x_sd^2) + (1 / sigma2)),
  #                   1 / ((beta[1]^2 / x_sd^2) + (1 / sigma2))
  # )
    rnorm(n_obs,
          (beta[2]*(y_est - beta[1]) / sigma2 + x_obs / x_sd^2) / (1 / x_sd^2 + beta[2]^2 / sigma2), #(x_obs / x_sd^2 + beta[2]*y_est / sigma2) / (1 / x_sd^2 + beta[2]^2 / sigma2),
          sqrt(1 / (1 / x_sd^2 + beta[2]^2 / sigma2))
    )

  # rnorm(n_obs,
  #       (beta[2]*y_est + beta[1]*x_sd^2/x_obs) / (beta[2]^2*x_sd^2/x_obs + 1),
  #       (sigma2*x_sd^2) / (beta[2]^2*x_sd^2/x_obs + 1)
  # )
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


points(apply(x_est_out,2,mean),apply(y_est_out,2,mean),pch = 21, col = NA, bg = rgb(0,0.9,0,0.33))

points(apply(fit_sd$BUGSoutput$sims.list$x_estimate,2,mean),
       apply(fit_sd$BUGSoutput$sims.list$y_estimate,2,mean),
       pch = 24, col = NA, bg = rgb(0.9,0,0.5,0.33))


