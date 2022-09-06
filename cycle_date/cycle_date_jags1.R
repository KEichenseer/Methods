regression_model_with_errors <- function() {

  ## Likelihood
  # dates
  for (i in 1:n_date){
    age_date[i] ~ dnorm(mean_date[i], 1/(sd_date[i]*sd_date[i]))
    age_date[i] <- alpha[i] + beta * age_date[i]
  }

  # cycles
  for (i in 1:n_date){
    age_date[i] ~ dnorm(mean_date[i], 1/(sd_date[i]*sd_date[i]))
    age_date[i] <- alpha[i] + beta * age_date[i]
  }


  ## Priors
  sigma ~ dnorm(0, 1/(100*100)); T(0,)
  tau <- 1 / (sigma * sigma) # precision = 1/standard deviation ^2
  alpha ~ dnorm(0, 1/(100*100))
  beta ~ dnorm(0, 1/(100*100))

}
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
