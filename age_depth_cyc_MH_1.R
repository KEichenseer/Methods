### Age Depth model with cycles and dates



alpha <- rep(NA_real_,n_iter)
delta <- rep(NA_real_,n_iter)

# intial values
alpha[1] = alpha_init
delta[1] = delta_init

for(d in seq_len(n_dates)) {
  date_age_new[d] <- alpha_new + delta_new*date_cycle_pos[d]
}


for (i in 2:length(n_iter)) {
  alpha_new = rnorm(1,alpha[i],alpha_step)
  delta_new = rnorm(1,delta[i],delta_step)

  for(d in seq_len(n_dates)) {
    date_age_new[d] <- alpha_new + delta_new*date_cycle_pos[d]
  }

  loglik_new <- dnorm(date_age_new, date_age_mean, date_age_sd)

  HR <- exp(loglik_new-loglik_old)
  accept <- runif(1)<HR
  if(accept) {
    alpha[i] <- alpha_new
    delta[i] <- delta_new
    loglik_old <- loglik_new
  }
}

