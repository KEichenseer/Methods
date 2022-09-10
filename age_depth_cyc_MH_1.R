### Age Depth model with cycles and dates
## 07/09/2022
## K Eichenseer

##
## Create test data
##
cycle_pos <- cumsum(c(0,50,20,20,20,20,50,30,30,30,30,30,30,30,50,50,50,20,20,
                      20,20,20,20,20,30,30,30,40,40,40,40,40,40))
delta_known <- 0.5
alpha_known <- 540

posd1 <- 3
posd2 <- 10
posd3 <-20

dates <- data.frame(
  height = c(cycle_pos[posd1]+0.5*diff(cycle_pos)[posd1],
             cycle_pos[posd2]+0.2*diff(cycle_pos)[posd2],
             cycle_pos[posd3]+0*diff(cycle_pos)[posd3]),
  mean = c(alpha_known-(posd1-0.5)*delta_known,
           alpha_known-(posd2-0.5)*delta_known,
           alpha_known-(posd3-0.5)*delta_known),
  sd = c(.5,.5,.75))

points_of_interest <- 400


##
## Plot test data
##
plot(0,0,type = "n", xlim = range(cycle_pos),
     ylim = c(alpha_known+2,alpha_known-delta_known*length(cycle_pos)-1),
     xlab = "height (m)", ylab = "age (Ma)")
points(dates$height,dates$mean, lwd = 2, cex = 1.3)
sapply(seq_len(nrow(dates)), function(x) points(rep(dates$height[x],2),
       dates$mean[x]+dates$sd[x]*c(-2,2), type = "l", lwd = 2))
#points(points_of_interest,rep(500,length(points_of_interest)),
#       pch = 8, col = "red", lwd = 1)
# cycles
for(i in 1:(length(cycle_pos)-1))
  points(seq(cycle_pos[i],cycle_pos[i+1],0.02),alpha_known+1+1*sin(
    (0.25*diff(cycle_pos)[i]+seq(0,diff(cycle_pos)[i],0.02))
    *2*pi/diff(cycle_pos)[i]),
    type = "l", col = "black")
points(cycle_pos,rep(542,length(cycle_pos)),pch = 3, cex = 0.8)

##
## define age model
##
age_model <- function(n_iter,dates,cycle_pos,
                      alpha_init = 540,
                      delta_init = 2,
                      alpha_step = 0.1,
                      delta_step = 0.025) {


  # calculate position of dates relative to cycles (height scale)
  # 0 would be at start of the first cycle
  n_dates <- nrow(dates)
  date_cycle_pos <- rep(NA_real_,n_dates)
  for(d in 1:n_dates) {
    c_ind <- max(which(cycle_pos<=dates$height[d]))
    date_cycle_pos[d] <-
    c_ind-1+(dates$height[d] - cycle_pos[c_ind]) / diff(cycle_pos)[c_ind]
  }

  # store alpha and delta
  alpha <- rep(NA_real_,n_iter)
  delta <- rep(NA_real_,n_iter)

  # intial values
  alpha[1] = alpha_init
  delta[1] = delta_init

  date_age_new <- rep(NA_real_,n_dates)
  for(d in seq_len(n_dates)) {
    date_age_new[d] <- alpha[1] - delta[1]*date_cycle_pos[d]
  }

  loglik_old <- sum(dnorm(date_age_new, dates$mean, dates$sd, log = T))

  for (i in 2:n_iter) {
    # propose new alphas and deltas
    alpha_new = rnorm(1,alpha[i-1],alpha_step)
    delta_new = rnorm(1,delta[i-1],delta_step)

    # calculate new date ages
    for(d in seq_len(n_dates)) {
      date_age_new[d] <- alpha_new - delta_new*date_cycle_pos[d]
    }

    # likelihood
    loglik_new <- sum(dnorm(date_age_new, dates$mean, dates$sd, log = T))

    # decide acceptance
    HR <- exp(loglik_new-loglik_old)
    accept <- runif(1)<HR
    if(accept) {
      alpha[i] <- alpha_new
      delta[i] <- delta_new
      loglik_old <- loglik_new
    } else {
      alpha[i] <- alpha[i-1]
      delta[i] <- delta[i-1]
    }
  }
  data.frame(alpha = alpha,
             delta = delta)
}

##
## test age model
##
n_iter = 25000
m1 <- age_model(n_iter = n_iter,dates = dates, cycle_pos = cycle_pos)

##
## Assess model output
##
assess_posterior = FALSE
if(assess_posterior) {
# assess chains
  par(mfrow = c(2,1), mar = c(4.25,4.25,1,1))
plot(m1$alpha)
plot(m1$delta)
# cross plot
par(mfrow = c(1,1))
plot(m1)
}

##
## function to calculate age of points of interest from the posterior
##
calc_age <- function(heights, cycle_pos, alpha, delta) {
  # calculate position of heights relative to cycles (height scale)
  # 0 would be at start of the first cycle
  heights_cycle_pos <- rep(NA_real_,length(heights))
  for(d in 1:length(heights)) {
    c_ind <- max(which(cycle_pos<=heights[d]))
    heights_cycle_pos[d] <-
      c_ind-1+(heights[d] - cycle_pos[c_ind]) / diff(cycle_pos)[c_ind]
    if(heights[d]==max(cycle_pos)) heights_cycle_pos[d] <- length(cycle_pos)-1

  }

  out <- matrix(nrow = length(alpha), ncol = length(heights))
  for(d in seq_len(length(heights))) {
    out[,d] <- alpha - delta*heights_cycle_pos[d]
  }

  # return output
  out
}

# add posterior age estimate of point of interest
burn_in <- 0.1*n_iter

pos1 <- calc_age(heights = points_of_interest,cycle_pos = cycle_pos, alpha = m1$alpha[burn_in:n_iter],
                 delta = m1$delta[burn_in:n_iter])
points(points_of_interest,apply(pos1,2,mean), col = "red", lwd = 2, pch = 1, cex = 1.3)
sapply(seq_len(length(points_of_interest)), function(x) points(
  rep(points_of_interest[x],2),
  apply(pos1,2,function(x) quantile(x, probs = c(0.025,0.975))),
  type = "l",lwd = 2, col = "red"))

# add overall posterior age - height line
points_of_interest <- seq(0,max(cycle_pos),1)
pos1 <- calc_age(heights = points_of_interest,cycle_pos = cycle_pos, alpha = m1$alpha[burn_in:n_iter],
                 delta = m1$delta[burn_in:n_iter])

points(points_of_interest,apply(pos1,2,mean), col = "red", type = "l")
points(points_of_interest,apply(pos1,2,function(x) quantile(x,probs=0.025)), type = "l", lty = 3)
points(points_of_interest,apply(pos1,2,function(x) quantile(x,probs=0.975)), type = "l", lty = 3)

points(dates$height,dates$mean, col = "dodgerblue", lwd = 2, pch = 1, cex = 1.3)
sapply(seq_len(nrow(dates)), function(x) points(
  rep(dates$height[x],2),dates$mean[x]+dates$sd[x]*c(-2,2),
  type = "l",col = "dodgerblue", lwd = 2))

legend("topleft", legend =
         c("dates", "point of interest", "cycle positions", "age - height curve"),
       pch = c(1,1,3,NA), lty = c(NA,NA,1,1), lwd = c(2,2,1,1),
       col = c("dodgerblue", "red", "black", "red"))
