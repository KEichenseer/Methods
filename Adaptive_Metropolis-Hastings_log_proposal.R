######################################################
######################################################
### A Metropolis-Hastings algorithm for Latitudinal 
### Temperature Gradients
###
### by Kilian Eichenseer
### February 2022
###
######################################################
######################################################
### FUNCTIONS

gradient <- function(x, coeff, sdy) { # sigma is labelled "sdy"
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  nu = coeff[5]
  return(A + max(c(K-A,0))/((1+(nu*exp(Q*(x-M))))^(1/nu)) + rnorm(length(x),0,sdy))
}

loglik <- function(x, y,  coeff, sdy) { 
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  nu = coeff[5]
  pred = A + max(c(K-A,0))/((1+(nu*exp(Q*(x-M))))^(1/nu))
  return(sum(dnorm(y, mean = pred, sd = sdy, log = TRUE)))
}

logprior <- function(coeff) {
  return(sum(c(
    dunif(coeff[1], -4, 40, log = TRUE),
    dunif(coeff[2], -4, 40, log = TRUE),
    dnorm(coeff[3], 45, 10, log = TRUE),
    dgamma(coeff[4], 0.2, 0.2, log = TRUE),
    dgamma(coeff[5], 1, 1, log = TRUE))))
}

logposterior <- function(x, y, coeff, sdy){
  return (loglik(x, y, coeff, sdy) + logprior(log_t_Q_nu(coeff)))
}

MH_propose <- function(coeff, proposal_sd){
  return(rnorm(5,mean = coeff, sd= c(.5,.5,.5,0.01,0.07)))
}

weighted_var <- function(x, weights, sum_weights) {
  sum(weights*((x-sum(weights*x)/sum_weights)^2))/(sum_weights)
}

log_t_Q_nu <- function(x) c(x[1:3],log(x[4:5]))
exp_t_Q_nu <- function(x) c(x[1:3],exp(x[4:5]))
log_jacobian <- function(x) sum(-log(x))

# Main MCMCM function
run_MCMC <- function(x, y, coeff_inits, sdy_init, nIter, proposal_sd_init = rep(10,5), nAdapt = 10000){
  ### Initialisation
  coefficients = array(dim = c(nIter,5)) # set up array to store coefficients
  coefficients[1,] = coeff_inits # initialise coefficients
  sdy = rep(NA_real_,nIter) # set up vector to store sdy
  sdy[1] = sdy_init # intialise sdy
  A_sdy = 3 # parameter for the prior on the inverse gamma distribution of sdy
  B_sdy = 0.1 # parameter for the prior on the inverse gamma distribution of sdy
  n <- length(y)
  shape_sdy <- A_sdy+n/2 # shape parameter for the inverse gamma
  sd_it <- 1
  coeff_sd <- array(NA_real_,dim = c(nAdapt,5))
  coeff_sd[1:3,] <- proposal_sd_init
  coeff_diff <- array(NA_real_,dim = c(nAdapt,5))
  allWeights <- exp((-(nAdapt-2)):0/500)
  
  ### The MCMC loop
  for (i in 2:nIter){ 
    
    ## 1. Gibbs step to estimate sdy
    sdy[i] = sqrt(1/rgamma(1,shape_sdy,B_sdy+0.5*sum((y-gradient(x,coefficients[i-1,],0))^2)))
    
    ## 2. Metropolis-Hastings step to estimate the regression coefficients
    proposal =exp_t_Q_nu(MH_propose(log_t_Q_nu(coefficients[i-1,]),coeff_sd[sd_it,])) # new proposed values
    if(any(proposal[c(4,5)] <= 0)) HR = 0 else # Q and nu need to be >0
      # Hasting's ratio of the proposal
      HR = exp(logposterior(x = x, y = y, coeff = proposal, sdy = sdy[i])+log_jacobian(coefficients[i-1,4:5]) -
                 logposterior(x = x, y = y, coeff = coefficients[i-1,], sdy = sdy[i])-log_jacobian(proposal[4:5]))
    #print(i)
    #print(proposal)
    #print(HR)
    #print(coefficients[i-1,])
    
    # accept proposal with probability = min(HR,1)
    if (runif(1) < HR){ 
      coefficients[i,] = proposal
      # if proposal is rejected, keep the values from the previous iteration
    }else{
      coefficients[i,] = coefficients[i-1,]
    }

    # Adaptation of proposal SD
    if(i < nAdapt){
    
    coeff_diff[i,] <- exp_t_Q_nu(log_t_Q_nu(coefficients[i,])-log_t_Q_nu(coefficients[i-1,]))
    if(i>=3) {
    weights = allWeights[(nAdapt-i+2):nAdapt-1]
    sum_weights = sum(weights)
    weighted_var_coeff <- apply(coeff_diff[2:i,], 2, function(f) weighted_var(
                            f, weights = weights, sum_weights = sum_weights))
   
    
    coeff_sd[i+1,] <- ifelse(weighted_var_coeff==0,
                           sqrt(coeff_sd[i,]^2/10),
                           sqrt(2.4^2 * weighted_var_coeff))
    }
    sd_it <- i
    }
  } # end of the MCMC loop
  
  ###  Function output
  output = list(data.frame(A = coefficients[,1],
                      K = coefficients[,2],
                      M = coefficients[,3],
                      Q = coefficients[,4],
                      nu = coefficients[,5],
                      sdy = sdy),
                coeff_sd)
  return(output)
}
# for plotting the 95 % CI shading
error_polygon <- function(x,en,ep,color) {
  polygon( c(x[1], x, x[length(x)], x[length(x)], rev(x), x[1]),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
}


### Data creation
set.seed(10)
sample_lat <- runif(10,0,90)
sample_data <- data.frame(
  x = sample_lat, 
  y = gradient(x = sample_lat, coeff = c(-2.2, 28, 39, 0.1, 1.2), sd = 2))


### Analysis
nIter <- 100000
system.time({m <- run_MCMC(x = sample_data$x, y = sample_data$y, 
              coeff_inits = c(10,20,30,0.4,0.4), sdy_init = 5, nIter = nIter,
              nAdapt = 3000)})
matplot(log10(m[[2]]))

m[[2]][3000,]
######################################################
######################################################
### PLOTS

######################################################
## Plots of the data


layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))

par(mar = c(4,4.25,1.25,0.75), las = 1, mgp = c(2.25,0.75,0), cex = 1.35)
latitude <- seq(0,90,by=0.2)
temperature <- gradient(x = latitude, coeff = c(-2.2, 28, 39, 0.1, 1.2), sd = 0)
plot(latitude, temperature, type = "l", lwd = 3, ylim = c(-4.5,32), 
     xlab = expression ("absolute latitude ("*degree*")"), yaxt = "n", 
     ylab = expression("temperature ("*degree~"C)"), yaxs = "i", xaxs = "i",
     xlim = c(0,90), main = "gradient", cex.main = 1)
axis(2,seq(-5,30,5),c(NA,0,NA,10,NA,20,NA,30))


plot(latitude, temperature, type = "l", lwd = 2, ylim = c(-4.5,32), 
     xlab = expression ("absolute latitude ("*degree*")"), yaxt = "n", 
     ylab = expression("temperature ("*degree~"C)"), yaxs = "i", xaxs = "i",
     xlim = c(0,90), main = "sample data", lty= 2, cex.main = 1)
points(sample_data$x, sample_data$y,  pch = 19, cex = 1.2, col = rgb(0,0,0,0.6), xpd = T)
axis(2,seq(-5,30,5),c(NA,0,NA,10,NA,20,NA,30))

######################################################
## Posterior checks

library(ggplot2)
library(dplyr)
library(reshape)
library(cowplot)
m <- m[[1]] %>% mutate(iteration = 1:nrow(.))
mm <- melt(m,id.vars = "iteration")
mm <- mm %>% mutate(trueval = rep(c(-2.2, 28, 39, 0.1, 1.2, 2), each = nrow(m)))

facet.labs = c("F", "K", "M", "Q", expression(nu), expression(sigma["y"]))
names(facet.labs) <- c("A", "K", "M", "Q", "nu", "sdy")

p1 <- ggplot(mm %>% filter(iteration %% 10 == 0)) + geom_line(aes(x = iteration, y = value),
                                                              colour = rgb(0,0.35,0.7,1))+
  facet_grid(variable ~ ., switch = "y",scales = "free_y",
             labeller = labeller(variable = facet.labs)) + 
  theme(legend.position = "none")+
  theme_bw(base_size = 14)+
  theme(strip.text.y.left = element_text(angle = 0))+
  scale_y_continuous(position = "right")+
  ylab(NULL)+
  ggtitle(expression("traceplot (showing every 10th"~"iteration)"))+ 
  theme(plot.title = element_text(size = 13, face = "bold",hjust = 0.5),
        axis.title=element_text(size=13), axis.text.y=element_blank(),)


p2 <- ggplot(mm %>% filter(iteration %% 10 == 0)) +
  geom_vline(aes(xintercept = trueval), colour = "black", size = 1,linetype = "solid")+
  
  stat_density(aes(x = value, y =..scaled..), 
               fill = rgb(0,0.35,0.7,0.55), size = 1.5)+
  facet_grid(variable ~ ., scales = "free_y")+ 
  theme(legend.position = "none")+
  theme_bw(base_size = 14)+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())+ 
  scale_y_continuous(position = "left",breaks = c(0,0.5,1))+
  scale_x_continuous(position = "top")+
  coord_flip()+
  ylab("scaled density")+
  xlab(NULL)+
  ggtitle(expression("density"))+ 
  theme(plot.title = element_text(size = 13, face = "bold",hjust = 0.5),
        axis.title=element_text(size=13))

plot_grid(p1, p2, ncol = 2, labels = c("", ""),
          rel_widths = c(7,2))

######################################################
## Gradient from posterior
burnin <- 10000
layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))
par(mar = c(4,4.25,1.25,0.75), las = 1, mgp = c(2.25,0.75,0), cex = 1.35)
latitude <- seq(0,90,by=0.5)
temperature <- gradient(x = latitude, coeff = c(-2.2, 28, 39, 0.1, 1.2), sdy = 0)
plot(latitude, temperature, type = "l", lwd = 2, ylim = c(-4.5,32), 
     xlab = expression ("absolute latitude ("*degree*")"), yaxt = "n", 
     ylab = expression("temperature ("*degree~"C)"), yaxs = "i", xaxs = "i", lty = 2,
     xlim = c(0,90), main = "8 draws from the posterior", cex.main = 1)
points(sample_data$x, sample_data$y,  pch = 19, cex = 1.1, col = rgb(0,0,0,0.5), xpd = T)
replicate(8, points(latitude, gradient(x = latitude, coeff = 
                           unlist(m[sample((burnin+1):nIter,1),1:5]), sdy = 0),
                           type = "l", col = rgb(0,0.35,0.7,0.33), lwd = 2))
axis(2,seq(-5,30,5),c(NA,0,NA,10,NA,20,NA,30))

# Calculate confidence intervals
sample_it <- sample((burnin+1):nIter,1000)
grad_025 <- sapply(1:length(latitude), function(f) quantile(
  apply(m[sample_it,1:5],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), 
  probs = 0.025))
grad_975 <- sapply(1:length(latitude), function(f) quantile(
  apply(m[sample_it,1:5],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), 
  probs = 0.975))

# Calculate confidence intervals
sample_it <- sample((burnin+1):nIter,1000)
grad_025 <- sapply(1:length(latitude), function(f) quantile(
  apply(m[sample_it,1:5],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), 
  probs = 0.025))
grad_975 <- sapply(1:length(latitude), function(f) quantile(
  apply(m[sample_it,1:5],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), 
  probs = 0.975))

plot(latitude, temperature, type = "n", lwd = 2, ylim = c(-4.5,32), 
     xlab = expression ("absolute latitude ("*degree*")"), yaxt = "n", 
     ylab = expression("temperature ("*degree~"C)"), yaxs = "i", xaxs = "i", lty = 2,
     xlim = c(0,90), main = "posterior median and 95% CI", cex.main = 1)

error_polygon <- function(x,en,ep,color) {
  polygon( c(x[1], x, x[length(x)], x[length(x)], rev(x), x[1]),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
}
error_polygon(latitude,grad_025,grad_975,rgb(0,0.35,0.7,0.33))
points(sample_data$x, sample_data$y,  pch = 19, cex = 1.1, col = rgb(0,0,0,0.5), xpd = T)

points(latitude, gradient(x=latitude, coeff =  apply(m[sample_it,1:5],2,median), 
                          sdy = 0), type = "l", lwd = 3, col = rgb(0,0.35,0.7,0.75))
points(latitude, temperature, type = "l", lwd = 2, lty = 2)
axis(2,seq(-5,30,5),c(NA,0,NA,10,NA,20,NA,30))

