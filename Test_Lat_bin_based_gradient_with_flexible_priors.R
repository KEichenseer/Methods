dat <- climate_proxy("D18O")
dat2 <- climate_bins(c("Maastrichtian", "Hettangian"))
head(dat)
class(dat$temperature)
library(dplyr)
test <- dat %>% filter(stage=="Maastrichtian")
plot(test$palaeolat,test$temperature)
uniloc <- unique(cbind(test$palaeolat,test$palaeolng))
test <- test %>% mutate(location = NA)
for(i in 1:nrow(uniloc)) test$location[which(test$palaeolat==uniloc[i,1] & test$palaeolng==uniloc[i,2])] <- i
table(test$location)

y <- lapply(1:nrow(uniloc), function(x) test$temperature[which(test$location==x)])
y[1:3]
x <- abs(uniloc[,1])

plot(x,sapply(y,mean))

nbin = nrow(uniloc)
prop_sd_yest <- matrix(0.001,nrow = nbin, ncol = nbin)
diag(prop_sd_yest) <- 0.1
prop_sd_coeff <- c(.1,.1,.1,0.001)

coeff_inits = c(10,30,45,0.01)
yest_inits = rep(25,nbin) #c(30,25,20,15,10,0,0)
sdyest_inits = rep(2,nbin)

npb <- sapply(y,length)
y[which(npb==0)] <- NA

system.time({maas3 <-  run_MCMC(nIter = 500000, x = x, yobs = y, prior_df = NA,
                 coeff_inits = coeff_inits,
                 sdy_init = 1, yest_inits = yest_inits,
                 sdyest_inits = sdyest_inits,
                 prop_sd_coeff=prop_sd_coeff, prop_sd_yest=prop_sd_yest,
                 nbin = nbin)})
burnin = 100000+1
nIter = 500000
mn <- maas3
plot(seq(0,90,0.1), gradient(seq(0,90,0.1), apply(mn[[1]][burnin:nIter,1:4],2,median), 0),
     ylim = c(-1,35), type = "l", lwd = 3, ylab = "Temperature", xlab = "Latitude")
### Calculate confidence intervals
latitude <- 0:90
sample_it <- sample((burnin+1):nIter,2000)
grad_025 <- sapply(1:length(latitude), function(f) quantile(apply(mn[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.025))
grad_975 <- sapply(1:length(latitude), function(f) quantile(apply(mn[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.975))
error_polygon(latitude,grad_025,grad_975,rgb(0,0,0,0.15))

for(i in 1:nbin) if(!(is.na(prior_df[i,1])) & prior_df[i,1]=="dunif")  lines(x[c(i,i)],c(prior_df[i,2],prior_df[i,3]), lwd = 7, col = rgb(0,0.9,0,0.5))
for(i in 1:nbin) {if(!(is.na(prior_df[i,1])) & prior_df[i,1]=="dnorm")  {beanplot::beanplot(rnorm(1000,prior_df[i,2],prior_df[i,3]), add = T,
                                                                                            at = x[i], maxwidth = 5, side = "second",
                                                                                            what = c(0,1,0,0), col = rgb(0,0.7,0.7,0.5), border = NA)}}


#sapply(1:nbin, function(a) points(c(x[a],x[a]),c(median(mn[[2]][burnin:nIter,a])+c(-2,2)*median(mn[[3]][burnin:nIter,a])), type = "l", col= "red"))

for(i in (1:nbin)) if(npb[i] !=0) points(rep(x[i],npb[i])+rnorm(npb[i],0,0.25),y[[i]], col = rgb(0,0.3,1,0.55), pch = 17)
points(x,apply(mn[[2]][burnin:nIter,],2,median), col = rgb(1,0,0,0.67), pch = 19, cex = 1.25)
sapply(1:nbin, function(a) points(c(x[a],x[a]),quantile(mn[[2]][burnin:nIter,a], probs = c(0.025,0.975)), type = "l", col= rgb(1,0,0,0.75), lwd =2))

legend("topright",legend=c("regression line", "temperature estimate","dO18 data", "coral reef range (uniform prior)", "additional info. (normal prior)"),
       cex = 0.8, col = c("black", "red", rgb(0,0.3,1,0.55),  rgb(0,0.9,0,0.5), rgb(0,0.7,0.7,0.5)), lwd = c(2,NA,NA,4,4),
       pch = c(NA,19,17,NA,NA), pt.cex = c(NA,1.25,1,NA,NA))


#### exploring the gamma distribution

a <- c(0.75,1,1.5,2)
b <- c(0.)
n <- c(1,2,3,10)
plot(sapply(y,length),sapply(y,var))
hist(1/rgamma(1000,a[1]+n[1]/2, b+0.5*(yobs[[j]]-yestimate[i-1,j])^2))
