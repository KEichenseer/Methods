dat <- climate_proxy("D18O")
dat2 <- climate_bins(c("Maastrichtian"))
head(dat)
class(dat$temperature)
library(dplyr)
test <- dat %>% filter(stage=="Zanclean")
plot(test$palaeolat,test$temperature)
uniloc <- unique(cbind(test$palaeolat,test$palaeolng))
test <- test %>% mutate(location = NA)
for(i in 1:nrow(uniloc)) test$location[which(test$palaeolat==uniloc[i,1] & test$palaeolng==uniloc[i,2])] <- i
table(test$location)

y <- lapply(1:nrow(uniloc), function(x) test$temperature[which(test$location==x)])
y[1:3]
x <- abs(uniloc[,1])



plot(x,sapply(y,mean))



npb <- sapply(y,length)
y[which(npb==0)] <- NA

### remove locations with 1 obs

y <- y[npb!=1]
x <- x[npb!=1]

npb <- npb[npb!=1]

nbin = length(npb)
prop_sd_yest <- matrix(0.05,nrow = nbin, ncol = nbin)
diag(prop_sd_yest) <- 1
prop_sd_coeff <- c(.5,.5,.5,0.05)

coeff_inits = c(10,30,45,0.01)
yest_inits = rep(25,nbin) #c(30,25,20,15,10,0,0)
sdyest_inits = rep(2,nbin)

nIter = 10000

system.time({zanc1 <-  run_MCMC(nIter = nIter, x = x, yobs = y,
                 coeff_inits = coeff_inits,
                 sdy_init = 1, yest_inits = yest_inits,
                 sdyest_inits = sdyest_inits,
                 prop_sd_coeff=prop_sd_coeff, prop_sd_yest=prop_sd_yest,
                 nbin = nbin)})
burnin = 5000+1
mn <- zanc1
latitude <- 0:90
sample_it <- sample((burnin+1):nIter,2000)
grad_025 <- sapply(1:length(latitude), function(f) quantile(apply(mn[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.025))
grad_975 <- sapply(1:length(latitude), function(f) quantile(apply(mn[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.975))


plot(seq(0,90,0.1), gradient(seq(0,90,0.1), apply(mn[[1]][burnin:nIter,1:4],2,median), 0),
     ylim = c(-1,35), type = "l", lwd = 3, ylab = "Temperature", xlab = "Latitude")
### Calculate confidence intervals

error_polygon(latitude,grad_025,grad_975,rgb(0,0,0,0.15))

for(i in 1:nbin) if(!(is.na(prior_df[i,1])) & prior_df[i,1]=="dunif")  lines(x[c(i,i)],c(prior_df[i,2],prior_df[i,3]), lwd = 7, col = rgb(0,0.9,0,0.5))
for(i in 1:nbin) {if(!(is.na(prior_df[i,1])) & prior_df[i,1]=="dnorm")  {beanplot::beanplot(rnorm(1000,prior_df[i,2],prior_df[i,3]), add = T,
                                                                                            at = x[i], maxwidth = 5, side = "second",
                                                                                            what = c(0,1,0,0), col = rgb(0,0.7,0.7,0.5), border = NA)}}

for(i in 1:5) points(latitude,gradient(latitude,unlist(mn[[1]][sample(burnin:nIter,1),1:4]),0), type = "l", lwd = 2, col = rgb(0,0,0,0.3))

#sapply(1:nbin, function(a) points(c(x[a],x[a]),c(median(mn[[2]][burnin:nIter,a])+c(-2,2)*median(mn[[3]][burnin:nIter,a])), type = "l", col= "red"))

loccol <- colorRampPalette(c(rgb(1,1,0,0.5),rgb(0,1,0,0.5),rgb(0,1,1,0.5),rgb(0,0,1,0.5)), alpha = TRUE)(22)

for(i in (1:nbin)) if(npb[i] !=0) points(rep(x[i],npb[i])+rnorm(npb[i],0,0.25),y[[i]],
                                         col = loccol[i], pch = 17)
points(x,apply(mn[[2]][burnin:nIter,],2,median), bg= loccol, lwd =1, col = rgb(0,0,0,0.5), pch = 21, cex = 1.45)

sapply(1:nbin, function(a) points(c(x[a],x[a]),quantile(mn[[2]][burnin:nIter,a], probs = c(0.05,0.95)),
                  type = "l", col= loccol[a], lwd =2))
points(x,sapply(y,mean), bg = loccol, pch = 23, cex = 1, col = rgb(0,0,0,0.5))

legend("topright",legend=c("regression line", "temperature estimate","dO18 data", "coral reef range (uniform prior)", "additional info. (normal prior)"),
       cex = 0.8, col = c("black", "red", rgb(0,0.3,1,0.55),  rgb(0,0.9,0,0.5), rgb(0,0.7,0.7,0.5)), lwd = c(2,NA,NA,4,4),
       pch = c(NA,19,17,NA,NA), pt.cex = c(NA,1.25,1,NA,NA))

plot(mn[[1]][,1],mn[[1]][,2])

plot(mn[[1]][,4])

#### exploring the gamma distribution

a <- c(0.75,1,1.5,2)
b <- c(0.)
n <- c(1,2,3,10)
plot(sapply(y,length),sapply(y,var))
hist(1/rgamma(1000,a[1]+n[1]/2, b+0.5*(yobs[[j]]-yestimate[i-1,j])^2))
