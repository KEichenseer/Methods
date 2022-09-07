scale_int <- function(x, int, na.rm = T){
  (x-min(x, na.rm = na.rm))/(max(x, na.rm = na.rm)-min(x, na.rm = na.rm)) * 
    (max(int, na.rm = na.rm)-min(int, na.rm = na.rm)) + min(int, na.rm = na.rm)}

dh = 50
dc = 2

x <- seq(0,dh,0.1)
exponent <- exp(rnorm(1,0,0.25))
xmod <- scale_int(seq(0,1,length.out=nx)^exponent,c(0,dh))

age <- dc*xmod/dh
points(x,age, type = "l", col = rgb(0,0.33,1,0.5))


plot(x,age, type = "l")
