library(tsdl)
library(MASS)
library(ggplot2)
library(ggfortify)
library(forecast)
library(MuMIn)
library(astsa)

# get data
xsales = tsdl[[358]]

# data information
length(xsales)
attr(xsales, "description")
attr(xsales, "source")

# set to ts
xts = ts(xsales, start = c(1965,1), end = c(1971,5), frequency = 12)

# view original plot
ts.plot(xts, main = "Raw Data")

# create training and test sets
xtrain = xts[c(1:72)]
xtest = xts[c(73:77)]

# mean and variance of xtrain
mean(xtrain)
var(xtrain)

# view plot of training set
ts.plot(xtrain, main = "Training Data")
xfit = lm(xtrain ~ as.numeric(1:length(xtrain)));abline(xfit, col = 'red')
abline(h = mean(xtrain), col = 'blue')

# view hist and acf of training data
hist(xtrain, col="light blue", main="Histogram of Company X's Sales")
acf(xtrain, lag.max = 40, main = "ACF of Company X's Sales")




# box cox transform
bcTransform = boxcox(xtrain ~ as.numeric(1:length(xtrain)), plotit = FALSE)
lambda = bcTransform$x[which(bcTransform$y == max(bcTransform$y))]
lambda
xtrain.bc = (1/lambda)*((xtrain^lambda)-1)

# log transform
xtrain.log = log(xtrain)

# sqrt transform
xtrain.sqrt = sqrt(xtrain)

# plot all transformations
par(mfrow = c(2,2))
ts.plot(xtrain, main = "Original")
ts.plot(xtrain.bc, main = "Box Cox")
ts.plot(xtrain.log, main = "Log")
ts.plot(xtrain.sqrt, main = "Square Root")

# plot all histograms
hist(xtrain, main = "Original")
hist(xtrain.bc, main = "Box Cox")
hist(xtrain.log, main = "Log")
hist(xtrain.sqrt, main = "Square Root")

# final transformed data
xtrans = xtrain.log

decomp = decompose(ts(as.ts(xtrans), frequency = 12))
plot(decomp)

# differencing lag 12 to remove seasonality
par(mfrow = c(1,1))
xd12 = diff(xtrans, 12)
ts.plot(xd12, main = "ln(X_t) diff at lag 12")
xd12.fit = lm(xd12~as.numeric(1:length(xd12)));abline(xd12.fit, col = "red")
abline(h=mean(xd12), col = "blue")
var(xd12)
mean(xd12)

# differencing at lag 1 to remove trend
xstat = diff(xd12, 1)
ts.plot(xstat, main = "ln(X_t) diff at lag 12 and lag 1")
xstat.fit = lm(xstat~as.numeric(1:length(xstat)));abline(xstat.fit, col = "red")
abline(h=mean(xstat), col = "blue")
var(xstat)
mean(xstat)

# acfs of different data
par(mfrow=c(2,2))
acf(xtrain, lag.max = 40, main = "ACF of X_t")
acf(xtrans, lag.max = 40, main = "ACF of ln(X_t)")
acf(xd12, lag.max = 40, main = "ACF of ln(X_t) diff at lag 12")
acf(xstat, lag.max = 40, main = "ACF of ln(X_t) diff at lag 12 and 1")

hist(xtrain, breaks = 16)
hist(xtrans, breaks = 16)
hist(xd12, breaks = 16)
hist(xstat, breaks = 16)

hist(xstat, breaks = 16, col = 'light blue', prob = TRUE)
curve(dnorm(x, mean(xstat), sqrt(var(xstat))), add = TRUE)




acf(xstat, lag.max = 40)
pacf(xstat, lag.max = 40)

# possible models
# SARIMA(2,1,0)x(0,1,1)
# SARIMA(2,1,0)x(0,1,2) no, principle of parsimony
# SARIMA(2,1,0)x(1,1,1) lower aic

A = arima(xtrans, order=c(2,1,0), seasonal = list(order = c(0,1,1), period = 12), method="ML")
B = arima(xtrans, order=c(2,1,0), seasonal = list(order = c(0,1,2), period = 12), method="ML")
C = arima(xtrans, order=c(2,1,0), seasonal = list(order = c(1,1,1), period = 12), method="ML")
#arima(xtrans, order=c(0,1,12), method = "ML")

#for (i in c(1,11)){
#  for (j in 0:2){
#    for (k in c(1)){
#      for (l in c(1)){
#        print(i); print(j); print(k); print(l); print(AICc(arima(xtrans,
#                                                                 order = c(i,0,j),
#                                                                 seasonal = list(order = #c(k,1,l),period = 12),method = "ML")))}}}}

# final model
arima(xtrans, order=c(2,1,0), seasonal = list(order = c(0,1,1), period = 12), method="ML")




# plot roots function
plot.roots <- function(ar.roots=NULL, ma.roots=NULL, size=2, angles=FALSE, special=NULL, sqecial=NULL,my.pch=1,first.col="blue",second.col="red",main=NULL)
{xylims <- c(-size,size)
omegas <- seq(0,2*pi,pi/500)
temp <- exp(complex(real=rep(0,length(omegas)),imag=omegas))
plot(Re(temp),Im(temp),typ="l",xlab="x",ylab="y",xlim=xylims,ylim=xylims,main=main)
abline(v=0,lty="dotted")
abline(h=0,lty="dotted")
if(!is.null(ar.roots))
{
  points(Re(1/ar.roots),Im(1/ar.roots),col=first.col,pch=my.pch)
  points(Re(ar.roots),Im(ar.roots),col=second.col,pch=my.pch)
}
if(!is.null(ma.roots))
{
  points(Re(1/ma.roots),Im(1/ma.roots),pch="*",cex=1.5,col=first.col)
  points(Re(ma.roots),Im(ma.roots),pch="*",cex=1.5,col=second.col)
}
if(angles)
{
  if(!is.null(ar.roots))
  {
    abline(a=0,b=Im(ar.roots[1])/Re(ar.roots[1]),lty="dotted")
    abline(a=0,b=Im(ar.roots[2])/Re(ar.roots[2]),lty="dotted")
  }
  if(!is.null(ma.roots))
  {
    sapply(1:length(ma.roots), function(j) abline(a=0,b=Im(ma.roots[j])/Re(ma.roots[j]),lty="dotted"))
  }
}
if(!is.null(special))
{
  lines(Re(special),Im(special),lwd=2)
}
if(!is.null(sqecial))
{
  lines(Re(sqecial),Im(sqecial),lwd=2)
}
}

# plot roots, check stationarity of AR

arcoef = polyroot(c(1, -0.5628, -0.0334))
macoef = polyroot(c(1, -0.5322))
plot.roots(arcoef, macoef)




# set model
arima.fit = A

# get residuals
res = residuals(arima.fit)
mean(res)
var(res)

# hist of residuals
hist(res, breaks = 20, col = "blue", prob = TRUE)
curve(dnorm(x, mean(res), sqrt(var(res))), add=TRUE)

# plot of residuals
ts.plot(res)
res.fit = lm(res~as.numeric(1:length(res)));abline(res.fit, col = 'red')
abline(h=mean(res), col='blue')

# qq plot
qqnorm(res)
qqline(res,col="blue")

# acf and pacf of residuals
acf(res^2, lag.max = 40)
pacf(res, lag.max = 40)

# shapiro test !!!!
shapiro.test(res)

# box-pierce test
Box.test(res, lag = 8, type=c("Box-Pierce"), fitdf = 3)

# ljung-box test
Box.test(res, lag = 8, type=c("Ljung-Box"), fitdf = 3)

# mcleod-li test !!!!
Box.test((res)^2, lag = 8, type=c("Ljung-Box"), fitdf = 0)

# yule-walker
ar(res, aic = TRUE, order.max = NULL, method = c("yule-walker"))

# VIEW ALL
fit.i <- sarima(xdata=xtrans, p=2, d=1, q=0, P=0, D=1, Q=1, S=12)




# forecast future values
fore.fit = A
forecast(fore.fit)

# graph with 5 future values
pred.A = predict(fore.fit, n.ahead = 5)
up.A = pred.A$pred + 2*pred.A$se
low.A = pred.A$pred - 2*pred.A$se

ts.plot(xtrans, xlim = c(1, length(xtrans)+5), ylim = c(min(xtrans),max(up.A)))
lines(up.A, col = "blue", lty = "dashed")
lines(low.A, col = "blue", lty = "dashed")
points((length(xtrans)+1):(length(xtrans)+5),pred.A$pred, col='red')

# original data with predictions
pred.orig = exp(pred.A$pred)
up.orig = exp(up.A)
low.orig = exp(low.A)

ts.plot(xtrain, xlim = c(1, length(xtrain)+5), ylim = c(min(xtrain),max(up.orig)))
lines(up.orig, col = 'blue', lty = 'dashed')
lines(low.orig, col = 'blue', lty = 'dashed')
points((length(xtrain)+1):(length(xtrain)+5), pred.orig, col = 'red')

xnum = xts[0:77]
ts.plot(xnum, xlim = c(67,77))
lines(up.orig, col = 'blue', lty = 'dashed')
lines(low.orig, col = 'blue', lty = 'dashed')
points((length(xtrain)+1):(length(xtrain)+5), pred.orig, col="red")
