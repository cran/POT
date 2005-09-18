## This is the demo file of package POT

data(ardieres)
plot(ardieres, ylab = expression(paste("Flood Discharge  ",m^3/s, sep="")),
     xlab = "Years", main = "River Ardières at Beaujeu")
cat("The threshold was initially to low. Select another threshold more suitable
for the GP assumtion...")
flows <- ardieres[,"flows"]
date <- ardieres[,"date"]
par(mfrow=c(1,2),ask = TRUE)
diplot(ardieres)
mrlplot(flows)
cat("A threshold around 6 should be raisonable...")
par(mfrow=c(2,2),ask = TRUE)
diplot(ardieres)
abline(v = 6)
mrlplot(flows)
abline(v = 6)
tcplot(flows, tlim =c(4,12))
cat("But, is the GP distribution suited to model our data")
par(mfrow=c(1,1))
lmomplot(flows, identify = FALSE)
cat("Now, fit the GP distribution to exceedances above this threshold.")
fitted <- fitgpd(flows, 6, 'mle',corr = TRUE)
readline("Press ENTER to continue")
cat("If, we want PWM estimates - unbiased one.")
fitgpd(flows, 6, 'pwmu')
readline("Press ENTER to continue")
cat("Try also, `fitgpd(flows, 6, 'pwmb')' and `fitgpd(flows, 6, 'moments')'
We can compute profile likelihood confidence interval...")
par(mfrow=c(1,2))
gpd.pfshape(fitted, range = c(-0.1, 0.8))
gpd.pfscale(fitted, range = c(2, 7))
cat("Of course, classical confidence interval can be computed...")
gpd.fiscale(fitted)
readline("Press ENTER to continue")
gpd.fishape(fitted)
readline("Press ENTER to continue")
par(mfrow=c(1,1),ask=TRUE)
cat("The same can also be obtained for return levels...")
mu <-  fitted$nhigh / diff(range(date))
gpd.pfrl(fitted, 10, mu, range=c(14, 35),
         main=expression(paste("95% Profile C.I. for the 10-year return level, "
             ,mu==1.5,sep="")))
cat("We can produce several graphics really usefull for diagnostic of our model
of Peaks over threshold...")
plot(ardieres, ylab = expression(paste("Flood Discharge  ",m^3/s, sep="")),
     xlab = "Years", main = "River Ardières at Beaujeu")
abline( h = 6, lty = 2)
par(mfrow=c(2,2))
plotgpd(fitted, mu)
