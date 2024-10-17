require(POT)
?clust
data(ardieres)
par(mfrow=c(1,2))
clust(ardieres, 4, 10 / 365, plot=TRUE)
