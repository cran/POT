library(POT)
set.seed(123)

mc <- simmc(100, alpha = 0.25)
mc <- qgpd(mc, 0, 1, 0.25)
##A first application when marginal parameter are estimated
f1 <- fitmcgpd(mc, 0)
str(f1)
f1
