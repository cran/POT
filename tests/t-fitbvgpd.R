library(POT)
set.seed(123)

x <- rbvgpd(50, model = "log", alpha = 0.5, mar1 = c(0, 1, 0.2))
f1 <- fitbvgpd(x, c(1,1)/2)
convassess(f1)
str(f1)
f1
