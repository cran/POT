library(POT)
set.seed(123)

x <- rgpd(30, 0, 1, 0.2)
f1 <- fitgpd(x, 0, "med")
convassess(f1)
str(f1)
f1
