##A demo for bivariate POT
x <- rbvgpd(100, "alog", alpha = 0.35, asCoef1 = 0.75,
            asCoef2 = 0.25, mar1 = c(0, 1, 0.25),
            mar2 = c(2, 3, -0.15))
Mlog <- fitbvgpd(x, c(0, 2), "log")
Malog <- fitbvgpd(x, c(0, 2), "alog")
Malog1 <- fitbvgpd(x, c(0,2), "alog", asCoef2 = 0.25)
anova(Malog, Mlog)
anova(Malog1, Mlog)
anova(Malog, Malog1)

layout(matrix(c(1,1,2,2,0,3,3,0), 2, byrow = TRUE))
plot(Mlog)
x11()
layout(matrix(c(1,1,2,2,0,3,3,0), 2, byrow = TRUE))
plot(Malog)

readline("Press ENTER to continue")
dev.off()
dev.off()
##Simulate a first order Markov Chain with extreme value dependence
##and GPD observations
mc <- simmc(32000, alpha = 0.3, asCoef1 = 0.3, asCoef2 = 0.1)
mc <- qgpd(mc, loc = 0, scale = 2, 0.25)
plot(mc, type = "l")
layout(matrix(c(1,1,2,2,0,3,3,0), 2, byrow = TRUE))
mrlplot(mc, u.range = c(16, 22), nt = 2000)
abline(v = 19.5, col = "blue")
tcplot(mc, u.range = c(16, 22), which = 1)
abline(v = 19.5, col = "blue")
tcplot(mc, u.range = c(16, 22), which = 2)
abline(v = 19.5, col = "blue")

Mcalog <- fitmcgpd(mc, 19.5, "alog")
par(mfrow=c(2,2))
plot(Mcalog, npy = 1)
