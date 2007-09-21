simmc <- function(n, alpha, model = "log", asCoef, asCoef1,
                  asCoef2, margins = "uniform"){

  if (missing(alpha))
    stop("``alpha must be present.")

  if ((model %in% c("alog", "anlog")) && (missing(asCoef1) ||
                                          missing(asCoef2)))
    stop("``asCoef1'' and ``asCoef2'' must be present.")

  if ((model == "amix") && missing(asCoef))
    stop("``asCoef'' must be present.")

  if ((model %in% c("log", "alog")) && ((alpha <= 0) ||
                                        (alpha > 1)))
    stop("``alpha'' must be in ]0, 1] with this model.")

  if ((model %in% c("nlog", "anlog")) && (alpha <= 0))
    stop("``alpha'' must be positive with this model.")

  if ((model == "mix") && ((alpha < 0) || (alpha > 1)))
    stop("``alpha'' must be in [0,1] with this model.")

  if ((model == "amix") &&
      ((alpha < 0) || (alpha + 2 * asCoef >1) ||
       (alpha + 3 * asCoef < 0)))
    stop("``alpha'' and ``asCoef'' are not valid. See the doc.")

  if ((model == "amix")){
    alpha <- alpha + 3 * asCoef
    asCoef <- - asCoef
  }

  alpha <- as.double(alpha)

  if (!missing(asCoef))
    asCoef <- as.double(asCoef)

  if (!missing(asCoef1) && !missing(asCoef2))
    asy <- as.double(c(asCoef1, asCoef2))
  
  evmc <- runif(n)
  nn <- as.integer(1)

  for (i in 2:n){
    evmc[c(i,i-1)] <-
      switch(model, log = .C("rbvlog", nn, alpha, sim = evmc[c(i,i-1)],
                      PACKAGE = "POT")$sim,
             alog = .C("rbvalog", nn, alpha, asy, sim = evmc[c(i,i-1)],
               PACKAGE = "POT")$sim,
             nlog = .C("rbvnlog", nn, alpha, sim = evmc[c(i,i-1)],
               PACKAGE = "POT")$sim,
             anlog = .C("rbvanlog", nn, alpha, asy, sim = evmc[c(i,i-1)],
               PACKAGE = "POT")$sim,
             mix = .C("rbvmix", nn, alpha, sim = evmc[c(i,i-1)],
               PACKAGE = "POT")$sim,
             amix = .C("rbvamix", nn, alpha, asCoef, sim = evmc[c(i,i-1)],
               PACKAGE = "POT")$sim)
  }

  switch(margins, frechet = -1/log(evmc), uniform = evmc,
         rweibull = log(evmc), gumbel = -log(-log(evmc)))
}
  
