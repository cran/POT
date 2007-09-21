##A function to perform anova for ``uvpot'' model selection. This is
##based on the Chisq statitics for two models M0 and M1 with M0 \in M1
anova.uvpot <- function(object, object2, ...){
  
  ##Check if object and object2 are fitted by MLE
  if (object$est != "MLE")
    stop("``object'' is not a MLE.")
  if (object2$est != "MLE")
    stop("``object2'' is not a MLE.")

  ##Check if object and object2 are nested
  n.estim <- length(fitted(object))
  n.estim2 <- length(fitted(object2))

  if (n.estim == n.estim2)
    stop("Models are not nested.")

  else{
    if (n.estim > n.estim2){
      M0 <- object2
      M1 <- object
      model0 <- deparse(substitute(object2))
      model1 <- deparse(substitute(object))
    }

    else{
      M0 <- object
      M1 <- object2
      model0 <- deparse(substitute(object))
      model1 <- deparse(substitute(object2))
    }
  }

  models <- c(model0, model1)
  Dev <- c(M0$deviance, M1$deviance)
  diffDev <- -diff(Dev)
  MDf <- c(length(fitted(M0)), length(fitted(M1)))
  Df <- diff(MDf)

  pvalue <- pchisq(diffDev, Df, lower.tail = FALSE)

  
  table <- data.frame(MDf, Dev, c(NA, Df), c(NA, diffDev),
                      c(NA, pvalue))

  dimnames(table) <- list(models, c("MDf", "Deviance", "Df",
                                    "Chisq", "Pr[>Chisq]"))

  structure(table, heading = "Analysis of Variance Table\n",
            class = c("anova", "data.frame"))
}
  
##A function to perform anova for ``bvpot'' model selection. This is
##based on the Chisq statitics for two models M0 and M1 with M0 \in M1
anova.bvpot <- function(object, object2, ..., half = FALSE){

  ##Check if object and object2 are nested.
  n.estim <- length(fitted(object))
  n.estim2 <- length(fitted(object2))

  if (n.estim > n.estim2){
    M0 <- object2
    M1 <- object
    model0 <- deparse(substitute(object2))
    model1 <- deparse(substitute(object))
  }
  
  else{
    M0 <- object
    M1 <- object2
    model0 <- deparse(substitute(object))
    model1 <- deparse(substitute(object2))
  }

  depFam0 <- M0$model
  depFam1 <- M1$model

  if ((depFam0 == depFam1) ||
      (paste("a", depFam0, sep="") == depFam1) ||
      (paste("a", depFam1, sep="") == depFam0)){

    models <- c(model0, model1)
    Dev <- c(M0$deviance, M1$deviance)
    diffDev <- -diff(Dev)
    if (half) diffDev <- 2 * diffDev
    MDf <- c(length(fitted(M0)), length(fitted(M1)))
    Df <- diff(MDf)
    
    pvalue <- pchisq(diffDev, Df, lower.tail = FALSE)
    
    
    table <- data.frame(MDf, Dev, c(NA, Df), c(NA, diffDev),
                        c(NA, pvalue))
    
    dimnames(table) <- list(models, c("MDf", "Deviance", "Df",
                                      "Chisq", "Pr[>Chisq]"))
    
    structure(table, heading = "Analysis of Variance Table\n",
              class = c("anova", "data.frame"))
  }

  else
    stop("Models are not nested.")
}
  
