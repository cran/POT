#############################################################################
#   Copyright (c) 2014 Mathieu Ribatet                                                                                                  
#                                                                                                                                                                        
#   This program is free software; you can redistribute it and/or modify                                               
#   it under the terms of the GNU General Public License as published by                                         
#   the Free Software Foundation; either version 2 of the License, or                                                   
#   (at your option) any later version.                                                                                                            
#                                                                                                                                                                         
#   This program is distributed in the hope that it will be useful,                                                             
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 
#   GNU General Public License for more details.                                                                                    
#                                                                                                                                                                         
#   You should have received a copy of the GNU General Public License                                           
#   along with this program; if not, write to the                                                                                           
#   Free Software Foundation, Inc.,                                                                                                              
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             
#                                                                                                                                                                         
#############################################################################

##A function to perform anova for ``uvpot'' model selection. This is
##based on the Chisq statitics for two models M0 and M1 with M0 \in M1
anova.uvpot <- function(object, object2, ...){
  
  if(!inherits(object, "uvpot"))
    stop("Use only with 'uvpot' objects")
  #do not test object2 as in anova.lm()
  
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
                                    "Chisq", "Pr(>Chisq)"))

  structure(table, heading = "Analysis of Variance Table\n",
            class = c("anova", "data.frame"))
}
  

  
