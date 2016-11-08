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

print.uvpot <- function(x, digits = max(3, getOption("digits") - 3), ...){
  
  if (!inherits(x, "uvpot"))
    stop("Use only with 'uvpot' objects")
  cat("Estimator:", x$est, "\n")

  if (x$est == "MGF")
    cat("Statistic:", x$stat, "\n")
  
  if (x$est == 'MLE'){
    cat(" Deviance:", x$deviance, "\n")
    cat("      AIC:", AIC(x), "\n")
  }

  if (x$est == 'MPLE'){
    cat("\nPenalized Deviance:", x$deviance, "\n")
    cat("     Penalized AIC:", AIC(x), "\n")
  }
  
  cat("\nVarying Threshold:", x$var.thresh, "\n")
  
  if(!x$var.thresh)
    x$threshold <- x$threshold[1]
  
  cat("\n  Threshold Call:", x$threshold.call, "\n")
  cat("    Number Above:", x$nat, "\n")
  cat("Proportion Above:", round(x$pat, digits), "\n")
  
  cat("\nEstimates\n") 
  print.default(format(x$fitted.values, digits = digits), print.gap = 2, 
                quote = FALSE)
  if(!is.null(x$std.err)) {
    cat("\nStandard Error Type:", x$std.err.type, "\n")
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  if(!is.null(x$var.cov)) {
    cat("\nAsymptotic Variance Covariance\n")
    print.default(format(x$var.cov, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  if(!is.null(x$corr)) {
    cat("\nCorrelation\n")
    print.default(format(x$corr, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  cat("\nOptimization Information\n")
  cat("  Convergence:", x$convergence, "\n")
  cat("  Function Evaluations:", x$counts["function"], "\n")
  if(!is.na(x$counts["gradient"]))
    cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
  if(!is.null(x$message)) cat("\nMessage:", x$message, "\n")
  cat("\n")
}

