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



tailind.test <- function(data, c = -0.1, emp.trans = TRUE,
                         chisq.n.class = 4){

  if (ncol(data) != 2)
    stop("data must be a matrix with two columns")

  if (emp.trans){
    ##We have to transform data to reverse exponential
    ##margins using empirical transformation
    n <- nrow(data)
    data <- apply(data, 2, rank) / (n + 1)
    data <- log(data)
  }

  else
    cat("``data'' is supposed to be reverse exponential distributed\n")

  ##The margin sum
  Msum <- apply(data, 1, sum)
  idx <- which(Msum > c)
  Msum <- Msum[idx]
  n.c <- length(Msum)
  V <- Msum / c


  ######################################################
  ##                                                  ##
  ##          The Neyman-Pearson Test                 ##
  ##                                                  ##
  ######################################################
  
  NP.stat <- -log(prod(2*V))
  if (n.c <= 170){
    NP.pval <- 1
    for (j in 1:(n.c-1))
      NP.pval <- NP.pval + (-2 * sum(log(V)))^j / prod(j:2)
  
    NP.pval <- exp(2*sum(log(V))) * NP.pval
  }

  else
    NP.pval <- pnorm((2 * sum(log(V)) + n.c) / sqrt(n.c))
  
  

  ######################################################
  ##                                                  ##
  ##        Fisher's Kappa Statistic Test             ##
  ##                                                  ##
  ######################################################

  U <- (1 - (1 - Msum) * exp(Msum)) /
    (1 - (1 - c) * exp(c))
  Usort <- c(0, sort(U))
  
  Fish.stat <- (n.c + 1) * max(c(diff(Usort), 1 - max(Usort)))

  
  Fish.pval <- 1

  for (j in 1:n.c)
    Fish.pval <- Fish.pval + (-1)^j * prod(2:(n.c + 1)) /
      prod(1:j) / prod(2:(n.c + 1 - j)) *
        max(0, 1 - j * Fish.stat / (n.c  + 1))^n.c 
  
  Fish.pval <- Fish.pval + (-1)^(n.c + 1) *
    max(0, 1 - Fish.stat)^n.c 
  Fish.pval <- 1 - Fish.pval

  
  ######################################################
  ##                                                  ##
  ##            Kolmogorov-Smirnov Test               ##
  ##                                                  ##
  ######################################################


  KStest <- ks.test(U, "punif")

  KS.pval <- KStest$p.value
  KS.stat <- KStest$statistic

  
  ######################################################
  ##                                                  ##
  ##                 Chi-Square Test                  ##
  ##                                                  ##
  ######################################################

  classes <- seq(0, 1, length.out = chisq.n.class + 1)

  m <- NULL
  for (bsup in classes[-1])
    m <- c(m, sum(U < bsup))

  m[-1] <- diff(m)
  p.theo <- diff(classes)

  if (any( (n.c * p.theo) < 6))
    warning("Classes for the Chi-Square test are ill-defined")
      
  ChiSq <- chisq.test(m, p = p.theo)
  ChiSq.pval <- ChiSq$p.value
  ChiSq.stat <- ChiSq$statistic
  

  p.values <- c(NP.pval, Fish.pval,KS.pval, ChiSq.pval)
  stats <- c(NP.stat, Fish.stat, KS.stat, ChiSq.stat)
  names(p.values) <- names(stats) <- c("NP", "Fish", "KS",
                                       "ChiSq")

  return(list(stats = cbind(stats, p.values), idx = idx))
  
}

