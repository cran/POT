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
      switch(model, log = .C(POT_do_rbvlog, nn, alpha, sim = evmc[c(i,i-1)])$sim,
             alog = .C(POT_do_rbvalog, nn, alpha, asy, sim = evmc[c(i,i-1)])$sim,
             nlog = .C(POT_do_rbvnlog, nn, alpha, sim = evmc[c(i,i-1)])$sim,
             anlog = .C(POT_do_rbvanlog, nn, alpha, asy, sim = evmc[c(i,i-1)])$sim,
             mix = .C(POT_do_rbvmix, nn, alpha, sim = evmc[c(i,i-1)])$sim,
             amix = .C(POT_do_rbvamix, nn, alpha, asCoef, sim = evmc[c(i,i-1)])$sim)
  }

  switch(margins, frechet = -1/log(evmc), uniform = evmc,
         rweibull = log(evmc), gumbel = -log(-log(evmc)))
}
  
