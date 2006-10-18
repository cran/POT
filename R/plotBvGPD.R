pickdep <- function(fitted, main, bound = TRUE, plot = TRUE,
                    ...){

  model <- fitted$model
  alpha <- fitted$param["alpha"]

  if (model == "log"){
    ##Logistic case :
    A <- function(w){
      if ( (w < 0) || (w > 1) )
        return(NaN)
      else
        ((1-w)^(1/alpha) + w^(1/alpha))^alpha
    }
  }
  
  if (model == "nlog"){
    ##Negative logistic case:
    A <- function(w){
      if ( (w < 0) || (w > 1) )
        return(NaN)
      else
        1 - ((1-w)^(-alpha) + w^(-alpha))^(-1/alpha)
    }
  }

  if (model == "alog"){
    ##Asymetric logistic case:
    asCoef1 <- fitted$param["asCoef1"]
    asCoef2 <- fitted$param["asCoef2"]
    A <- function(w){
      if ( (w < 0) || (w > 1) )
        return(NaN)
      else
        (1 - asCoef1)*(1-w) + (1 - asCoef2) * w +
          ( (asCoef1 * (1-w))^(1/alpha) +
           (asCoef2 * w)^(1/alpha) )^alpha
    }
  }

  if (model == "anlog"){
    ##Asymetric negatif logistic case:
    asCoef1 <- fitted$param["asCoef1"]
    asCoef2 <- fitted$param["asCoef2"]
    A <- function(w){
      if ( (w < 0) || (w > 1) )
        return(NaN)
      else
        1 - ( ((1-w)/asCoef1)^(-alpha) +
             (w/asCoef2)^(-alpha) )^(-1/alpha)
    }
  }

  if (model == "mix"){
    ##Mixed model:
    A <- function(w){
      if ( (w < 0) || (w > 1) )
        return(NaN)
      else
        1 - w * (1-w) * alpha
    }
  }

   if (model == "amix"){
    ##Asymetric Mixed model:
     asCoef <- fitted$param["asCoef"]
    A <- function(w){
      if ( (w < 0) || (w > 1) )
        return(NaN)
      else
        1 - (alpha + asCoef) * w + alpha * w^2 +
          asCoef * w^3
    }
  }

  if (plot){
    if (missing(main))
        main <- "Pickands' Dependence Function"
    
    plot(A, ylim = c(0.5, 1), main = main, type = "n")
    
    if (bound){
      lines(x= c(0,1), y = c(1,1), col = "grey", ...)
      lines(x = c(0,0.5,1), y = c(1, 0.5, 1), col = "grey", ...)
    }
    plot(A, ylim = c(0.5, 1), add = TRUE)
  }
  
  attributes(A) <- list(model = model)
  invisible(A)
}


bvretlev <- function(fitted, p = seq(0.75,0.95,0.05), main,
                     n = 5000, only.excess = FALSE, ...){
  
  if (all(class(fitted) != "bvpot"))
    stop("``fitted'' should be an object of class ``bvpot''.")

  if (missing(main))
    main <- "Bivariate Return Level Plot"

  w <- c(0, seq(0,1, length.out = n), 1)
  ##The Pickands' dependence function
  A <- pickdep(fitted, plot = FALSE)(w)
  y1 <- y2 <- matrix(NA, ncol = length(p), nrow = n + 2)
  colnames(y1) <- colnames(y2) <- p
  rownames(y1) <- rownames(y2) <- 1:(n+2)

  for (i in 1:length(p)){
    z1 <- - A / (w * log(p[i]))
    z2 <- w / (1 - w) * z1
    
    ##We have to transform frechet observation to original
    ##scale
    y1[,i] <- ((1-exp(-1/z1)) / fitted$pat[1])^
    (-fitted$param["shape1"]) - 1
    y1[,i] <- fitted$threshold[1] + fitted$param["scale1"] /
      fitted$param["shape1"] * y1[,i]
    y2[,i] <- ((1-exp(-1/z2)) / fitted$pat[2])^
    (-fitted$param["shape2"]) - 1
    y2[,i] <- fitted$threshold[2] + fitted$param["scale2"] /
      fitted$param["shape2"] * y2[,i]
  }

  data1 <- fitted$data[,1]
  data2 <- fitted$data[,2]

  if (only.excess){
    idx <- which(data1 > fitted$threshold[1] |
                 data2 > fitted$threshold[2])
    data1 <- data1[idx]
    data2 <- data2[idx]
  }
  
  plot(data1, data2, main = main, ...)
  lines(y1, y2)

  ##Add the marginal threshold axis
  abline(v = fitted$threshold[1])
  abline(h = fitted$threshold[2])

  invisible(list(y1 = y1, y2 = y2))

}

plot.bvpot <- function(x, mains, which = 1:2,
                       ask = nb.fig < length(which) &&
                       dev.interactive(), ...){

  if (!is.numeric(which) || any(which < 1) || any(which > 2)) 
        stop("`which' must be in 1:2")

  if(missing(mains))
    mains <- c("Pickands' Dependence Function",
               "Bivariate Return Level Plot")

  else
    if (length(mains) != 2){
      warning("``mains'' must be of length two. Passing to default titles.")
      mains <- c("Pickands' Dependence Function",
                 "Bivariate Return Level Plot")
    }


  show <- rep(FALSE, 2)
  show[which] <- TRUE
  nb.fig <- prod(par("mfcol"))
  
  if (ask){
    op <- par(ask = TRUE)
    on.exit(par(op))
  }

  if (show[1])
    pickdep(x, main = mains[1], ...)

  if (show[2])
    bvretlev(x, main = mains[2], ...)

}
