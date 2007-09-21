##A function to transform GPD observations to a unit Frechet scale
gpd2frech <- function(x, loc = 0, scale = 1, shape = 0, pat = 1){
  z <- pgpd(x, loc, scale, shape, lower.tail = FALSE)
  z <- -1 / log(1 - pat * z)
  return(z)
}

##A function to transform unit Frechet observations to GPD ones
frech2gpd <- function(z, loc = 0, scale = 1, shape = 0, pat = 1){
  ##First, convert to uniform(0,1)
  x <- exp(-1/z) / pat
  x <- qgpd(x, loc, scale, shape)
  return(x)  
}
