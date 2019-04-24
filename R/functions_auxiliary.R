# functions_auxiliary.R: auxiliary functions required by the main functions
# bbm/R/functions_auxiliary.R

# Copyright: European Union, 2018
# Authors: Leire Ibaibarriaga (AZTI), <libaibarriaga@azti.es>
#          Sonia Sanchez (AZTI) <ssanchez@azti.es>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# periods {{{
periods <- function(indices) {
  
  if (class(indices)!='FLIndices')
    stop("'indices' must be of class FLIndices")
  
  findex <- unlist(lapply( indices, function(x) (x@range[['startf']]+x@range[['endf']])/2)) 
  # checkings
  if(any(findex<0 | findex>=1))
    stop("'findex' must be in the interval [0,1).")
  per   <- sort(unique(c(findex,0,1))) # periods' limits
  nper  <- length(per)-1               # number of periods
  f     <- diff(per)                   # fraction of the year corresponding to each period (i.e. length of each period)
  
  indexper <- match(findex, per)
  names(indexper) <- names(findex)
  
  out <- list( nper=nper, f=f, indexper=indexper)
  
  return(out)
  
}
# }}} end periods

# AllEqual {{{
AllEqual <- function(x) { ### numeric, character vector, or time series of type ts
  res <- FALSE
  x <- na.omit(as.vector(x))
  if (length(unique(x)) == 1 | length(x) == 0) res <- TRUE
  return(res)
  ### The function returns TRUE if all values are equal and FALSE if it contains different values.
}
# }}} end AllEqual
