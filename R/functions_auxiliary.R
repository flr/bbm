#-------------------------------------------------------------------------------  
# auxiliary functions:
# 
# Created: Sonia Sanchez - 2018-05-29 12:09:38
# Changed: 
#------------------------------------------------------------------------------- 

# functions_auxiliary.r: auxiliary functions required by the main functions
# bbm/R/functions_auxiliary.r

# Copyright: European Union, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


# AllEqual {{{
AllEqual <- function(x) { ### numeric, character vector, or time series of type ts
  res <- FALSE
  x <- na.omit(as.vector(x))
  if (length(unique(x)) == 1 | length(x) == 0) res <- TRUE
  return(res)
  ### The function returns TRUE if all values are equal and FALSE if it contains different values.
}
# }}} end AllEqual

