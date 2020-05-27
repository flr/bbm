#-------------------------------------------------------------------------------  
# periods function:
# 
# Created: Sonia Sanchez - 2018-05-29 12:09:38
# Changed: 
#------------------------------------------------------------------------------- 

# peridos_fun.r: function to estimate the seasons
# bbm/R/peridos_fun.r

# Copyright: European Union, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


#' @title Function to estimate the seasons
#'
#' @description Given the fraction of the year when each of the indices are observed, this function estimates the number of seasons, 
#' the fraction of the year of each season and the season when each of the indices are observed.
#'
#' @name periods
#' @rdname periods
#' @aliases periods
#'
#' @param findicesB A named \code{vector} with fraction of the year when each of the total biomass indices are observed.
#' @param findicesP A named \code{vector} with fraction of the year when each of the proportion of recruits indices are observed.
#' 
#' @return A \code{list} with four numeric elements: \code{nper}, with the number of seasons; 
#'         \code{f} with the fraction of the year corresponding to each of the seasons;
#'         \code{perindicesB} with the season in which the index in total biomass is observed, 
#'          this object is a named vector with one element per index; and 
#'         \code{perindicesP} with the season in which the index of proportion of recruits is observed, 
#'          this object is a named vector with one element per index.
#'
#' @author Leire Ibaibarriaga & Sonia Sanchez.
#' @seealso \link{}
#' @keywords calcPop
#'
#'  
#' @examples 
#' 
#' # Load data
#' data(ane)
#' 
#' # Generate population estimates, given some estimated parameters
#' per <- periods( findicesB=unlist(lapply( indicesB.ane, function(x) mean(range(x)[c('startf','endf')]))), 
#'                 findicesP=unlist(lapply( indicesP.ane, function(x) mean(range(x)[c('startf','endf')]))))
#' class(per)
#' 
#' per$nper
#' per$f
#' per$perindicesB
#' per$perindicesB
#' 

periods <- function(findicesB=NULL, findicesP=NULL) {
  
  if(any(findicesB<0 | findicesB>=1))
    stop("'findicesB' must be in the interval [0,1).")
  if(any(findicesP<0 | findicesP>=1))
    stop("'findicesP' must be in the interval [0,1).")

  per   <- sort(unique(c(findicesB, findicesP, 0, 1))) # periods' limits
  nper  <- length(per)-1               # number of periods
  f     <- diff(per)                   # fraction of the year corresponding to each period (i.e. length of each period)
  
  perindicesB <- match(findicesB, per)
  names(perindicesB) <- names(findicesB)

  perindicesP <- match(findicesP, per)
  names(perindicesP) <- names(findicesP)
  
  out <- list( nper=nper, f=f, perindicesB=perindicesB, perindicesP=perindicesP)
  
  return(out)
  
}

