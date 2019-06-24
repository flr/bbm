#-------------------------------------------------------------------------------  
# calcPop function: 
# function to calculate population dynamics (SSB and SSB1) from the B0, R and catch data
# 
# Created: Leire Ibaibarriaga - 2009-05-07
# Changed: 2018-05-29 11:44:12 (ssanchez)
#------------------------------------------------------------------------------- 

# calcPop_fun.r - function to estimate abundance in biomass for recruits and adults (given some information on growth, periods and an BBMpar object)
# bbm/R/calcPop_fun.r

# Copyright: European Union, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


# calcPop {{{

#' @title Function to estimate abundances
#'
#' @description This function estimates abundances in mass at age (for recruits and adults) by period, given information on growth, periods duration, catches 
#' and some additional values.
#'
#' @name calcPop
#' @rdname calcPop
#' @aliases calcPop
#'
#' @param g     A \code{vector} with information on instantaneous rate of biomass decrease, g = M - G, for recruits (rec) and adults (adult).
#' @param f     A \code{vector} with fraction of the year corresponding to each of the periods (determined by the period of the year of the different indices).
#' @param catch An \code{FLQuant} with catch information (for recruits and adults).
#' @param inits An \code{FLPar} with parameter values for the parameters required by the \code{bbm} function.
#' 
#' @return A \code{list} with two elements: \code{stock.bio}, with information on estimated stock (\code{FLQuant}); 
#'         and \code{ok}, an indicator on whether estimated parameters are valid (i.e. positive, \code{ok==TRUE}) or not (\code{ok==FALSE}).
#'
#' @author Leire Ibaibarriaga & Sonia Sanchez.
#' @seealso \link{FLPar}, \link{bbmFLPar}
#' @keywords calcPop
#'
#'  
#' @examples 
#' 
#' # Load required libraries
#' library(bbm)
#' 
#' # Load data
#' data(ane)
#' 
#' # Generate population estimates, given some estimated parameters
#' findicesB.ane <- unlist(lapply( indicesB.ane, function(x) mean(range(x)[c('startf','endf')])))
#' findicesP.ane <- unlist(lapply( indicesP.ane, function(x) mean(range(x)[c('startf','endf')])))
#' bioAge <- calcPop(g=control.ane@g, 
#'                   f=periods( findicesB=findicesB.ane, findicesP=findicesP.ane)$f, 
#'                   catch=catch.ane, inits=inits.ane)
#' class(bioAge)
#' 
#' # Check if valid ouput (i.e. positive biomass values)
#' bioAge$ok
#' 
#' # Estimates
#' bioAge$stock
#'  
#' 

calcPop <- function( g, f, catch, inits){
  
  # checkings
  
  if (!class(catch) %in% c('FLQuant'))
    stop( "'catch' argument must be of class 'FLQuant'.")
   
  if (class(inits)!='FLPar')
    stop( "'inits' argument must be of class 'FLPar'.")
  
  if( length(g) != 2 | is.null(names(g)) | sum(!names(g) %in% c("rec", "adult"))>0 )
    stop("Check growth parameters: g=c(rec=numeric(1), adult=numeric(1))")
  
  if(length(f)!=dim(catch)[4])
    stop("'f' length has be equal to the number of seasons in 'catch'.")
  
  if (any(dim(catch)[c(3,5)] !=1))
    stop("'catch' should have length 1 for dimensions unit and area.")
  
  # should we add more checkings on the dimensions?? iter, number of years, number of periods?
  
  # dimensions
  
  nyr  <- dim(catch)[2]  # number of years
  nper <- length(f)      # number of periods

  # definition 
  grec <- g['rec']    # instantaneous rate of biomass decrease for recruits
  gadu <- g['adult']  # instantaneous rate of biomass decrease for adults
  
  # create object to be filled in
  
  stock.bio <- FLQuant(dim=dim(catch), dimnames=dimnames(catch))
  stock.bio <- expand(stock.bio, season=1:(nper+1))
  
  # fill-in recruits and initial biomass
  
  stock.bio[1,,,1,,] <- inits[paste("R",dimnames(catch)$year,sep="_"), ]
  stock.bio[2,1,,1,,] <- inits["B0", ]

  # first year
  
  for (p in 2:(nper+1)) { # j=1; j<=nstep
    stock.bio[1,1,,p,,] = stock.bio[1,1,,p-1,,]*exp(-grec*f[p-1]) - catch[1,1,,p-1,,]*exp(-grec*f[p-1]/2.0)
    stock.bio[2,1,,p,,] = stock.bio[2,1,,p-1,,]*exp(-gadu*f[p-1]) - catch[2,1,,p-1,,]*exp(-gadu*f[p-1]/2.0)
  }
  
  # Rest of the years
  
  for (yr in 2:nyr){
    stock.bio[2,yr,,1,,] <- quantSums(stock.bio[,yr-1,,nper+1,,]) 
    for (p in 2:(nper+1)) { 
      stock.bio[1,yr,,p,,] = stock.bio[1,yr,,p-1,,]*exp(-grec*f[p-1]) - catch[1,yr,,p-1,,]*exp(-grec*f[p-1]/2.0)
      stock.bio[2,yr,,p,,] = stock.bio[2,yr,,p-1,,]*exp(-gadu*f[p-1]) - catch[2,yr,,p-1,,]*exp(-gadu*f[p-1]/2.0)
    }
  }
  
  # check if there are negative values
  if (any(stock.bio<0)) {
    ok <- FALSE
    if (any(stock.bio[,,,-(nper+1),,]<0)) {
      warning("Negative biomass values. Please check inputs.")
    } else
      warning("Negative biomass values at the end of the year (not shown in the output). Please check inputs.")
  } else{
    ok <- TRUE
  }

  stock.bio <- stock.bio[,,,-(nper+1),,]
  
  return( list(stock.bio=stock.bio, ok=ok))
  
}
