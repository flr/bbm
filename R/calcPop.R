#-------------------------------------------------------------------------------  
# calcPop function: 
# function to calculate population dynamics (SSB and SSB1) from the B0, R and catch data
# 
# Created: Leire Ibaibarriaga - 2009-05-07
# Changed: 2018-05-29 11:44:12 (ssanchez)
#------------------------------------------------------------------------------- 

# calcPop.r - function to estimate abundance in biomass for recruits and adults (given some information on growth, periods and an BBMpar object)
# bbm/R/calcPop.r

# Copyright: European Union & AZTI, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


# calcPop {{{

#' Function to estimate abundances in mass at age by period.
#'
#' This function estimates abundances in mass at age (for recruits and adults) by period, given information on growth, periods duration, catches 
#' and some additional values.
#'
#' @name calcPop
#' @rdname calcPop
#' @aliases calcPop
#'
#' @param g     A \code{vector} with information on instantaneous rate of biomass decrease, g = M - G, for recruits (rec) and adults (adult).
#' @param f     A \code{vector} with fraction of the year corresponding to each of the periods (determined by the period of the year of the different indices).
#' @param catch An \code{FLQuant} with catch information (for recruits and adults).
#' @param inits An \code{BBMpar} with initial values for the parameters required by the BBM function.
#' 
#' @return A \code{list(stock, ok)} with information on estimated stock (\code{FLQuant}) 
#'         and an indicator on wheter estimated parameters are valid (i.e. positive, \code{ok==TRUE}) or not (\code{ok==FALSE}).
#'
# @section Methods:
# Methods exist for various calculations based on the output class (\code{BBMpar}). For details: \code{?BBMpar}.
#
#' 
#' @author Leire Ibaibarriaga & Sonia Sanchez.
#' @seealso \link{BBMpar}, \link{BBM.control}
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
#' bioAge <- calcPop(g=control.ane@g, f=periods(indices.ane)$f, catch=catch.ane, inits=inits.ane)
#' class(bioAge)
#' 
#' # check if valid ouput (i.e. positive biomass values)
#' bioAge$ok
#' 
#' # estimates
#' bioAge$stock
#'  
#' 

calcPop <- function( g, f, catch, inits){
  
  # checkings
  
  if (!class(catch) %in% c('FLQuant'))
    stop( "'catch' argument must be of class 'FLQuant'.")
   
  if (class(inits)!='BBMpar')
    stop( "'inits' argument must be of class 'BBMpar'.")
  
  # check names of slot g
  if( length(g) != 2 | is.null(names(g)) | sum(!names(g) %in% c("rec", "adult"))>0 )
    stop("Check growth parameters: g=c(rec=numeric(1), adult=numeric(1))")
  
  if(length(f)!=dim(catch)[4])
    stop("'f' length has be equal to the number of seasons in 'catch'.")
  
  
  # ABUNDANCES
  
  nyr  <- dim(catch)[2]  # number of years
  nper <- length(f)      # number of periods
  
  grec <- g['rec']    # instantaneous rate of biomass decrease for recruits
  gadu <- g['adult']  # instantaneous rate of biomass decrease for adults
  
  Crec <- catch[1,1:nyr,1,1:nper,1, drop=T]  # recruits catch [nyr,nper]
  Cadu <- catch[2,1:nyr,1,1:nper,1, drop=T]  # adults catch [nyr,nper]
  
  R <- exp(inits@logR)
  Brec <- Badu <- Btot <- matrix(nrow=nyr, ncol=nper+1) # abundances [nyr,nper]
  dimnames(Brec) <- dimnames(Badu) <- dimnames(Btot) <- list(year=dimnames(catch)$year,season=1:(nper+1))
  
  # First year
  
  Brec[1,1] <- R[1]
  Badu[1,1] <- exp(inits@logB0)
  Btot[1,1] = Brec[1,1] + Badu[1,1]
  
  for (p in 2:(nper+1)) { # j=1; j<=nstep
    Brec[1,p] = Brec[1,p-1]*exp(-grec*f[p-1]) - Crec[1,p-1]*exp(-grec*f[p-1]/2.0)
    Badu[1,p] = Badu[1,p-1]*exp(-gadu*f[p-1]) - Cadu[1,p-1]*exp(-gadu*f[p-1]/2.0)
    Btot[1,p] = Brec[1,p] + Badu[1,p]
  }
  
  # Rest of the years
  
  for (yr in 2:nyr){ #i=1; i<nyr
    
    Brec[yr,1] <- R[yr]
    Badu[yr,1] <- Btot[yr-1,nper]
    Btot[yr,1] = Brec[yr,1] + Badu[yr,1]
    
    for (p in 2:(nper+1)) { # j=1; j<=nstep
      Brec[yr,p] = Brec[yr,p-1]*exp(-grec*f[p-1]) - Crec[yr,p-1]*exp(-grec*f[p-1]/2.0)
      Badu[yr,p] = Badu[yr,p-1]*exp(-gadu*f[p-1]) - Cadu[yr,p-1]*exp(-gadu*f[p-1]/2.0)
      Btot[yr,p] = Brec[yr,p] + Badu[yr,p];
    }
  }
  
  if (any(Brec<0) | any(Badu<0)) {
    ok <- FALSE
    if (any(Brec[,-(nper+1)]<0) | any(Badu[,-(nper+1)]<0)) {
      warning("Negative biomass values. Please check inputs.")
    } else 
      warning("Negative biomass values at the end of the year (not shown in the output). Please check inputs.")
  } else
      ok <- TRUE

  stk <- catch*NA
  stk[1,] <- Brec[,-(nper+1)]
  stk[2,] <- Badu[,-(nper+1)]

  return( list(stock=stk, ok=ok))

}
