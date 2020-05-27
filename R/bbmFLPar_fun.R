#-------------------------------------------------------------------------------  
# bbmFLPar function: 
# function to generate an FLPar object with the input parameters of the bbm function
# 
# Created: Sonia Sanchez - 2018-06-29 11:27:48
# Changed: 
#------------------------------------------------------------------------------- 

# bbmFLPar_fun.r - function to generate an FLPar object with the input parameters required by the bbm function
# bbm/R/bbmFLPar_fun.r

# Copyright: European Union, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


# bbmFLPar {{{

#' @title Function to generate an FLPar for bbm function
#'
#' @description This function generates an \code{FLPar} object with the input parametes required by the \code{bbm} function,
#'              given information on the years, the indices names and the number of iterations. 
#'
#' @name bbmFLPar
#' @rdname bbmFLPar
#' @aliases bbmFLPar
#' 
#' @param x        Input numeric values for the parameters. If not set, then NA value is assigned.
#' @param years    Names of the years used for fitting.
#' @param namesB   Names of the indices of total biomass.
#' @param namesP   Names of the indices of proportion of recruits in biomass.
#' @param niter    Number of iterations.
#' @param logscale Logical, if TRUE the parameters are in the scale used by the bbm model, 
#'                 otherwise, all parameters are in the linear scale. This is used for example, 
#'                 in the vcov matrix retuned from \code{bbm} function.
#' 
#' @return An \code{FLPar} with the appropriate format for the \code{bbm} function.
#'
#' @author Leire Ibaibarriaga & Sonia Sanchez.
#' @seealso \link{FLPar}, \link{bbm}
#' @keywords bbmFLPar
#'
#'  
#' @examples 
#' 
#' # Load data
#' data(ane)
#' 
#' years.ane <- dimnames(catch.ane)$year
#' niter.ane <- dim(catch.ane)[6]
#' namesB.ane <- names(indicesB.ane)
#' namesP.ane <- names(indicesP.ane)
#' 
#' # Generate population estimates, given some estimated parameters
#' pars <- bbmFLPar( years=years.ane, namesB=namesB.ane, namesP=namesP.ane, niter=niter.ane)
#' class(pars)
#' pars
#'  


bbmFLPar <- function( x=NULL, years, namesB, namesP, niter=1, logscale=FALSE) {

  if (logscale==FALSE) {
    out <- FLPar( NA, dimnames=list(params=c(paste('q',namesB,sep='_'),paste('psi',namesB,sep='_'),
                                             paste('xi',namesP,sep='_'),'B0',
                                             paste('R',years,sep='_'),'mur','psir'), iter=1:niter), units='NA')
  } else {
    out <- FLPar( NA, dimnames=list(params=c(paste('logq',namesB,sep='_'),paste('logpsi',namesB,sep='_'),
                                             paste('xi',namesP,sep='_'),'logB0',
                                             paste('logR',years,sep='_'),'mur','logpsir'), iter=1:niter), units='NA')
  }
  
  if (!is.null(x)) out[] <- x
  
  return(out)
  
}

