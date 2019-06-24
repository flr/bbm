#-------------------------------------------------------------------------------  
# Bay of Biscay anchovy data.
# Created: Sonia Sanchez - 2018-06-21 14:36:46
# Changed: 
#------------------------------------------------------------------------------- 

# BoBane_data.R - Data from Ibaibarriaga et al. (2008)
# bbm/R/BoBane_data.R

# Copyright: European Union, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# DATA

#' @title bbm package datasets
#' 
#' @description Data of Bay of Biscay anchovy from Ibaibarriaga et al. (2008).
#' 
#' Example dataset for the classes defined in bbm package.
#' At the moment there is only one dataset of Bay of Biscay anchvoy (ICES Subarea 8), 
#' with information on catches together with surveys observations and timing.
#' These are the same data as used in Ibaibarriaga et al. (2008).
#' 
#' Dataset can be loaded by issuing the \code{data} command, like in
#' \code{data(ane)}. and has defined the following objects:
#' 
#' \itemize{
#'    \item{\code{catch.ane}, \code{\link{FLQuant}}}
#'         {Catch information in biomass for recruits (age 1 individuals) and adults.}
#'    \item{\code{indicesB.ane}, \code{\link{FLIndices}}}
#'         {Observations from surveys of total biomass, along with survey timing.}
#'    \item{\code{indicesP.ane}, \code{\link{FLIndices}}}
#'         {Observations from surveys of proportions of recruits in biomass, along with survey timing.}
#'    \item{\code{control.ane}, \code{\link{bbmControl}}}
#'         {Object with information on instantaneous rate of biomass decrease (g = M - G, \code{vector}) 
#'         and information on whether bbm parameters are fixed or not (param.fix, \code{\link{FLPar}}).}        
#'    \item{\code{inits.ane}, \code{\link{FLPar}}}
#'         {Initial values for the parameters of the bbm assessment model.}
#' }
#' 
#' 
#' @docType data
#' @keywords datasets
# @format An object of class CLASS with
#' @name datasets
#' @rdname datasets
#' 
#' @aliases ane
#' @seealso \link{bbm}, \linkS4class{FLQuant}, \linkS4class{FLQuants}, \linkS4class{FLStock}, 
#'          \linkS4class{FLIndices}, \linkS4class{bbmControl}, \linkS4class{FLPar}, 
#'          \link{bbmFLPar}
#' @references Ibaibarriaga et al. (2008). "A Two-Stage Biomass Dynamic Model for Bay of Biscay Anchovy: A Bayesian Approach." 
#'             ICES J. of Mar. Sci. 65: 191-205.
#' @keywords datasets
#' @examples
#'
#' data(ane)
#' ls()
#'
#' is(catch.ane)
#' catch.ane
#' 
#' is(inits.ane)
#' inits.ane
#' 
#' is(control.ane)
#' control.ane
#' 
#' is(indicesB.ane)
#' indicesB.ane
#' is(indicesP.ane)
#' indicesP.ane
#' 

NULL
