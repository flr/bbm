#-------------------------------------------------------------------------------  
# calcIndices function: 
# function to calculate Indices (indicesB and indicesP) from the B0, R, q and catch data
# 
# Created: Leire Ibaibarriaga - 2009-05-07
# Changed: 2018-05-29 11:44:12 (ssanchez)
#------------------------------------------------------------------------------- 

# calcIndices_fun.r - function to calculate indicesB and indicesP (given some information on growth, periods, catch and an FLPar object)
# bbm/R/calcIndices_fun.r

# Copyright: European Union, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


# calcIndices {{{

#' @title Function to estimate fitted indices
#'
#' @description This function estimates fitted indices, given information on growth, periods duration, catches and estimated parameters. 
#'
#' @name calcIndices
#' @rdname calcIndices
#' @aliases calcIndices
#'
#' @param g         A \code{vector} with information on instantaneous rate of biomass decrease, g = M - G, for recruits (rec) and adults (adult).
#' @param findicesB A \code{vector} with fraction of the year corresponding to each of the indicesB.
#' @param findicesP A \code{vector} with fraction of the year corresponding to each of the indicesP.
#' @param catch     An \code{FLQuant} with catch information (for recruits and adults).
#' @param inits     An \code{FLPar} with parameter values for the parameters required by the \code{bbm} function.
#' 
#' @return A \code{list} with information on estimated indices (\code{FLQuants}). 
#'         The list has 2 elements: \code{indicesB} (for indices in biomass) and \code{indicesP} (for proportions of recruits in biomass). 
#'
#' @section Methods:
#' Methods exist for various calculations based on the output class (\code{FLPar}). For details: \code{?FLPar}.
#' 
#' 
#' @author Leire Ibaibarriaga & Sonia Sanchez.
#' @seealso \link{FLPar}, \link{bbmFLPar}, \link{bbmControl}
#' @keywords calcIndices
#'
#'  
#' @examples 
#' 
#' # Load data
#' data(ane)
#' 
#' args(calcIndices)
#' 
#' indices <- calcIndices( g=control.ane@g, 
#'                         findicesB=unlist(lapply( indicesB.ane, function(x) mean(range(x)[c('startf','endf')]))),
#'                         findicesP=unlist(lapply( indicesP.ane, function(x) mean(range(x)[c('startf','endf')]))),
#'                         catch=catch.ane, inits=inits.ane )
#' 
#' class(indices)
#' slotNames(indices)
#' 

calcIndices <- function( g, findicesB, findicesP, catch, inits){
  
  # checkings
  
  if (!class(catch) %in% c('FLQuant'))
    stop( "'catch' argument must be of class 'FLQuant'.")
   
  if (class(inits)!='FLPar')
    stop( "'inits' argument must be of class 'FLPar'.")
  
  if( length(g) != 2 | is.null(names(g)) | sum(!names(g) %in% c("rec", "adult"))>0 )
    stop("Check growth parameters: g=c(rec=numeric(1), adult=numeric(1))")
  
  if (any(dim(catch)[c(3,5)] !=1))
    stop("'catch' should have length 1 for dimensions unit and area.")
  
  aux <- periods(findicesB=findicesB, findicesP=findicesP)
  nper <- aux$nper
  f    <- aux$f
  perindicesB <- aux$perindicesB
  perindicesP <- aux$perindicesP 
  
  stock.bio <- calcPop(g=g, f=f, catch=catch, inits=inits)$stock.bio
  
  flq <- FLQuant(dimnames=list(age="all", year=dimnames(catch)$year, iter=dimnames(catch)$iter))
  indicesB <- indicesP <- FLQuants()
  for (i in 1:length(findicesB)){
    indicesB[[i]]   <- flq
    indicesB[[i]][] <- inits[grep(paste("^q_",names(findicesB)[i],sep=""),rownames(inits)), ] * quantSums(stock.bio[,,,perindicesB[i],,])
  }
  for (i in 1:length(findicesP)){
    indicesP[[i]]   <- flq
    indicesP[[i]][] <- stock.bio[1,,,perindicesP[i],,] / quantSums(stock.bio[,,,perindicesP[i],,])
  }
  names(indicesB) <- names(findicesB)
  names(indicesP) <- names(findicesP)
  
  return( list(indicesB=indicesB, indicesP=indicesP))
}

