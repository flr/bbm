#-------------------------------------------------------------------------------  
# simBBM function.
# Created: Leire Ibaibarriaga - 2018-05-29 13:10:22
# Changed: 2018-06-01 11:59:17 (ssanchez)
#------------------------------------------------------------------------------- 

# simBBM.r - function to generate initial values for BBM
# bbm/R//simBBM.r

# Copyright: European Union & AZTI, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.



# simBBM {{{

#' Function to generate indices.
#' 
#' This function generates indices from a population (with only 2 age classes, recruits and adults), given information on catches at age,   and survey indices.
#'
#' @name simBBM
#' @rdname simBBM
#' @aliases simBBM simBBM-methods
#'
#' @param catch An \code{FLQuant} with catch information (for recruits and adults).
#' @param g     A \code{vector} with information on instantaneous rate of biomass decrease, g = M - G, for recruits (rec) and adults (adult).
#' @param inits An \code{BBMpar} with initial values.
#' @param f     A \code{vector} with fraction of the year corresponding to each of the periods (determined by the period of the year of the different indices).
#' 
#' @return An object of class FLInidices.
#'
# @section Methods:
# Methods exist for various calculations based on the output class (\code{BBMpar}). For details: \code{?BBMpar}.
#
#' 
#' @author Leire Ibaibarriaga & Sonia Sanchez.
#' @seealso \link{BBM}, \link{BBMpar}
#' @keywords BBMpar
#'
#'  
#' @examples 
#' 
#' # Load required libraries
#' library(bbm)
#' library(ggplotFL)
#' 
#' # Load data
#' data(ane)
#' 
#' # Generate intitial values
#' indices <- simBBM( catch=catch.ane, g=control.ane@g, inits=inits.ane, 
#'                    findex=unlist(lapply( indices.ane, function(x) (x@range[['startf']]+x@range[['endf']])/2)))
#' class(indices)
#' 
#' # Run assessment
#' run  <- BBM(catch.ane, indices=indices.ane, control=control.ane, inits=inits.ane)
#' run2 <- BBM(catch.ane, indices=indices, control=control.ane, inits=inits.ane)
#' 
#' # Plot assessed population
#' \dontrun{
#'   plot(FLQuants(bc=quantSums(run@stock.bio)[,,,1,], alt=quantSums(run2@stock.bio)[,,,1,]))
#' }
#' 


simBBM <- function(catch, g, inits, findex=NULL){
 
  nyear  <- dim(catch)[2]
  nindex <- length(inits@logq)
  
  if (is.null(findex)){
    findex <- rep(0, nindex)
  }else{
    # checkings
    if(any(findex<0 | findex>=1))
      stop("'findex' must be in the interval [0,1).")
    if (length(findex) != nindex){
      stop("findex does not have the proper dimension")
    }
    if(!is.null(names(findex))){
      index.nam <- names(findex)
    }else{
      index.nam <- paste("index", 1:nindex, sep="")
    }
  }
  
  per <- sort(unique(c(findex, 0, 1)))
  nper <- length(per) - 1
  f <- diff(per)
  indexper <- match(findex, per)

  if (dim(catch)[4]==1 & nper>1)  { # extend to the adequate number of seasons
    catch <- expand( catch, season=1:nper)
    for (s in nper:1)
      catch[,,,s,,] <- catch[,,,1,,] * f[s]
  }
  
  pop <- calcPop(g, f, catch, inits) 
  
  if(!pop$ok){
    stop("These values lead to negative biomasses")
  }

  flqa <- FLQuant( dimnames=list(age=ac(1:2), year=dimnames(catch)$year))
  
  p.flq <- FLQuant( dimnames=list(age='all', year=dimnames(catch)$year))

  indices <- FLIndices() # how do I create an empty FLINdices of a given length?
  
  for (i in 1:nindex){
    btot.flq <- rlnorm(1, inits@logq[i] + log( quantSums(pop$stock[,,,indexper[i],,])), 1/sqrt(exp(inits@logpsi[i])))
    p.flq <- FLQuant( dimnames=list(age='all', year=dimnames(catch)$year))
    p.flq[]    <- rbeta(nyear, exp(inits@xi[i])*pop$stock[1,,,indexper[i],,], exp(inits@xi[i])*pop$stock[2,,,indexper[i],,])
    flqa[1,,,,,] <- btot.flq * p.flq
    flqa[2,,,,,] <- btot.flq * (1-p.flq)
    indices[[i]] <- FLIndex(index=flqa, name=index.nam[i])
    indices[[i]]@range[c("startf", "endf")] <- findex[i]
  }
  return(indices)
}  
  
