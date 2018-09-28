#-------------------------------------------------------------------------------  
# initsBBM function.
# Created: Leire Ibaibarriaga - 2018-05-29 13:10:22
# Changed: 2018-06-01 11:59:17 (ssanchez)
#------------------------------------------------------------------------------- 

# initsBBM.r - function to generate initial values for BBM
# bbm/R//initsBBM.r

# Copyright: European Union & AZTI, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# initsBBM {{{

#' initsBBM class for generating intial values.
#' 
#' Function to generate initial parameters required for the BBM function, given information on catches, survey indices 
#' and instantaneous rate of biomass decrease, g.
#'
#' @name initsBBM
#' @rdname initsBBM
#' @aliases initsBBM initsBBM-methods
#'
#' @param object  An \code{FLQuant} with catch information (for recruits and adults) or an \code{FLStock}.
#' @param indices Abundance indices in biomass for same age classes (class: FLQuants).
#' @param g       A \code{vector} with information on instantaneous rate of biomass decrease, g = M - G, for recruits (rec) and adults (adult).
#' 
#' @return An object of class BBMpar.
#'
#' @section Methods:
#' Methods exist for various calculations based on the output class (\code{BBMpar}). For details: \code{?BBMpar}.
#'
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
#' inits <- initsBBM(object=catch.ane, indices=indices.ane, g=control.ane@g)
#' class(inits)
#' 
#' # Run assessment
#' run <- BBM(object=catch.ane, indices=indices.ane, control=control.ane, inits=inits)
#' 
#' # Plot assessed population
#' \dontrun{
#'   plot(stock(run))
#' }
#'
 
initsBBM <- function(object, indices, g){
  
  catch <- object; rm(object)
  
  # DATA CHECKS
  
  check <- BBM.checks( object=catch, indices=indices)
  
  if( length(g) != 2 | is.null(names(g)) | sum(!names(g) %in% c("rec", "adult"))>0 ) # check names of slot g
    stop("Check growth parameters: g=c(rec=numeric(1), adult=numeric(1))")
  
  # number of periods and fraction of each period
  
  aux <- periods(indices) 
  nper     <- aux$nper
  f        <- aux$f
  indexper <- aux$indexper
  
  n <- dim(catch)[2] # number of years
  
  yrs <- dimnames(catch)$year
  
  if (dim(catch)[4]==1 & nper>1)  { # extend to the adequate number of seasons
    catch <- expand( catch, season=1:nper)
    for (s in nper:1)
      catch[,,,s,,] <- catch[,,,1,,] * f[s]
  }
  
  nsurv <- length(indices)
  
  # we assume indices and catch have the same number of years. this could be immproved
  
  # back-calculated recruitment
  
  Rback <- matrix(NA, ncol=nsurv, nrow=n)
  for (i in 1:nsurv){
    aux <- indices[[i]]@index[1,,,,,] * exp(g['rec'] * f[indexper[i]-1]) + seasonSums(catch[1,,,1:(indexper[i]-1),,])* exp(g['rec'] * f[indexper[i]-1]/2)
    Rback[,i] <- aux[1,,1,1,1,1,drop=T]
  }
  colnames(Rback) <- names(indices)
  rownames(Rback) <- dimnames(indices[[i]]@index)$year
  
  logR <- log(apply(Rback, 1, mean, na.rm=T))
  logR[is.na(logR)] <- median(logR, na.rm=T) # if there is no data on recruits from any index in a given year, we assume a median recruitment
  
  B0back <- matrix(NA, ncol=nsurv, nrow=n)
  for (i in 1:nsurv){
    aux <- indices[[i]]@index[2,,,,,] * exp(g['adult'] * f[indexper[i]-1]) + seasonSums(catch[2,,,1:(indexper[i]-1),,])* exp(g['adult'] * f[indexper[i]-1]/2)
    B0back[,i] <- aux[1,,1,1,1,1,drop=T]
  }
  colnames(B0back) <- names(indices)
  rownames(B0back) <- dimnames(indices[[i]]@index)$year
  
  if(any(!is.na(B0back[1,]))){
    logB0 <- as.numeric(log(apply(B0back, 1, mean, na.rm=T))[1])   # if there is adults data in the first year, back calculate adults to get B0
  }else{
    logB0 <- median(log(apply(B0back, 1, mean, na.rm=T)), na.rm=T) # otherwise, use a median of the back-calculated adults from the whole time series
  }
  
  inits <- new( "BBMpar", logR=logR, logB0=logB0)
  
  pop <- calcPop(g, f, catch, inits)
  
  if (!pop$ok){
    stop("Not appropriate initial values.")
  }
  
  mur <- mean(logR)
  psir <- n/sum( (logR-mur)^2 )
  
  logq <- rep(NA, nsurv)
  psi <- rep(NA, nsurv)
  xi <- rep(NA, nsurv)
  for (i in 1:nsurv){
    logq[i] <- mean(log(quantSums(indices[[i]]@index)) - log(quantSums(pop$stock)[,,,indexper[i],,]), na.rm=T)
    psi[i] <- 1/mean( (log(quantSums(indices[[i]]@index)) - logq[i] - log(quantSums(pop$stock)[,,,indexper[i],,]))^2 , na.rm=T)
    
    fun <- function(x){ # first derivative of the log likelihood wrt to xi for each survey. 
      pobs <- indices[[i]]@index[1,,,,,]/quantSums(indices[[i]]@index)
      p <- pop$stock[1,,,indexper[i],,]/quantSums(pop$stock[,,,indexper[i],,])
      out <- digamma(exp(x)) - digamma(exp(x)*p)*p - digamma(exp(x)*(1-p)) * (1-p) + p*log(pobs) + (1-p) * log(1-pobs)
      out <- sum(out, na.rm=T)
      return(out)
    }
    xi[i] <- uniroot(fun, interval=c(-50, 50))$root
  }
  
  return(new( "BBMpar", logq=logq, logpsi=log(psi), xi=xi, logB0=logB0, logR=logR, mur=mur, logpsir=log(psir)))
  
}


# initsBBM(FLStock)
setMethod('initsBBM', signature(object='FLStock'),
          function(object, indices, g)
          {
            
            # check the FLStock has only two age classess
            if(dim(object)[1]!=2)
              stop("'object' must have only two age classes: recruits and adults")
            
            # calculate total catches at age (i.e. for recruits and adults)
            catch <- object@catch.n*object@catch.wt
            
            # initsBBM(catches='FLQuant', indices='FLIndices', g='BBM.control@g')
            return(initsBBM(catch, indices=indices, g=g))
            
          }
)

