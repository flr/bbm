#-------------------------------------------------------------------------------  
# simIndices function.
# Created: Leire Ibaibarriaga - 2018-05-29 13:10:22
# Changed: 2018-06-01 11:59:17 (ssanchez)
#------------------------------------------------------------------------------- 

# simIndices_method.r - method to generate initial values for BBM
# bbm/R/simIndices_method.r

# Copyright: European Union, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.



# simIndices {{{

#' @title Method simIndices
#' 
#' @description  Function to generate indices of total biomass and of proportion of recruits in biomass, 
#' given information on catches at age, instantaneous rate of biomass decrease (g = M - G ) and values for the bbm parameters.
#'
#' @name simIndices
#' @rdname simIndices
#' @aliases simIndices
#'
#' @param object    An \code{FLQuant} with catch information (for recruits and adults).
#' @param g         A \code{vector} with information on instantaneous rate of biomass decrease, g = M - G, for recruits (rec) and adults (adult).
#' @param inits     An \code{FLPar} with initial values.
#' @param findicesB A \code{vector} with period of the year of the different indices in total biomass. 
#'                  Optional parameter, if not set, then the survey is assumed to occur at the beggining of the year.
#' @param findicesP A \code{vector} with period of the year of the different indices of proportion of recruits in mass. 
#'                  Optional parameter, if not set, then the survey is assumed to occur at the beggining of the year.
#' 
#' @return A list with indices in biomass (Btot) and indices in proportion of recruits (Prec), 
#'         both elements of the list are FLIndices.
#'
#' @section Methods: 
#' Methods exist for various calculations based on the output class (\code{FLPar}). For details: \code{?FLPar}.
#
#' 
#' @author Leire Ibaibarriaga & Sonia Sanchez.
#' @seealso \link{bbm}, \linkS4class{FLQuant}, \linkS4class{FLQuants}, \link{FLIndices}, \linkS4class{bbmControl}, \linkS4class{FLPar}, \link{bbmFLPar}
#' @keywords FLPar methods
#'
#'  
#' @examples 
#' 
#' # Load data
#' data(ane)



#' @rdname simIndices
#' @aliases simIndices,FLQuant-method
#' @examples
#'
#' # Case:  object='FLQuant'
#' indices1 <- simIndices( catch.ane, g=control.ane@g, inits=inits.ane, 
#'                         findicesB=unlist(lapply( indicesB.ane, function(x) mean(range(x)[c('startf','endf')]))),
#'                         findicesP=unlist(lapply( indicesP.ane, function(x) mean(range(x)[c('startf','endf')]))) )
#' class(indices1)
#' slotNames(indices1)
#' 

# simIndices(object='FLQuant',...) {{{
setMethod('simIndices', signature(object='FLQuant'),
          function(object, g, inits, findicesB=NULL, findicesP=NULL)
          {
            
            catch <- object; rm(object)
            
            # check input objects
            
            bbmChecks( catch, findicesB=findicesB, findicesP=findicesP, inits=inits)
            
            nyear <- dim(catch)[2]
            niter <- dim(catch)[6]
            
            nindicesBtot <- length(iter(inits, 1)[sapply(dimnames(inits)$params, FUN=function(x) unlist(strsplit(x,split="_"))[1]) %in% "q"])
            if (is.null(findicesB)){
              findicesB <- rep(0, nindicesBtot)
            }else{
              # checkings
              if(any(findicesB<0 | findicesB>=1))
                stop("'findicesB' must be in the interval [0,1).")
              if (length(findicesB) != nindicesBtot){
                stop("findicesB does not have the proper dimension")
              }
              if(!is.null(names(findicesB))){
                namBtot <- names(findicesB)
              }else{
                namBtot <- paste("index", 1:nindicesBtot, sep="")
              }
            }
            
            nindicesP <- length(iter(inits, 1)[sapply(dimnames(inits)$params, FUN=function(x) unlist(strsplit(x,split="_"))[1]) %in% "xi"])
            if (is.null(findicesP)){
              findicesP <- rep(0, nindicesP)
            }else{
              # checkings
              if(any(findicesP<0 | findicesP>=1))
                stop("'findicesP' must be in the interval [0,1).")
              if (length(findicesP) != nindicesP){
                stop("findicesP does not have the proper dimension")
              }
              if(!is.null(names(findicesP))){
                namPrec <- names(findicesP)
              }else{
                namPrec <- paste("index", 1:nindicesP, sep="")
              }
            }
            
            per <- sort(unique(c(findicesB, findicesP, 0, 1)))
            nper <- length(per) - 1
            f <- diff(per)
            indexperBtot <- match(findicesB, per)
            indexperPrec <- match(findicesP, per)
            
            if (dim(catch)[4]==1 & nper>1)  { # extend to the adequate number of seasons
              catch <- expand( catch, season=1:nper)
              for (s in nper:1)
                catch[,,,s,,] <- catch[,,,1,,] * f[s]
            }
            
            pop <- calcPop(g, f, catch, inits) 
            
            if(!pop$ok){
              stop("These values lead to negative biomasses")
            }
            
            indicesB <- indicesP <- FLIndices() # how do I create an empty FLINdices of a given length?
            
            btot.flq <- p.flq <- FLQuant( dimnames=list(age='all', year=dimnames(catch)$year), iter=niter)
            
            inits.it <- lapply( 1:niter, function(x) prepareFLPar(iter(inits,x)) )
            
            for (i in 1:nindicesBtot){ #! PENDIENTE: necesario separar Btot y Prec
              for (it in 1:niter) {
                btot.flq[,,,,,it] <- rlnorm(1, inits.it[[it]]$logq[i] + log( quantSums(pop$stock[,,,indexperBtot[i],,it])), 1/sqrt(exp(inits.it[[it]]$logpsi[i])))
              }
              indicesB[[i]] <- FLIndex(index=btot.flq)
              indicesB[[i]]@range[c("startf", "endf")] <- findicesB[i]
            }
            names(indicesB) <- namBtot
            
            for (i in 1:nindicesP){ #! PENDIENTE: necesario separar Btot y Prec
              for (it in 1:niter) {
                1
                p.flq[,,,,,it] <- rbeta(nyear, exp(inits.it[[it]]$xi[i])*pop$stock[1,,,indexperPrec[i],,it], exp(inits.it[[it]]$xi[i])*pop$stock[2,,,indexperPrec[i],,it])
              }
              indicesP[[i]] <- FLIndex(index=p.flq)
              indicesP[[i]]@range[c("startf", "endf")] <- findicesP[i]
            }
            names(indicesP) <- namPrec
            
            indices <- list( Btot=indicesB, Prec=indicesP)
            
            return(indices)
            
          }
)
# }}}

#---------------------
#---------------------

#' @rdname simIndices
#' @aliases simIndices,FLQuant-method
#' @examples
#'
#' # Case:  object='FLStock'
#' stock <- FLStock(catch.n=catch.ane, catch.wt=catch.ane*0+1)
#' units(stock@catch.wt) <- ''
#' stock@catch <- quantSums(stock@catch.n*stock@catch.wt)
#' 
#' indices2 <- simIndices( stock, g=control.ane@g, inits=inits.ane, 
#'                         findicesB=unlist(lapply( indicesB.ane, function(x) mean(range(x)[c('startf','endf')]))),
#'                         findicesP=unlist(lapply( indicesP.ane, function(x) mean(range(x)[c('startf','endf')]))) )
#' class(indices2)
#' 
#' # Run assessment with the alternative indices
#' run  <- bbm(catch.ane, indicesB=indicesB.ane, indicesP=indicesP.ane, control=control.ane, inits=inits.ane)
#' run1 <- bbm(catch.ane, indicesB=indices1$Btot, indicesP=indices1$Prec, control=control.ane, inits=inits.ane)
#' run2 <- bbm(catch.ane, indicesB=indices2$Btot, indicesP=indices2$Prec, control=control.ane, inits=inits.ane)
#' 
#' # Plot assessed populations
#' plot( FLQuants( bc=quantSums(run@stock.bio)[,,,1,], alt1=quantSums(run1@stock.bio)[,,,1,], alt2=quantSums(run2@stock.bio)[,,,1,]))
#' 

# simIndices(object='FLStock',...)
setMethod('simIndices', signature(object='FLStock'),
          function(object, g, inits, findicesB=NULL, findicesP=NULL)
          {
            
            # check the FLStock has only two age classess
            if(dim(object)[1]!=2)
              stop("'object' must have only two age classes: recruits and adults")
            
            # calculate total catches at age (i.e. for recruits and adults)
            catch <- object@catch.n*object@catch.wt
            
            # bbm(object='FLQuant', indicesB='FLQuant', indicesP='FLQuant', findicesB='vector', findicesP='vector', control='bbmControl', inits='FLPar')
            return(simIndices(catch, g=g, inits=inits, findicesB=findicesB, findicesP=findicesP))
            
          }
)

