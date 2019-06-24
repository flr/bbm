#-------------------------------------------------------------------------------  
# createInits method.
# Created: Leire Ibaibarriaga - 2018-05-29 13:10:22
# Changed: 2018-06-01 11:59:17 (ssanchez)
#------------------------------------------------------------------------------- 

# createInits_method.r - method to generate initial values for bbm 
# bbm/R/createInits_method.r

# Copyright: European Union, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


# createInits {{{

#' @title Method createInits
#' 
#' @description Function for generating intial values for the parameters of the \code{bbm} function, 
#'              given information on catches, survey indices and instantaneous rate of biomass decrease, g.
#'
#' @name createInits
#' @rdname createInits
#' @aliases createInits
#'
#' @param object   An \code{FLQuant} with catch information (for recruits and adults) or an \code{FLStock}.
#' @param indicesB Abundance indices in total biomass (element of class: \code{FLQuants}, \code{FLQuant} or \code{FLIndices}) from surveys. 
#'                 Please assign a survey name to each index.
#' @param indicesP Percentage of recruits in biomass (element of class: \code{FLQuants}, \code{FLQuant} or \code{FLIndices}) from surveys. 
#'                 Please assign a survey name to each proportion.
#' @param findicesB A \code{vector} with fraction of the year corresponding to each of the indicesB.
#' @param findicesP A \code{vector} with fraction of the year corresponding to each of the indicesP.
#' @param g        A \code{vector} with information on instantaneous rate of biomass decrease, g = M - G, for recruits (rec) and adults (adult).
#' 
#' @return An object of class FLPar.
#'
#' @section Methods:
#' Methods exist for various calculations based on the output class (\code{FLPar}). For details: \code{?FLPar}.
#'
#' 
#' @author Leire Ibaibarriaga & Sonia Sanchez.
#' @seealso \link{bbm}, \linkS4class{FLQuant}, \linkS4class{FLQuants}, \link{FLIndices}, \link{FLPar}, \link{bbmFLPar}
#' @keywords FLPar methods
#'
#'  
#' @examples 
#' 
#' # Load data
#' data(ane)
#' 

 
#' @rdname createInits
#' @aliases createInits,FLQuant,FLQuants,FLQuants-method
#' @examples
#'
#' # Case: object='FLQuant'; indicesB=indicesP='FLQuants'
#' inits1 <- createInits( catch.ane,
#'                        indicesB=lapply( indicesB.ane, function(x) x@index), 
#'                        indicesP=lapply( indicesP.ane, function(x) x@index), 
#'                        findicesB=unlist(lapply( indicesB.ane, function(x) mean(range(x)[c('startf','endf')]))),
#'                        findicesP=unlist(lapply( indicesP.ane, function(x) mean(range(x)[c('startf','endf')]))), 
#'                        g=control.ane@g )
#' class(inits1)

# createInits(object='FLQuant', indicesB='FLQuants', indicesP='FLQuants',...) {{{
setMethod('createInits', signature(object='FLQuant', indicesB='FLQuants', indicesP='FLQuants'),
          function(object, indicesB, indicesP, findicesB, findicesP, g)
            {

            catch <- object; rm(object)

            # DATA CHECKS

            bbmChecks( catch, indicesB=indicesB, indicesP=indicesP, findicesB=findicesB, findicesP=findicesP)
            
            if( length(g) != 2 | is.null(names(g)) | sum(!names(g) %in% c("rec", "adult"))>0 ) # check names of slot g
              stop("Check growth parameters: g=c(rec=numeric(1), adult=numeric(1))")
            
            # number of periods and fraction of each period
            
            aux <- periods(findicesB, findicesP)
            nper     <- aux$nper
            f        <- aux$f
            perindicesB <- aux$perindicesB
            perindicesP <- aux$perindicesP
            
            
            nyr   <- dim(catch)[2] # number of years
            niter <- dims(catch)$iter
            
            yrs <- dimnames(catch)$year
            
            if (dim(catch)[4]==1 & nper>1)  { # extend to the adequate number of seasons
              catch <- expand( catch, season=1:nper)
              for (s in nper:1)
                catch[,,,s,,] <- catch[,,,1,,] * f[s]
            }
            
            nsurvB <- length(indicesB)
            nsurvP <- length(indicesP)
            
            namesB <- names(indicesB)
            namesP <- names(indicesP)
            
            # back-calculated recruitment
            
            Rback <- array(NA, dim=c(nyr, nsurvB, niter), dimnames=list(year=dimnames(indicesB[[1]])$year,survey=namesB,iter=1:niter))
            for (i in 1:nsurvB) {
              surv <- namesB[i]
              if( any(namesP %in% surv)) {
                Rback[,i,] <- indicesB[[surv]]*indicesP[[surv]] * exp(g['rec'] * f[perindicesP[surv]-1]) + 
                             seasonSums(catch[1,,,1:(perindicesP[surv]-1),,]) * exp(g['rec'] * f[perindicesP[surv]-1]/2)
              }
            }
            
            R <- apply(Rback, c(1,3), mean, na.rm=T)
            R[is.nan(R)] <- NA
            for (it in 1:niter) R[is.na(R[,it]),it] <- median(R[,it], na.rm=T) # if there is no data on recruits from any index in a given year, we assume a median recruitment
            
            B0back <- array(NA, dim=c(nyr, nsurvB, niter), dimnames=list(year=dimnames(indicesB[[1]])$year,survey=namesB,iter=1:niter))
            for (i in 1:nsurvB) {
              surv <- namesB[i]
              if( any(namesP %in% surv)) {
                B0back[,i,] <- indicesB[[surv]]*(1-indicesP[[surv]]) * exp(g['adult'] * f[perindicesP[surv]-1]) + 
                               seasonSums(catch[2,,,1:(perindicesP[surv]-1),,]) * exp(g['adult'] * f[perindicesP[surv]-1]/2)
              }
            }
            
            B0 <- numeric(niter)
            for (it in 1:niter)
              if(any(!is.na(B0back[1,,it]))){
              B0[it] <- mean( B0back[1,,it], na.rm=T)   # if there is adults data in the first year, back calculate adults to get B0
            }else{
              B0[it] <- median(apply(B0back[,,it], 1, mean, na.rm=T), na.rm=T) # otherwise, use a median of the back-calculated adults from the whole time series
            }
            
            inits <- bbmFLPar( years=yrs, namesB=namesB, namesP=namesP, niter=niter)
            
            parnames <- sapply(dimnames(inits)$params, FUN=function(x) unlist(strsplit(x,split="_"))[1])
            
            inits[parnames %in% "R"] <- R
            inits["B0"]              <- B0
            
            pop <- suppressWarnings(calcPop(g, f, catch, inits))
            
            if (pop$ok==FALSE){
              # check <- logical(niter)
              for (i in 1:niter) {
                popit <- suppressWarnings(calcPop(g, f, catch[,,,,,i], inits[,i]))
                if (popit$ok==TRUE) {
                  next()
                } else {
                  k <- 1
                  while (popit$ok==FALSE) {
                    inits[parnames %in% "R",i] <- R[,i]*(1+k*0.05)
                    inits["B0",i]              <- B0[i]*(1+k*0.05)
                    popit <- suppressWarnings(calcPop(g, f, catch[,,,,,i], inits[,i]))
                    k <- k+1
                  }
                  pop$stock.bio[,,,,,i] <- popit$stock.bio
                }
              }
              pop$ok <- TRUE
            }
            
            logR <- log(R)
            
            inits["mur"]  <- apply( logR, 2, mean)
            inits["psir"] <- nyr/apply( (logR[,it]-t(matrix( inits["mur"], ncol=nyr, nrow=niter)))^2, 2, sum)

            for (i in 1:nsurvB) {
              inits[parnames %in% "q"][i,]   <- exp( apply(log(quantSums(indicesB[[i]])) - log(quantSums(pop$stock)[,,,perindicesB[i],,]), 6, mean, na.rm=T))
              inits[parnames %in% "psi"][i,] <- 1/apply( (log(quantSums(indicesB[[i]])) - log(inits[parnames %in% "q"][i,]) - 
                                                           log(quantSums(pop$stock)[,,,perindicesB[i],,]))^2, 6, mean, na.rm=T)
            }
            
            for (i in 1:nsurvP) {
              fun <- function(x, iter){ # first derivative of the log likelihood wrt to xi for each survey.
                pobs <- indicesP[[i]][,,,,,iter]
                p <- pop$stock[1,,,perindicesP[i],,iter]/quantSums(pop$stock[,,,perindicesP[i],,iter])
                out <- digamma(exp(x)) - digamma(exp(x)*p)*p - digamma(exp(x)*(1-p)) * (1-p) + p*log(pobs) + (1-p) * log(1-pobs)
                out <- sum(out, na.rm=T)
                return(out)
              }
              for (it in 1:niter) inits[parnames %in% "xi"][i,it] <- uniroot(fun, interval=c(-50, 50), iter=it)$root
            }
            
            return(inits)

            }
)
# }}}

#---------------------
#---------------------

#' @rdname createInits
#' @aliases createInits,FLStock,ANY,ANY-method
#' @examples
#'
#' # Case: object='FLQuant'; indicesB=indicesP='ANY'
#' 
#' stock <- FLStock(catch.n=catch.ane, catch.wt=catch.ane*0+1)
#' units(stock@catch.wt) <- ''
#' stock@catch <- quantSums(stock@catch.n*stock@catch.wt)
#' 
#' inits2 <- createInits( stock, indicesB=indicesB.ane, indicesP=indicesP.ane, 
#'                        g=control.ane@g )
#' class(inits2)

# createInits(object='FLStock',...)
setMethod('createInits', signature(object='FLStock', indicesB='ANY', indicesP='ANY'),
          function(object, indicesB, indicesP, findicesB=NULL, findicesP=NULL, g)
          {

            # check the FLStock has only two age classess
            if(dim(object)[1]!=2)
              stop("'object' must have only two age classes: recruits and adults")

            # calculate total catches at age (i.e. for recruits and adults)
            catch <- object@catch.n*object@catch.wt

            # createInits(object='FLQuant', indicesB='FLQuant', indicesP='FLQuant', findicesB='vector', findicesP='vector', g='vector')
            return(createInits(catch, indicesB=indicesB, indicesP=indicesP, findicesB=findicesB, findicesP=findicesP, g=g))

          }
)

#---------------------
#---------------------

#' @rdname createInits
#' @aliases createInits,FLQuant,FLIndices,FLIndices-method
#' @examples
#'
#' # Case: object='FLQuant'; indicesB=indicesP='FLIndices'
#' 
#' inits3 <- createInits( catch.ane,
#'                        indicesB=indicesB.ane, indicesP=indicesP.ane, 
#'                        g=control.ane@g )
#' class(inits3)

# createInits(indicesB='FLIndices',indicesP='FLIndices',...)
setMethod('createInits', signature(object='FLQuant', indicesB='FLIndices', indicesP='FLIndices'),
          function(object, indicesB, indicesP, findicesB=NULL, findicesP=NULL, g)
          {

            # Calculate f:
            findicesB <-  unlist(lapply( indicesB, function(x) mean(range(x)[c('startf','endf')])))
            findicesP <-  unlist(lapply( indicesP, function(x) mean(range(x)[c('startf','endf')])))

            # From FLIndices to FLQuants:
            idxs      <- indicesB
            idxs.nam  <- names(indicesB)
            idxsP     <- indicesP
            idxsP.nam <- names(indicesP)
            indicesB <- indicesP <- FLQuants()
            for (i in 1:length(idxs))  indicesB[[i]] <- idxs[[idxs.nam[i]]]@index
            names(indicesB) <- idxs.nam
            for (i in 1:length(idxsP)) indicesP[[i]] <- idxsP[[idxsP.nam[i]]]@index
            names(indicesP) <- idxsP.nam

            # createInits(object='FLQuant', indicesB='FLQuant', indicesP='FLQuant', findicesB='vector', findicesP='vector', g='vector')
            return(createInits(object, indicesB=indicesB, indicesP=indicesP, findicesB=findicesB, findicesP=findicesP, g=g))

          }
)

#---------------------
#---------------------

#' @rdname createInits
#' @aliases createInits,FLStock,FLQuant,FLQuant-method
#' @examples
#'
#' # Case: object='FLQuant'; indicesB=indicesP='FLQuant'
#' inits4 <- createInits( catch.ane,
#'                        indicesB=indicesB.ane[[1]]@index, indicesP=indicesP.ane[[1]]@index, 
#'                        findicesB=c( depm=(indicesB.ane[[1]]@range[['startf']]+indicesB.ane[[1]]@range[['endf']])/2),
#'                        findicesP=c( depm=(indicesP.ane[[1]]@range[['startf']]+indicesP.ane[[1]]@range[['endf']])/2), 
#'                        g=control.ane@g )
#' class(inits4)
#' 
#' # Run assessment (with the different initial values)
#' run0 <- bbm(catch.ane, indicesB=indicesB.ane, indicesP=indicesP.ane, control=control.ane, inits=inits.ane)
#' run1 <- bbm(catch.ane, indicesB=indicesB.ane, indicesP=indicesP.ane, control=control.ane, inits=inits1)
#' run2 <- bbm(catch.ane, indicesB=indicesB.ane, indicesP=indicesP.ane, control=control.ane, inits=inits2)
#' run3 <- bbm(catch.ane, indicesB=indicesB.ane, indicesP=indicesP.ane, control=control.ane, inits=inits3)
#' namdel <- c("q_acoustic","psi_acoustic","xi_acoustic") # we will take only one of the indices --> need to delete the parameters related to other indices
#' control <- control.ane
#' control@param.fix <- control@param.fix[dimnames(control@param.fix)$params[!dimnames(control@param.fix)$params %in% namdel],]
#' run4 <- bbm(catch.ane, indicesB=indicesB.ane[[1]]@index, indicesP=indicesP.ane[[1]]@index, 
#'             findicesB=c( depm=(indicesB.ane[[1]]@range[['startf']]+indicesB.ane[[1]]@range[['endf']])/2),
#'             findicesP=c( depm=(indicesP.ane[[1]]@range[['startf']]+indicesP.ane[[1]]@range[['endf']])/2), 
#'             control=control, inits=inits4)
#' 
#' # Plot assessed populations
#' biomass <- FLQuants()
#' runs <- paste("run",0:4,sep="")
#' names(runs) <- c('bc','run1','run2','run3','only_depm')
#' for (i in 1:length(runs)) biomass[[i]] <- quantSums(stock.bio(get(runs[i])))
#' names(biomass) <- names(runs)
#' plot( biomass)
#' 

# createInits(indicesB='FQuant',indicesP='FQuant',...)
setMethod('createInits', signature(object='FLQuant', indicesB='FLQuant', indicesP='FLQuant'),
          function(object, indicesB, indicesP, findicesB, findicesP, g)
          {

            if (is.null(findicesB)) {
              stop("'findicesB' is missing.")
            } else if (length(findicesB)!=1 | is.null(names(findicesB)) ) {
              stop("'findicesB' must be a named vector of length 1.")
            } else {
              indicesB <- FLQuants(indicesB)
              names(indicesB) <- names(findicesB)
            }

            if (is.null(findicesP)) {
              stop("'findicesP' is missing.")
            } else if (length(findicesP)!=1 | is.null(names(findicesP)) ) {
              stop("'findicesP' must be a named vector of length 1.")
            } else {
              indicesP <- FLQuants(indicesP)
              names(indicesP) <- names(findicesP)
            }

            # createInits(object='FLQuant', indicesB='FLQuant', indicesP='FLQuant', findicesB='vector', findicesP='vector', g='vector')
            return(createInits(object, indicesB=indicesB, indicesP=indicesP, findicesB=findicesB, findicesP=findicesP, g=g))

          }
)

