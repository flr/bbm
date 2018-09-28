#-------------------------------------------------------------------------------  
# BBMfit class.
# Created: Sonia Sanchez - 2018-05-22 13:49:23
# Changed: 2018-06-01 10:24:36 (ssanchez)
#------------------------------------------------------------------------------- 

# Copyright: European Union, 2018
# Authors:
#   Leire Ibaibarriaga <libaibarriaga@azti.es>
#   Sonia Sanchez (AZTI) <ssanchez@azti.es>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# BBMfit {{{

#' BBMfit class for storing the output of the BBM function.
#'
#' This includes abundance estimates in biomass (for recruits and adults) and information on the model fitting.
#' 
#' @name BBMfit
#' @rdname BBMfit
#' @docType class
#' @aliases BBMfit BBMfit-methods BBMfit-class
#' @section Slots:
#'  \describe{
#'     \item{stock.bio  }{: Estimated stock biomass for recruits and adults in the different seasons, where seasons are dertermined by the index times. 
#'                          \code{FLQuant} with two age classes: recruits and adults.  }
#'     \item{params     }{: Estimated BBM parameters, \code{FLPar}. #! BASTARIA CON PONER O ESCALA NORMAL Y DESCARTA INFO EN ESCALA LOG? }
#'     \item{params.se  }{: Standard errors in parameters' estimates. \code{FLPar}. #! IDEM A PARAM}
#'     \item{data       }{: Input data. \code{list}.
#'     \item{method     }{: Optimisation library and optimisation method. \code{"TMB - method"}.}
#'     \item{nopar      }{: Number of parameters, \code{numeric}.}
#'     \item{noobs      }{: Number of observations, \code{numeric}.}
#'     \item{nlogl      }{: Negative log-likelihood, \code{numeric}.}
#'     \item{counts     }{: A two-element integer vector giving the number of calls to fn (i.e."bbm" function) and to gr (i.e. the optimisation method) respectively. 
#'                        This excludes those calls needed to compute the Hessian, if requested, and any calls to fn to compute a finite-difference approximation to the gradient.}
#'     \item{convergence}{: An integer code. 0 indicates successful completion. For other possible error codes see \code{?optim}.}
#'     \item{message    }{: Character string giving any additional information returned by the optimizer, or "".}
#'     \item{info       }{: Matrix with information on libraries' versions and platform used and object creation time.}
#' }
# @section Validity: 
#  \describe{
#     \item{XXX}{XXX}
# }
# You can inspect the class validity function by using
#    \code{getValidity(getClassDef('BBMfit'))}
#
# @section Accessors:
# All slots in the class have accessor and replacement methods defined that
# allow retrieving and substituting individual slots.
#
# The values passed for replacement need to be of the class of that slot.
# A numeric vector can also be used when replacing FLQuant slots, and the
# vector will be used to substitute the values in the slot, but not its other
# attributes.
#'
#' @section Constructor:
#' A construction method exists for this class that can take named arguments for
#' any of its slots. All slots are then created to match the requirements of the
#' class validity. If \code{dim} or \code{dimnames} are provided, this is used
#' for sizing the stock.bio slot.
#'
#' @section Methods:
#' Methods exist for various calculations based on values stored in the class:
#'
#' \describe{
#'     \item{BBMfit   }{: Creation of a new BBMfit object.}
#'     \item{stock    }{: Extracts information on estimated abundances in biomass for recruits and adults. Returns an \code{FLQuant}.}
#'     \item{+        }{: Updates an \code{FLStock} with new information on the BBM assessment.}
#'     \item{residuals}{: XXX.}
#'     \item{xxx}{: XXX.}
#'     \item{xxx}{: XXX.}
#'     \item{xxx}{: XXX.}
#' }
#' 
#' 
#' @author Leire Ibaibarriaga & Sonia Sanchez
#' @seealso \link{BBM}
#' @keywords classes
#' 
#' @examples
#' 
#' # Load data
#' data(ane)
#' 
#' # Run assessment (output is of class BBMfit)
#' run <- BBM(catch.ane, indices=indices.ane, control=control.ane, inits=inits.ane)
#' class(run)
#' run
#' 
#' # Available methods for the class
#' # - stock()
#' stock(run)
#' run@stock.bio
#' 
#' # - residuals
#' residuals(run)
#' 
#' # - AIC and BIC
#' AIC(run)
#' BIC(run)
#' 
#'


setClass("BBMfit",
         representation(
           "list",
           stock.bio   = "FLQuant",
           params      = "BBMpar",
           params.se   = "BBMpar",
           data        = "list",
           method      = "character",
           nopar       = "numeric",
           nobs        = "numeric",
           nlogl       = "numeric",
           counts      = "integer",
           convergence = "integer",
           message     = "character",
           info        = "matrix"
          ),
         prototype=prototype(
           stock.bio   = FLQuant(),
           params       = BBMpar(),
           params.se    = BBMpar(),
           data        = list(),
           method      = "TMB",
           nopar       = as.numeric(NA),
           nobs        = as.numeric(NA),
           nlogl       = as.numeric(NA),
           counts      = as.integer(NA),
           convergence = as.integer(NA),
           message     = "",
           info        = matrix( t(data.frame(FLbbm.version=packageDescription("bbm")$Version,
                                    FLCore.version=packageDescription("FLCore")$Version,
                                    TMB.version=packageDescription("TMB")$Version,
                                    R.version=R.version$version.string,
                                    platform=R.version$platform,
                                    run.date=Sys.time())), 
                                 dimnames=list(c("FLbbm.version","FLCore.version","TMB.version","R.version","platform","run.date"),""))
          ),
         validity=function(object){
           
           # Everything is fine
           return(TRUE)
          
          }
) # }}}

# BBMfit()   {{{
setMethod('BBMfit', signature(object='missing'),
          function(object, dim="missing", dimnames="missing", ...)
          {
            args <- list(...)

            if (missing(dim) && missing(dimnames)) {
              flq <- FLQuant()
            } else if (missing(dim) | missing(dimnames)) {
              if (missing(dim))
                flq <- FLQuant(dimnames=dimnames)
              if (missing(dimnames))
                  flq <- FLQuant(dim=dim)
            } else {
              flq <- FLQuant(dim=dim, dimnames=dimnames)
            }
            
            res <- new("BBMfit", stock.bio = flq)
            
            # Load given slots
            if(length(args) > 0)
              for(i in names(args))
                slot(res, i) <- args[[i]]
            
            return(res)
          }
) # }}}

# stock	{{{
setMethod('stock', signature(object='BBMfit'),
          function(object, ...)
            return(object@stock.bio)
) # }}}

# "+"   {{{
setMethod("+", signature(e1="FLStock", e2="BBMfit"),
          function(e1, e2) {
            if(validObject(e1) & validObject(e2)) {
              # si dim season 1 en el objeto e1, 
              # entonces cogemos solo lo del 1er season, info principio de anno, de e2
              # y si no, comprobar que los 2 objetos tengan las mismas seasons
              # XXX tambien hay que comprobar las iteraciones (seguramente mergeFLStock lo haga)
             
              # stoc info
              e2 <- FLStock( stock.n=stock(e2))
              stock.wt(e2) <- 1; units(stock.wt(e2)) <-''
              stock(e2) <- quantSums(stock.n(e2)*stock.wt(e2))
              
              # catch info
              catch.n(e2)[1,,,,,]  <- run@data$Crec
              catch.n(e2)[2,,,,,]  <- run@data$Cadu
              catch.wt(e2) <- 1; units(catch.wt(e2)) <-''
              catch(e2) <- quantSums(catch.n(e2)*catch.wt(e2))
              
              # harvest
              harvest(e2) <- catch.n(e2)/stock.n(e2)
              
              
              # recalcular harvest: captura/biomasa
              return(FLCore:::mergeFLStock(e1, e2))
            } else 
                stop("Input objects are not valid: validObject == FALSE")
          }
) # }}}

# residuals   {{{
setMethod("residuals", signature(object="BBMfit"),
          function(object) {

            stock.bio <- stock(object)
            pars      <- object@params
            
            nindex   <- object@data$nindex
            nyr      <- dim(stock.bio)[1]
            nper     <- dim(stock.bio)[4]-1
            indexper <- object@data$indexper
            
            # should we check that the nper from stock.bio and periods(indices) match??
            # if we apply this function to the class bbmFit I think it is not necessary
            
            residuals.Btot <- vector("list", length = nindex)
            residuals.P    <- vector("list", length = nindex)
            
            if(!is.null(names(run@data$indexper))){
              index.nam <- names(run@data$indexper)
            }else{
              index.nam <- paste("index", 1:nindex, sep="")
            }
            
            names(residuals.Btot) <- names(residuals.P) <- index.nam
            
            dnm <- dimnames(stock.bio)
            dnm$age    <- "all"
            dnm$season <- "1"
            
            for (i in 1:nindex){
              
              Btot <- Btot.obs <- P <- Pobs <- FLQuant( dimnames=dnm)
              
              B.yrs <- object@data$idxBobs[object@data$idxBobs[,1]==i,2]
              Btot              <- quantSums(stock.bio[,,,indexper[i],]) 
              Btot.obs[,B.yrs,] <- object@data$Bobs[object@data$idxBobs[,1]==i]
              
              P.yrs <- object@data$idxPobs[object@data$idxPobs[,1]==i,2]
              P             <- stock.bio[1,,,indexper[i],]/Btot[,,] 
              Pobs[,P.yrs,] <- object@data$Pobs[object@data$idxPobs[,1]==i]
              
              residuals.Btot[[i]] <- (log(Btot.obs) - pars@logq[i] - log(Btot)) * sqrt(exp(pars@logpsi[i]))
              residuals.P[[i]]    <- (Pobs - P) * sqrt( (1 + exp(pars@xi[i])) / ( P*(1-P) ))
              
            }
            
            return(list(residuals.Btot=residuals.Btot, residuals.P=residuals.P))
            
          }
) # }}}

# AIC	{{{
setMethod('AIC', signature(object='BBMfit'),
          function(object, ...){
            return(2*object@nopar - 2*object@nlogl)
          }
) # }}}

# BIC	{{{
setMethod('BIC', signature(object='BBMfit'),
          function(object, ...){
            return(object@nopar * log(object@nobs) - 2 * object@nlogl)
          }
) # }}}
