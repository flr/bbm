#-------------------------------------------------------------------------------  
# bbmFit class.
# Created: Sonia Sanchez - 2018-05-22 13:49:23
# Changed: 2018-06-01 10:24:36 (ssanchez)
#------------------------------------------------------------------------------- 

# bbmFit_class.R - output of bbm function
# bbm/R/bbmFit_class.R

# Copyright: European Union, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# bbmFit

#' @title bbmFit class
#' 
#' @description \code{bbmFit} class is used for storing the output of the \code{bbm} function. 
#'              This includes abundance estimates in biomass (for recruits and adults) and information on the model fit.
#' 
#' @name bbmFit
#' @rdname bbmFit
#' @docType class
#' @aliases bbmFit bbmFit-class input,bbmFit-method convergence,bbmFit-method message,bbmFit-method fitSumm,bbmFit-method
#'          params,bbmFit-method params.se,bbmFit-method vcov,bbmFit-method stock.bio,bbmFit-method indicesB,bbmFit-method 
#'          indicesP,bbmFit-method
#'          plot,bbmFit-method
#'
#' @section Slots:
#'  \describe{
#'     \item{input      }{ Input data. \code{list}, Containing the following information: 
#'                          catch, indicesB, indicesP, perindicesB, perindicesP, control, f and nper.}
#'     \item{convergence}{ Convergence code, \code{vector(niter)}. Where 0 indicates successful completion. 
#'                          For other possible error codes see \code{?optim}.}
#'     \item{message    }{ Character string giving any additional information returned by the optimizer, or \code{""}.}
#'     \item{fitSumm    }{ Fit summary (with information on 'nlogL', 'nobs', 'nopar'), \code{array[3,niter]}.}
#'     \item{params     }{ Estimated parameters in \code{bbm} function, \code{FLPar[npar,niter]}, in linear scale.}
#'     \item{params.se  }{ Standard errors in parameters' estimates. \code{FLPar[npar,niter]}, in linear scale.}
#'     \item{vcov       }{ Variance-covariance matrix, \code{array[npar,npar,niter]}.}
#'     \item{stock.bio  }{ Estimated stock biomass for recruits and adults in the different seasons, 
#'                          where seasons are dertermined by the index times. 
#'                          \code{FLQuant} with two age classes: recruits and adults. }
#'     \item{indicesB   }{ Estimates of surveys' total abundances in biomass, \code{FLQuants}.}
#'     \item{indicesP   }{ Estimates of surveys' percentage of recruits in biomass, \code{FLQuants}.}
#' }
#'
#'
#' @section Validity: 
#' \describe{
#'     \item{Dimensions}{ \itemize{ \item{\code{age}:} 
#'                                       {stock.bio must be an \code{FLQuant} with only 2 age classes (recruits and adults) and 
#'                                        each index in \code{indicesB} must be an \code{FLQuant} with only 1 age class ('all')}
#'                                  \item{\code{year, unit, season, area}:} 
#'                                       {equal for \code{stock.bio}, \code{indicesB} and \code{indicesP}}
#'                                  \item{\code{iter}:} 
#'                                       {equal for \code{stock.bio}, \code{convergence}, \code{fitSumm}, 
#'                                        \code{indicesB} and \code{indicesP}}
#'                                } 
#'                      }
#'      \item{Parameters}{Same number of parameters required in \code{params}, \code{params.se} and \code{vcov}}
#' }
#' You can inspect the class validity function by using
#'   \code{getValidity(getClassDef('bbmFit'))}
#'   
#' @section Accessors:
#' All slots in the class have accessor methods defined that allow retrieving individual slots.
#' 
#' @section Constructor:
#' A construction method exists for this class that can take named arguments for
#' any of its slots. All slots are then created to match the requirements of the
#' class validity. If \code{years}, \code{niter}, \code{namesB} or \code{namesP} are provided, 
#' this is used for sizing and naming the different slots.
#'
#' @section Methods:
#' Methods exist for various calculations based on values stored in the class:
#'
#' \describe{
#'     \item{+        }{ Updates an \code{FLStock} with new information on the BBM assessment.}
#'     \item{residuals}{ Calculates Pearson residuals, returns an object of class \code{bbmFitresiduals}.}
#'     \item{logLik   }{ Method to extract Log-Likelihood, returns an object of class \code{logLik}.}
#'     \item{AIC      }{ Method to calculate Akaike's 'An Information Criterion' (AIC) of a \code{bbmFit} object 
#'                       from the value of the obtained log-likelihood stored in its \code{logLik} slot.}
#'     \item{BIC      }{ Method to calculate the Bayesian information criterion (BIC), 
#'                       also known as Schwarz's Bayesian criterion of a \code{bbmFit} object 
#'                       from the value of the obtained log-likelihood stored in its \code{logLik} slot.}
#'     \item{iter     }{ Extracts a subset of the iterations contained in a \code{bbmFit} object.}
#'     \item{plot     }{ One plot for estimated abundances and one extra plot for each of the surveys 
#'                       with the fitting of total biomass and proportion of recruits.}
#' }
#' 
#' 
#' @author Leire Ibaibarriaga & Sonia Sanchez
#' @seealso \link{bbm}, \link{bbmFitresiduals}, \link{logLik}, \link{bbmFLPar}, \link{plot}
#' @keywords classes
#' 
#' @examples
#' 
#' # Load data
#' data(ane)
#' 
#' # Generate an object of bbmFit class (different alternatives)
#' new("bbmFit")                # empty object
#' slotNames(bbmFit())          # slots
#' 
#' # bbmFit: setting dimensions for stock.bio
#' bbmFit( stock.bio = FLQuant(dim=c(2,20,1,3,1,1), dimnames=list(age=1:2, year=1980:1999))) 
#' 
#' # bbmFit: params class - FLPar with specific parameters for bbm function
#' bbmFit( params=bbmFLPar(years=dimnames(catch.ane)$year, 
#'                         namesB=names(indicesB.ane), namesP=names(indicesP.ane))) 
#' 
#' # Run assessment (output is of class bbmFit)
#' run <- bbm(catch.ane, indicesB=indicesB.ane, indicesP=indicesP.ane, control=control.ane, inits=inits.ane)
#' class(run)
#' run
#' 
#' # Plot
#' plot(run)
#' 

setClass("bbmFit",
         representation(
           inputs      = "list",
           convergence = "numeric",
           message     = "character",
           fitSumm     = "array",
           params      = "FLPar",
           params.se   = "FLPar",
           vcov        = "array",
           stock.bio   = "FLQuant",
           indicesB    = "FLQuants",
           indicesP    = "FLQuants"
          ),
         prototype=prototype(
           inputs      = list(),
           convergence = new('numeric'),
           message     = "",
           fitSumm     = new('array'),
           params      = new('FLPar'),
           params.se   = new('FLPar'),
           vcov        = new('array'),
           stock.bio   = new('FLQuant'),
           indicesB    = new('FLQuants'),
           indicesP    = new('FLQuants')
          ),
         validity=function(object){
           
           # All FLQuant objects must have same dimensions
           # (except for ages, where stock.bio has two age classes, recruits and adults, 
           # and indicesB and indicesP only one)
           
           if (dim(object@stock.bio)[1]!=2)
             stop("'stock.bio' must be an FLQuant with only 2 age classes (recruits and adults).")
           
           if ( sum(unlist(lapply( object@indicesB, function(x) dim(x)[1])) != 1)>0 | 
                (sum(unlist(lapply( object@indicesB, function(x) dim(x)[1])) != 1)==0 & sum(unlist(lapply( object@indicesB, function(x) dimnames(x)[1])) !="all")>0) )
             stop("Each index in 'indicesB' slot must be an FLQuant with only 1 age class ('all'), as should give total biomass estimated by the survey.")
           
           if ( sum(unlist(lapply( object@indicesP, function(x) dim(x)[1])) != 1)>0 | 
                (sum(unlist(lapply( object@indicesP, function(x) dim(x)[1])) != 1)==0 & sum(unlist(lapply( object@indicesP, function(x) dimnames(x)[1])) !="all")>0) )
             stop("Each index in 'indicesP' slot must be an FLQuant with only 1 age class ('all'), as should give proportion of recruits in biomass estimated by the survey.")
           
           Dim <- dim(object@stock.bio)[-1]
           if ((sum(dim(object@indicesB)[-1] != Dim) + sum(dim(object@indicesP)[-1] != Dim))>=1)
             return("stock.bio, indicesB and indicesP slots must have same dimensions except the first one (age)")
           
           
           # Same number of iterations for the rest of the slots 
           niter <- dim(object@stock.bio)[6]
           
           if (length(object@convergence) != niter)
             stop("convergence slot must have length equal to the number of iterations in stock.bio.")
           
           if (dim(object@fitSumm)[2] != niter)
             stop("fitSumm slot must have length equal to the number of iterations in stock.bio.")
           
           # same number of parameters
           npar <- dim(object@params)[1]
           
           if (dim(object@params.se)[1] != npar | dim(object@vcov)[1] != npar)
             stop("Different number of parameters in slots params, params.se and vcov.")
           
           
           # Check the names of the parameters
           
           nam <- unlist(getNamesFromParnames(dimnames(object@params)$params))
           if ( sum(unlist(getNamesFromParnames(dimnames(object@params.se)$params)) != nam)>0 | 
                sum(unlist(getNamesFromParnames(dimnames(object@vcov)[[1]],logscale=TRUE)) != nam)>0 )
             stop("Parameter names in params, params.se and vcov are unconsistent.")
           
           
           # Everything is fine
           return(TRUE)
          
          }
)


#---------------------
#---------------------

# -- CREATOR

# bbmFit()   {{{
#' @rdname bbmFit
#' @aliases bbmFit,missing-method
setMethod('bbmFit', signature(object='missing'),
          function(object, years="missing", niter="missing", namesB="missing", namesP="missing", ...)
          {
            
            args <- list(...)
            
            # Create FLQuants with correct dimensions
            
            if (missing(niter)) niter <- 1
            
            if (missing(years)) {
              flqa <- FLQuant( dimnames=list(age=1:2), iter=niter)
            } else {
              flqa <- FLQuant( dimnames=list(age=1:2, year=years), iter=niter)
            }
            
            res <- new("bbmFit")
            
            # Load given slots
            
            if(length(args) > 0)
              for(i in names(args))
                slot(res, i) <- args[[i]]
            
            # Resize slots with correct dimensions (if not given)
            
            if (is.null(args$convergence)) 
              res@convergence <- numeric(niter)*NA
            
            if (is.null(args$fitSumm))
              res@fitSumm     <- array(dim = c(3,niter), dimnames = list(c("nlogl","nobs","nopar"),iter=1:niter))
            
            if (missing(namesB) | missing(namesP) | missing(years))
              # get names from any given slot with this info
              if (is.null(args$params) & is.null(args$params.se) & is.null(args$vcov)) {
                if (missing(namesB)) namesB <- "index1"
                if (missing(namesP)) namesP <- "index1"
                if (missing(years))  years  <- 1
              } else {
                if (!is.null(args$params)) {
                  nam <- getNamesFromParnames(dimnames(res@params)$params)
                } else if (!is.null(args$params.se)) {
                  nam <- getNamesFromParnames(dimnames(res@params.se)$params)
                } else {
                  nam <- getNamesFromParnames(dimnames(res@params.se)$params)
                }
                if (missing(namesB)) namesB <- nam$namesB
                if (missing(namesP)) namesP <- nam$namesP
                if (missing(years))  years  <- nam$years
              }

            if (is.null(args$params))
              res@params      <- bbmFLPar( years=years, namesB=namesB, namesP=namesP, niter=niter)

            if (is.null(args$params.se))
              res@params.se   <- bbmFLPar( years=years, namesB=namesB, namesP=namesP, niter=niter)
            
            if (is.null(args$vcov)){
              params <- c(paste('logq',namesB,sep='_'),paste('logpsi',namesB,sep='_'),
                          paste('xi',namesP,sep='_'),'logB0',
                          paste('logR',years,sep='_'),'mur','logpsir')
              res@vcov <- array( dim = c(length(params), length(params), niter), dimnames = list( params, params) )
            }
              
            if (is.null(args$stock.bio))
              res@stock.bio   <- flqa
              
            if (is.null(args$indicesB)){
              flq <- quantSums(flqa)
              for(i in 1:length(namesB))
                res@indicesB[[i]] <- flq
              names(res@indicesB) <- namesB
            }

            if (is.null(args$indicesP)) {
              flq <- quantSums(flqa)
              for(i in 1:length(namesP))
                res@indicesP[[i]] <- flq
              names(res@indicesP) <- namesP
            } 
            
            # check object validity
            validObject(res)
            
            return(res)
            
          }
) # }}}


#---------------------
#---------------------

getNamesFromParnames <- function(parnam, logscale=FALSE) { # function to get namesB, namesP and years from parameter names
  
  namesP <- unlist(lapply( strsplit(parnam[grep("^xi_",parnam)],split="_"), function(x) x[2]))
  
  if (logscale==FALSE) {
    namesB <- unlist(lapply( strsplit(parnam[grep("^q_",parnam)],split="_"), function(x) x[2]))
    years  <- unlist(lapply( strsplit(parnam[grep("^R_",parnam)],split="_"), function(x) x[2]))
  } else {
    namesB <- unlist(lapply( strsplit(parnam[grep("^logq_",parnam)],split="_"), function(x) x[2]))
    years  <- unlist(lapply( strsplit(parnam[grep("^logR_",parnam)],split="_"), function(x) x[2]))
  }
  
  return(list(namesB=namesB, namesP=namesP, years=years))
}

#---------------------
#---------------------


# -- ACCESSORS {{{

#' @rdname bbmFit-class
setMethod("inputs", "bbmFit", function(object) object@inputs)

#' @rdname bbmFit-class
setMethod("convergence", "bbmFit", function(object)  object@convergence)

#' @rdname bbmFit-class
setMethod("message", "bbmFit", function(object)  object@message)

#' @rdname bbmFit-class
setMethod("fitSumm", "bbmFit", function(object)  object@fitSumm)

#' @rdname bbmFit-class
setMethod("params", "bbmFit", function(object)  object@params)

#' @rdname bbmFit-class
setMethod("params.se", "bbmFit", function(object)  object@params.se)

#' @rdname bbmFit-class
setMethod("vcov", "bbmFit", function(object)  object@vcov)

#' @rdname bbmFit-class
setMethod("stock.bio", "bbmFit", function(object)  object@stock.bio)

#' @rdname bbmFit-class
setMethod("indicesB", "bbmFit", function(object)  object@indicesB)

#' @rdname bbmFit-class
setMethod("indicesP", "bbmFit", function(object)  object@indicesP)

# }}}

#---------------------
#---------------------

# -- METHODS

#' @rdname bbmFit
#' @aliases +,FLStock,bbmFit-method
#' @examples
#'
#' stock <- FLStock(catch.n=catch.ane, catch.wt=catch.ane*0+1)
#' units(stock@catch.wt) <- ''
#' stock@catch <- quantSums(stock@catch.n*stock@catch.wt)
#' 
#' newst <- stock + run # we must sum to the bbmFit object not to stock.bio(run)
#' 

# "+"   {{{
setMethod("+", signature(e1="FLStock", e2="bbmFit"),
          function(e1, e2) {
            if(validObject(e1) & validObject(e2)) {
              
              #! At the moment only valid if both objects have exactly the same dimensions
              #  It would be possible to do the following: 
              #  - if only one season in e1 --> take 1st season in e2
              #  - else check that both objects have the same number of seasons
              
              # stock info
              e2 <- FLStock( stock.n=stock.bio(e2))
              stock.wt(e2) <- 1; units(stock.wt(e2)) <-''
              stock(e2) <- quantSums(stock.n(e2)*stock.wt(e2))
              
              # catch info
              catch.n(e2) <- run@inputs$catch
              catch.wt(e2) <- 1; units(catch.wt(e2)) <-''
              catch(e2) <- quantSums(catch.n(e2)*catch.wt(e2))
              
              # harvest
              harvest(e2) <- catch.n(e2)/stock.n(e2)
              
              
              # recalcular harvest: captura/biomasa
              return(FLCore:::mergeFLStock(e1, e2)) #! XXX check if this is the appropriate function
            } else 
                stop("Input objects are not valid: validObject == FALSE")
          }
) # }}}


#' @rdname bbmFit
#' @aliases residuals,bbmFit-method
#' @examples
#'
#' # calculate residuals
#' residuals(run)

# residuals   {{{
setMethod("residuals", signature(object="bbmFit"),
          function(object) {

            nindicesB <- length(object@indicesB)
            nindicesP <- length(object@indicesP)
            
            residuals.B <- residuals.P <- FLQuants()
            # residuals.B <- vector("list", length = nindicesB)
            # names(residuals.B) <- names(object@indicesB)
            # residuals.P <- vector("list", length = nindicesP)
            # names(residuals.P) <- names(object@indicesP)

            for (i in 1:nindicesB){
              residuals.B[[i]] <- (log(object@inputs$indicesB[[i]]) - log(object@indicesB[[i]])) * sqrt( object@params[ grep("^psi_", rownames(object@params))[i], ])
            }
            names(residuals.B) <- names(object@indicesB)
                        
            for (i in 1:nindicesP){
              residuals.P[[i]] <- (object@inputs$indicesP[[i]] - object@indicesP[[i]]) * 
                                  sqrt( (1 + exp(object@params[ grep("^xi_", rownames(object@params))[i], ])) / ( object@indicesP[[i]]*(1-object@indicesP[[i]]) ))
            }
            names(residuals.P) <- names(object@indicesP)
            
            res <- new("bbmFitresiduals")
            res@residuals.B <- residuals.B
            res@residuals.P <- residuals.P
            
            return(res)
            
          }
) # }}}

#' @rdname bbmFit
#' @aliases logLik,bbmFit-method
#' @examples
#'
#' # log-Likelihood
#' logLik(run)

setMethod("logLik", signature(object = 'bbmFit'),
          function(object, ...)
          {
            dim2 <- length(dim(object @ fitSumm))
            if (dim2 == 1) {
              val <- -1 * unname(object @ fitSumm["nlogl"])
              attr(val, "nobs") <- unname(object @ fitSumm["nobs"])
              attr(val, "df") <- unname(object@fitSumm["nopar"])
            } else if (dim2 == 2) {
              val <- -1 * unname(object @ fitSumm["nlogl",])
              attr(val, "nobs") <- unname(object @ fitSumm["nobs",])
              attr(val, "df") <- unname(object@fitSumm["nopar",])
            }
            class(val) <- "logLik"
            val
          }) #! Copied from a4aFit-class.R, but strange output when more than 1 iteration


#' @rdname bbmFit
#' @aliases AIC,bbmFit-method
#' @examples
#'
#' # AIC and BIC
#' AIC(run)

# AIC	{{{
setMethod('AIC', signature(object='bbmFit'),
          function(object, ...){
            return(AIC(logLik(object))) 
            # return(2*object@fitSumm["nopar",] + 2*object@fitSumm["nlogl",])
          }
) # }}}

#' @rdname bbmFit
#' @aliases BIC,bbmFit-method
#' @examples
#'
#' BIC(run)

# BIC	{{{
setMethod('BIC', signature(object='bbmFit'),
          function(object, ...){
            return(BIC(logLik(object))) 
            # return(object@fitSumm["nopar",]*log(object@fitSumm["nobs",]) + 2*object@fitSumm["nlogl",])
          }
) # }}}

#' @rdname bbmFit
#' @aliases iter,bbmFit-method
#' @param obj The object to be subset
#' @param it  Iteration(s) to be extracted

setMethod("iter", "bbmFit", function(obj, it){
  
  obj@inputs$catch <- iter(obj@inputs$catch, it)
  obj@inputs$indicesB <- iter(obj@inputs$indicesB, it)
  obj@inputs$indicesP <- iter(obj@inputs$indicesP, it)
  obj@inputs$control@param.fix <- iter(obj@inputs$control@param.fix, it)
  
  obj@convergence <- obj@convergence[it]
  obj@fitSumm <- obj@fitSumm[,it, drop=FALSE]
  obj@params    <- iter(obj@params, it)
  obj@params.se <- iter(obj@params.se, it)
  obj@vcov <- obj@vcov[,,it]
  obj@stock.bio <- iter(obj@stock.bio, it)
  obj@indicesB  <- iter(obj@indicesB, it)
  obj@indicesP  <- iter(obj@indicesP, it)

  obj
  
})


# # sd {{{
# setMethod('sd', signature(x='bbmFit', na.rm='missing'),
#           function(x)
#           {
#             if(dim(params(x))[2] == 1)
#             {
#               res <- as.vector(t(sqrt(apply(x@vcov, 3, diag))))
#               names(res) <- dimnames(x@vcov)[[1]] #dimnames(params(x))$params
#               return(res)
#             }
#             else
#               return(sd(params(x)))
#           }
# ) # }}}
# 
# logxse <- res/sqrt(x@fitSumm['nobs',])
# xse    <- exp(logxse) - 1
# x@params.se
# 
# # cv {{{
# setMethod('cv', signature(x='bbmFit'),
#           function(x, ...)
#           {
#             if(dim(params(x))[1] == 1)
#               return(sd(x)/apply(params(x), 2, mean))
#             else
#               return(sd(params(x))/mean(params(x)))
#           }
# ) # }}}

#' @name plot
#' @docType methods
#' @rdname plot-methods
#' @aliases plot plot,bbmFit,missing-method
#' @description Method to plot bbm fitting, if applied to a \code{bbmFit} object. The output are several plots. 
#'              Firstly, the estimated abundance and next the fitted versus observed indices (with one plot for each survey). 
#'              Note that in plots related to indices the yaxis doesn't has a scale. 
#'              The visual is about the difference between the two lines, not about the value of each line, which in any case would be very difficult to assess visually.
#' @param x   A \code{bbmFit} object with the fitted values and observed values.
#' @param y   Missing.
#' @param ... Additional argument list.
#' @return If \code{class(x)=='bbmFit'}, a \code{plot} with estimated abundances and one extra \code{plot} for each survey with fitted and observed indices.
#' 
#' @examples
#' 
#' # Load data
#' data(ane)
#' 
#' # Run the assessment
#' obj <- bbm( catch.ane, indicesB=indicesB.ane, indicesP=indicesP.ane, control=control.ane, inits=inits.ane)
#' 
#' # Plot the output
#' plot(obj)
#' 

setMethod("plot", c("bbmFit", "missing"), function(x, y, ...){
  
  # Plot estimated abundance
  
  biomass <- stock.bio(x)
  dimnames(biomass)$season <- paste( "season", dimnames(biomass)$season, sep=" ")
  
  bioplot <- ggplot(biomass, mapping=aes(year, data)) + 
    geom_line(aes(group=age, colour=factor(age))) + 
    facet_grid(~season) + 
    labs(x="", y=paste( "biomass (", units(biomass), ")", sep="")) + 
    scale_colour_discrete(name = "age")
  
  print(bioplot)
  
  
  # Plot indices fit
  
  indicesBfit <- x@indicesB
  indicesPfit <- x@indicesP
  indicesBobs <- x@inputs$indicesB
  indicesPobs <- x@inputs$indicesP
  dfx <- rbind(cbind(as.data.frame(indicesBfit), idx = "total.biomass"), 
               cbind(as.data.frame(indicesPfit), idx = "prop.recruits"))
  dfy <- rbind(cbind(as.data.frame(indicesBobs), idx = "total.biomass"), 
               cbind(as.data.frame(indicesPobs), idx = "prop.recruits"))
  dfx$src <- "fit"
  dfy$src <- "obs"
  df0 <- rbind(dfx, dfy)
  df0$src <- factor(df0$src, levels = c("obs", "fit"))
  df0$idx <- factor(df0$idx, levels = c("total.biomass", "prop.recruits"))
  args <- list()
  args$x <- data ~ factor(year)|factor(idx)
  args$type <- c("l")
  args$groups <- quote(src)
  args$ylab <- ""
  args$xlab <- ""
  args$scales <- list(y = list(relation = "free", draw=FALSE))
  args$auto.key <- list(lines = TRUE, points = FALSE, columns = 2)
  args$par.settings <- list(superpose.line = list(col = c("gray70", 
                                                          "black"), lty = 1, lwd = c(2, 1)), strip.background = list(col = "gray90"), 
                            strip.border = list(col = "black"))
  if (length(x@indicesB) > 1) {
    for (i in unique(names(x@indicesB), names(x@indicesP))) {
      args$data <- subset(df0, qname == i)
      args$main <- paste(i, " fitted and observed indices", 
                         sep = ":")
      print(do.call("xyplot", args))
    }
  }
  else {
    args$data <- df0
    args$main <- "fitted and observed indices"
    do.call("xyplot", args)
  }
  
})

