# bbmFit-class.R - Class ad methods for storing the results of running bbm()
# bbm/R/bbmFit-class.R

# Copyright European Union, 2018
# Authors:
#   Leire Ibaibarriaga <libaibarriaga@azti.es>
#   Sonia Sanchez (AZTI) <ssanchez@azti.es>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# bbmFit {{{

#' bbmFit class for storing the output of the bbm function.
#'
#' This includes abundance estimates in biomass (for recruits and adults) and
#' information on the model fitting.
#' 
#' @name bbmFit
#' @rdname bbmFit
#' @docType class
#' @aliases bbmFit bbmFit-methods bbmFit-class
#' @slot stock.bio Estimated stock biomass for recruits and adults in the different seasons, where seasons are dertermined by the index times. An *FLQuant* object with two age classes: "recruits" and "adults".
#' @slot params Estimated *bbm* parameters, *FLPar*.
#' @slot params.se Standard errors in parameters' estimates. *FLPar*.
#' @slot data Input data after processing, *list*.
#' @slot method Optimisation library and optimisation method, *TMB - method*.
#' @slot nopar Number of estimated parameters, *numeric*.
#' @slot noobs Number of observations, *numeric*.
#' @slot nlogl Negative log-likelihood, *numeric*.
#' @slot counts A two-element integer vector giving the number of calls to fn (i.e."bbm" function) and to gr (i.e. the optimisation method) respectively. This excludes those calls needed to compute the Hessian, if requested, and any calls to fn to compute a finite-difference approximation to the gradient.
#' @slot convergence An integer code. 0 indicates successful completion. For other possible error codes see \code{?optim}.
#' @slot message Character string giving any additional information returned by the optimizer.
#' @slot info Matrix with information on libraries' versions and platform used and object creation time.
#' @section Accessors:
#' All slots in the class have accessor and replacement methods defined that
#' allow retrieving and substituting individual slots.
#'
#' The values passed for replacement need to be of the class of that slot.
#' A numeric vector can also be used when replacing FLQuant slots, and the
#' vector will be used to substitute the values in the slot, but not its other
#' attributes.
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
#'     \item{bbmFit   }{: Creation of a new bbmFit object.}
#'     \item{stock    }{: Extracts information on estimated abundances in biomass for recruits and adults. Returns an \code{FLQuant}.}
#'     \item{+        }{: Updates an \code{FLStock} with new information on the BBM assessment.}
#'     \item{residuals}{: XXX.}
#'     \item{xxx}{: XXX.}
#'     \item{xxx}{: XXX.}
#'     \item{xxx}{: XXX.}
#' }
#' @author Leire Ibaibarriaga & Sonia Sanchez
#' @seealso \link{bbm}
#' @keywords classes
#' @examples
#' # Load data
#' data(ane)
#' 
#' # Run assessment (output is of class bbmFit)
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

setClass("bbmFit",
  representation(
    stock  = "FLQuant",
    params = "FLPar",
    params.se = "FLPar",
    data = "list",
    call = "character",
    logLik = "logLik",
    counts = "integer",
    convergence = "integer",
    message = "character",
    info = "matrix"),
  prototype=prototype(
    stock = FLQuant(NA, dimnames=list(age=c("recruits", "adults"))),
    params  = FLPar(logq.recruits=NA, logq.adults=NA, logpsi.recruits=NA,
      logpsi.adults=NA, xi.recruits=NA, xi.adults=NA,
      logB0=NA, mur=NA, logpsir=NA),
    params.se  = FLPar(logq.recruits=NA, logq.adults=NA, logpsi.recruits=NA,
      logpsi.adults=NA, xi.recruits=NA, xi.adults=NA,
      logB0=NA, mur=NA, logpsir=NA),
    data        = list(),
    call = character(1),
    logLik = "logLik",
    counts      = as.integer(NA),
    convergence = as.integer(NA),
    message     = "",
    info        = run.info(c("FLCore", "TMB", "bbm"))),
  validity=function(object){
    # params and parms.se have equal dimnames        

    # Everything is fine
    return(TRUE)
}) # }}}

# bbmFit()   {{{
setMethod('bbmFit', signature(object='missing'),
  function(object, dim="missing", dimnames="missing", ...) {
  
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

    res <- new("bbmFit", stock.bio = flq)

    # Load given slots
    if(length(args) > 0)
      for(i in names(args))
        slot(res, i) <- args[[i]]
            
    return(res)
  }
) # }}}

# accessors {{{

# stock
setMethod("stock", signature(object="bbmFit"),
  function(object, ...) {
    return(object@stock)
  } 
) 

setReplaceMethod("stock", signature(object="bbmFit", value="FLQuant"),
  function(object, value) {
    object@stock <- value
    return(object)
  } 
) 

# params
setMethod("params", signature(object="bbmFit"),
  function(object, ...) {
    return(object@params)
  } 
) 

setReplaceMethod("params", signature(object="bbmFit", value="FLPar"),
  function(object, value) {
    object@params <- value
    return(object)
  } 
)

# params.se
setMethod("params.se", signature(object="bbmFit"),
  function(object, ...) {
    return(object@params.se)
  } 
) 

setReplaceMethod("params.se", signature(object="bbmFit", value="FLPar"),
  function(object, value) {
    object@params.se <- value
    return(object)
  } 
)

# }}}

# "+"   {{{
setMethod("+", signature(e1="FLStock", e2="bbmFit"),
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
setMethod("residuals", signature(object="bbmFit"),
  function(object, type=c("biomass", "recruits")) {

  # EXTRACT elements
  stock <- stock(object)
  params <- object@params
  nindex   <- object@data$nindex
  nyr      <- dim(stock.bio)[1]
  nper     <- dim(stock.bio)[4] - 1
  indexper <- object@data$indexper
  


  residuals.Btot <- vector("list", length = nindex)
  residuals.P    <- vector("list", length = nindex)
 
  # SET names  
  if(!is.null(names(run@data$indexper))) {
    index.nam <- names(run@data$indexper)
  } else {
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
  
  residuals.Btot[[i]] <- (log(Btot.obs) - pars@logq[i] - log(Btot)) *
    sqrt(exp(pars@logpsi[i]))
  residuals.P[[i]]    <- (Pobs - P) * sqrt( (1 + exp(pars@xi[i])) / ( P*(1-P) ))
  }
   
            return(list(residuals.Btot=residuals.Btot, residuals.P=residuals.P))
            
          }
) # }}}

# AIC	{{{
setMethod('AIC', signature(object='bbmFit'),
          function(object, ...){
            return(2*object@nopar - 2*object@nlogl)
          }
) # }}}

# BIC	{{{
setMethod('BIC', signature(object='bbmFit'),
          function(object, ...){
            return(object@nopar * log(object@nobs) - 2 * object@nlogl)
          }
) # }}}
