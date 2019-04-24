# BBMass.R - function to simulate a two-stage biomass-based model
# bbm/R/BBMass.R

# Copyright: European Union & AZTI, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# BBM {{{

#' Function to simulate a two-stage biomass-based model.
#' 
#' BBM function simulates a two-stage biomass model. Specifically a maximum likelihood version of a generalisation of the model in Ibaibarriaga et al. (2008).
#' In Ibaibarriaga et al. (2008) only two indices ocurring at the same moment of the year were considered. 
#' However, in this case is generalised to n indices that can occur in different times of the year.
#'
#' @name BBM
#' @rdname BBM
#' @aliases BBM BBM-methods
#'
#' @param object  An \code{FLQuant} with catch information or an \code{FLStock}. These objects must only have two age classes (recruits and adults) 
#'                and the number of seasons should be 1 or the number of seasons determined by the period of the year of the different indices.
#' @param indices Abundance indices in biomass for same age classes (\code{FLIndices}).
#' @param control An \code{BBM.control} with control arguments.
#' @param inits   An \code{BBMpar} with initial values.
#' 
#' 
#' @return An object of class BBMfit.
#' 
#' @author Leire Ibaibarriaga & Sonia Sanchez.
#' @seealso \link{BBMfit}.
#' @keywords BBM
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
#' # Run assessment
#' run <- BBM(object=catch.ane, indices=indices.ane, control=control.ane, inits=inits.ane)
#' 
#' is(run)
#' names(run)
#' 
#' # Plot assesses population
#' \dontrun{
#'   plot(run@stock.bio)
#' }
#' 

load("../data/ane.RData")
catch <- catch.ane
indices <- lapply(indices.ane, index)
indices <- indices.ane
control <- control.ane
inits <- inits.ane


setMethod("bbm", signature(object="FLQuant", indices="FLIndices"),
  function(catch, indices, control, inits) {

  # DATA CHECKS
  # check <- BBM.checks( object, indices, control, inits)
  
  # TODO Check that there are no NA's in catch
  if (any(is.na(catch))){
    stop("NA's not allowed in catch")
  }
  
  # SET timing
  aux <- periods(indices)
  nper <- aux$nper
  f    <- aux$f
  indexper <- aux$indexper
  
  # DIMENSIONS and NAMES
  dm <- dim(catch)
  yrs <- dimnames(catch)$year
  idxs   <- names(indices)
  nindex <- length(indices)
  
  # EXTEND catch to no. seasons
  if (dm[4] == 1 & nper > 1)  {
    catch <- expand(catch, season=1:nper)
    catch[,,,nper:1,,] <- catch[,,,1,,] * f[s]
  }

  # - control parameters
  g <- control@g
  param.fix <- control@param.fix
  
  # DATA transformation for input
  Crec <- catch[1, seq(1, dm[2]),, 1:nper,,, drop=T]  # recruits catch
  Cadu <- catch[2, seq(1, dm[2]),, 1:nper,,, drop=T]  # adults catch

  Bobs    <- NULL
  idxBobs <- NULL
  Pobs    <- NULL
  idxPobs <- NULL

  for (i in 1:nindex){
    aux <- quantSums(indices[[i]]@index)[,,,,,,drop=T]
    aux <- aux[!is.na(aux)]
    Bobs <- c(Bobs, aux)
    idxBobs <- rbind( idxBobs, cbind( rep(i, length(aux)), match(names(aux), yrs)))
    aux <- indices[[i]]@index[1,,,,,]/quantSums(indices[[i]]@index)
    aux <- aux[,,,,,,drop=T]
    aux <- aux[!is.na(aux)]
    Pobs <- c(Pobs, aux)
    idxPobs <- rbind( idxPobs, cbind( rep(i, length(aux)), match(names(aux), yrs)))
  }
  
  # create index
  dat <- list(grec=g['rec'], gadu=g['adult'], 
    nindex=nindex, 
    Crec=Crec, Cadu=Cadu,
    Bobs=Bobs, Pobs=Pobs,
    idxBobs=idxBobs, idxPobs=idxPobs,
    indexper = indexper,
    f = as.double(f))
  
  # INITIAL PARAMETERS
  inits <- as.list(inits)
  
  # MODEL
  
  # model definition
  if (sum(unlist(as.list(param.fix))) == 0) {
    model <- MakeADFun(dat, inits, DLL="bbm", silent=TRUE)
  } else {
    map <- list()
    for (sl in slotNames(param.fix)[!slotNames(param.fix) %in% ".Data"]) {
      if(sum(slot( param.fix, sl))>0) 
        map[[sl]] <- slot( param.fix, sl)
    }
    map <- lapply( map, function(x) {
      x[x==1] <- NA
      return(as.factor(x))
    })
    model <- MakeADFun(dat, inits, DLL="bbm", map=map)
  }

  model$hessian <- TRUE
  
  # model fit: minimise the negative log-likelihood (nll)
  fit <- optim(model$par, model$fn, model$gr, control=list(maxit=1e+9))
  
  rep   <- summary(sdreport(model))
  param <- rownames(rep)
  nopar <- length(param)
  rownames(rep) <- NULL
  colnames(rep) <- c("Est","SE")
  
  rep <- as.data.frame(rep)
  rep$param <- param
  
  # parameter estimates
  
  browser()
  
  for (sl in slotNames(param.fix)[!slotNames(param.fix) %in% ".Data"]) {
    slot(params, sl) <- slot( param.fix, sl)
    if(sum(slot( params, sl))>0){
      slot(params, sl)[slot( params, sl)==0] <- rep$Est[rep$param==sl]
    }else{
      slot(params,sl) <- rep$Est[rep$param==sl]
    } 
  }

  # SE of parameter estimates
  
  params.se <- new("BBMpar")
  
  for (sl in slotNames(param.fix)[!slotNames(param.fix) %in% ".Data"]) {
    slot(params.se, sl) <- slot( param.fix, sl)
    if(sum(slot( params.se, sl))>0){
      slot(params.se, sl)[slot( params.se, sl)==0] <- rep$SE[rep$param==sl]
      slot(params.se, sl)[slot( params.se, sl)==1] <- 0.0
    }else{
      slot(params.se,sl) <- rep$SE[rep$param==sl]
    } 
  }
  
  # compute biomasses (this part can be moved to a function)
  stock.fit <- calcPop( g=g, f=f, catch=object, inits=params)
  
  if (stock.fit$ok == FALSE)
    stop("Returned parameter estimates lead to negative biomassess.")
      
  stock.bio <- stock.fit$stock
  
  
  # OUTPUT
  
  out <- bbmFit(stock.bio=stock.bio) # estimated biomass
  
  # - parameters and their se
  out@params       <- params
  out@params.se    <- params.se
  
  # - input data
  out@data        <- dat
  
  # - optimisation method
  out@method       <- paste("TMB", model$method, sep=" - ")
  
  # - fitting summary (number of parameters, negative log likelihood, )
  out@nopar       <- nopar
  out@nobs        <- nrow(idxBobs) + nrow(idxPobs)
  out@nlogl       <- fit$value
  out@counts      <- fit$counts
  out@convergence <- fit$convergence
  out@message     <- ifelse( is.null(fit$message), "", fit$message)
  
  #Populate the info slot as last step
  info <- data.frame(FLbbm.version=packageDescription("bbm")$Version,
                     FLCore.version=packageDescription("FLCore")$Version,
                     TMB.version=packageDescription("TMB")$Version,
                     R.version=R.version$version.string,
                     platform=R.version$platform,
                     run.date=Sys.time())
  out@info <- t(info)
  colnames(out@info) <- ""
    
  return(out)

})  # }}}

# BBM.checks {{{
BBM.checks <- function( object, indices, control=NULL, inits=NULL) {
  
  # Checking that argument's class is correct
  
  if (!class(object) %in% c('FLQuant'))
    stop( "'object' argument must be of class 'FLQuant' or 'FLStock'.")
  
  if (!class(indices) %in% c('FLIndices'))
    stop( "'indices' argument must be an 'FLIndices' object.")
  
  if (!is.null(control) & class(control)!='BBM.control')
    stop( "'control' argument must be of class 'BBM.control'.")
  
  if (!is.null(inits) & class(inits)!='BBMpar')
    stop( "'inits' argument must be of class 'BBMpar'.")
  
  
  # Age dimension
  
  if (dim(object)[1]!=2)
    stop("'object' argument must have only 2 age classes (one for recruits and the other for adults).")
  
  if ( sum(unlist(lapply( indices, function(x) dim(x)[1]) !=2))>0)
    stop("Each index in 'indices' argument must have only 2 age classes (one for recruits and the other for adults).")
  
  # Unit, area and iter dimensions
  
  if (sum(dim(object)[c(3,5,6)]>1)>0)
    warning("'object' should have length 1 for dimensions unit, area and iter. Therefore the 1st value will be used.")
  if (sum(unlist(lapply( indices, function(x) dim(x@index)[c(3,5,6)]))>1)>0)
    warning("'indices' should have length 1 for dimensions unit, area and iter. Therefore the 1st value will be used.")
  
  
  # Season dimension
  
  aux <- periods(indices)
  nper <- aux$nper
  # f    <- aux$f
  # indexper <- aux$indexper
  
  if (dim(object)[4]!=1 & dim(object)[4]!=nper)
    stop(paste("'object' argument must have 1 or ", nper, "seasons (being this last value determined by the period of the year of the different indices)."))
  
  if (sum(unlist(lapply( indices, function(x) dim(x)[4]))!=1)>0)
    stop("Indices must have only one season.")
  
  
  # Year dimension
  
  nyear <- dim(object)[2]
  # (indices and object should have same range of years, otherwise add NA's in the indices object)
  if ( !AllEqual(unlist(lapply( indices, function(x) dims(x)$minyear))) | 
       !AllEqual(unlist(lapply( indices, function(x) dims(x)$minyear))) |
       dims(object)$minyear != min(unlist(lapply( indices, function(x) dims(x)$minyear))) |
       dims(object)$maxyear < max(unlist(lapply( indices, function(x) dims(x)$maxyear))) )
    stop("Check year ranges in the input objects, as years must be equal in 'object' and 'indices'.")
  
  
  # check units
  
  ui <- lapply( indices, function(x) units(x)$index )
  if (any(units(object)!=ui))
    stop("Units from 'object' and 'indices' must be equal.")
  
  #! check dimensions in param.fix
  if (!is.null(control)) {
    nindex <- length(indices)
    if (length(control@param.fix@logq) != nindex | length(control@param.fix@logpsi) != nindex | length(control@param.fix@xi) != nindex)
      stop("'logq', 'logpsi' and 'xi' in 'control@param.fix' must have one value for each index in indices.")
    if (length(control@param.fix@logR) != nyear)
      stop(paste("'control@param.fix@logR' must be of length ", nyear, ", with one value for each year.", sep=""))
    if ( any(!unlist(as.list(control@param.fix)) %in% c(0,1)))
      stop("Only 0 and 1 values are accepted in control@param.fix")
  }
  
  
  
} # }}}

# BBM(FLStock, FLIndices)
setMethod('BBM', signature(object='FLStock'),
          function(object, indices, control, inits)
          {
            
            # check the FLStock has only two age classess
            if(dim(object)[1]!=2)
              stop("'object' must have only two age classes: recruits and adults")
            
            # calculate total catches at age (i.e. for recruits and adults)
            catch <- object@catch.n*object@catch.wt
            
            # BBM(catches='FLQuant', indices='FLIndices', control='BBM.control', inits='BBMpar')
            return(BBM(catch, indices=indices, control=control, inits=inits))
            
          }
)

