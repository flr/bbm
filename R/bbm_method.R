#-------------------------------------------------------------------------------  
# bbm function.
# Created: Sonia Sanchez - 2018-05-02 11:01:14
# Changed: 2018-06-01 13:24:31 (ssanchez)
#          2018-06-11 09:08:47 (ssanchez) - separated indices into indicesB and indicesP
#------------------------------------------------------------------------------- 

# bbm_method.R - function to simulate a two-stage biomass-based model
# bbm/R/bbm_method.R

# Copyright: European Union, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


# bbm {{{

#' @title Method bbm
#' 
#' @description Function to fit a two-stage biomass-based model.
#' 
#'
#' @name bbm
#' @rdname bbm
#' @aliases bbm bbm-method
#'
#' @param object   An \code{FLQuant} with catch information (in mass) or an \code{FLStock}. These objects must only have two age classes (recruits and adults) 
#'                 and the number of seasons should be 1 or the number of seasons determined by the timing of the different indices.
#' @param indicesB Abundance indices in total biomass (element of class: \code{FLQuants}, \code{FLQuant} or \code{FLIndices}) from surveys. Please assign a survey name to each index.
#' @param indicesP Percentage of recruits in biomass (element of class: \code{FLQuants}, \code{FLQuant} or \code{FLIndices}) from surveys. Please assign a survey name to each proportion.
#' @param findicesB A \code{vector} with fraction of the year corresponding to each of the indicesB.
#' @param findicesP A \code{vector} with fraction of the year corresponding to each of the indicesP.
#' @param control  A \code{bbmControl} with control arguments.
#' @param inits    An \code{FLPar} with initial values.
#' 
#' 
#' @return An object of class \code{bbmFit}.
#' 
#' @author Leire Ibaibarriaga & Sonia Sanchez.
#' @seealso \link{bbmFit}, \linkS4class{FLQuant}, \linkS4class{FLQuants}, \linkS4class{bbmControl}, \linkS4class{FLPar}, \link{bbmFLPar}
#' @keywords bbm methods
#'  
#' @examples
#' 
#' # Load required libraries
#' library(bbm)
#' 
#' # Load data
#' data(ane)

 
#' @rdname bbm
#' @aliases bbm,FLQuant,FLQuants,FLQuants-method
#' @examples
#'
#' # Case:  object='FLQuant'; indicesB=indicesP='FLQuants'; control='bbmControl'; inits='FLPar'
#' 
#' run <- bbm( catch.ane, 
#'             indicesB=lapply( indicesB.ane, function(x) x@index), 
#'             indicesP=lapply( indicesP.ane, function(x) x@index), 
#'             findicesB=unlist(lapply( indicesB.ane, function(x) mean(range(x)[c('startf','endf')]))),
#'             findicesP=unlist(lapply( indicesP.ane, function(x) mean(range(x)[c('startf','endf')]))),
#'             control=control.ane, inits=inits.ane)
#'             
#' is(run)
#' slotNames(run)
#' 
#' # Run the assessment with fixed parameters
#' ctrl <- control.ane
#' ctrl@param.fix['q_depm'] <- 1
#' 
#' run0 <- bbm(catch.ane, indicesB=indicesB.ane, indicesP=indicesP.ane, control=ctrl, inits=inits.ane)
#' params(run0)['q_depm']; inits.ane['q_depm']
#' params.se(run0)['q_depm']
#' 

# bbm(object='FLQuant', indicesB='FLQuants', indicesP='FLQuants',...) {{{
setMethod('bbm', signature(object='FLQuant', indicesB='FLQuants', indicesP='FLQuants'),
          function(object, indicesB, indicesP, findicesB, findicesP, control, inits)
          {

            # check input objects
            
            bbmChecks( object, indicesB, indicesP, findicesB, findicesP, control, inits)
            
            # compute periods and indices timings
            
            aux <- periods(findicesB=findicesB, findicesP=findicesP)
            nper <- aux$nper
            f    <- aux$f
            perindicesB <- aux$perindicesB
            perindicesP <- aux$perindicesP 
            
            # catch
            
            catch <- object
            rm(object)
            
            # if necessary, extend catch to the adequate number of periods
            
            if (dim(catch)[4]==1 & nper>1)  { 
              catch <- expand( catch, season=1:nper)
              for (s in nper:1)
                catch[,,,s,,] <- catch[,,,1,,] * f[s]
            }


            niter <- dims(catch)$iter
            
            out <- bbmFit(years=dimnames(catch)$year, niter=niter, namesB=names(indicesB), namesP=names(indicesP))
            out@inputs <- list(catch=catch, indicesB=indicesB, indicesP=indicesP, 
                               perindicesB=perindicesB, perindicesP=perindicesP, 
                               control=control, f=f, nper=nper)     # this needs to be included only once       
                        
            for (it in 1:niter){
              
              # prepare input data to TMB
              dat.it <- prepareInput( catch=iter(catch, it), indicesB=iter(indicesB, it), indicesP=iter(indicesP, it), 
                                   perindicesB=perindicesB, perindicesP=perindicesP, control=control, f=f, nper=nper)
              # initial parameters
              
              inits.it <- prepareFLPar(iter(inits,it))
              
              param.fix <- prepareFLPar( control@param.fix, param.fix=T)
              
              # model definition, accounting for fixed parameters 

              if (sum(unlist(param.fix))==0) {   #! check this works ok
                model <- MakeADFun(dat.it, inits.it, DLL="bbm")
              } else {
                map <- list()
                for (sl in names(param.fix)) {
                  if(sum(param.fix[[sl]])>0) 
                    map[[sl]] <- param.fix[[sl]]
                }
                map <- lapply( map, function(x) {x[x==1] <- NA; return(as.factor(x))})
                model <- MakeADFun(dat.it, inits.it, DLL="bbm", map=map)
              }
              
              # model fit: minimise the negative log-likelihood (nll)
              
              opt <- optim(model$par, model$fn, model$gr, hessian=T, control=list(maxit=1e+6))
              
              # other results
              
              out@convergence[it] <- opt$convergence 
              out@message <- ifelse(it>1 | is.null(opt$message), "", opt$message)  
              out@fitSumm["nlogl", it] <- opt$value
              out@fitSumm["nobs", it] <- length(dat.it$Bobs) + length(dat.it$Cobs)
              out@fitSumm["nopar", it] <- sum(unlist(param.fix)==0)  ## or length(res$par.fixed)
              
              # sdreport
              
              res <- sdreport(model, getJointPrecision=T)
              
              #optres$par.fixed  # parameters in the fit scale
              
              if (sum(unlist(param.fix))==0){
                out@vcov[,,it] <- res$cov.fixed  # vcov of the parameters in the fit scale
              }else{
                out@vcov[,,it] <- 0
                diag(out@vcov[,,it]) <- 1
                idx <- (unlist(param.fix)!=1)
                out@vcov[idx,idx,it] <- res$cov.fixed
              }
              
              # summary of derived parameters
              
              optres <- summary(res, select="all")
              idx <- unlist(sapply(c("^q$","^psi$","^xi$","^B0$","^R$","^mur$","^psir$"), 
                                   function(x) grep(x, rownames(optres))), use.names=F) # ^matches the beginning of the string and $ the end
              # out@params[unlist(param.fix)!=1, it] <- optres[idx,"Estimate"] # estimated parameters
              # out@params[unlist(param.fix)==1, it] <- unlist(inits.it)[unlist(param.fix)==1] # fixed parameters
              out@params[, it] <- optres[idx,"Estimate"] # estimated and fixed parameters
              # out@params.se[unlist(param.fix)!=1, it] <- optres[idx,"Std. Error"] # estimated parameters
              # out@params.se[unlist(param.fix)==1, it] <- 0 # fixed parameters
              out@params.se[, it] <- optres[idx,"Std. Error"] # estimated and fixed parameters
              
              # #! check in several examples that these internal values are the same as obtained using the calcPop function
              # #! Once that is checked, we can remove the adreport in the cpp code for Bpred and Ppred
              # 
              # kk1 <- calcPop(control@g, f, catch, out@params)
              # kk1$stock.bio[1,,,1,,] - optres["Brec.col(0)" ==rownames(optres),"Estimate"] # vector with estimates (multiplied by q)  
              # kk1$stock.bio[1,,,2,,] - optres["Brec.col(1)" ==rownames(optres),"Estimate"] # vector with estimates (multiplied by q)  
              # kk1$stock.bio[1,,,3,,] - optres["Brec.col(2)" ==rownames(optres),"Estimate"] # vector with estimates (multiplied by q)  
              # kk1$stock.bio[2,,,1,,] - optres["Badu.col(0)" ==rownames(optres),"Estimate"] # vector with estimates (multiplied by q)  
              # kk1$stock.bio[2,,,2,,] - optres["Badu.col(1)" ==rownames(optres),"Estimate"] # vector with estimates (multiplied by q)  
              # kk1$stock.bio[2,,,3,,] - optres["Badu.col(2)" ==rownames(optres),"Estimate"] # vector with estimates (multiplied by q)  
              # 
              # kk2 <- calcIndices(control@g, findicesB, findicesP, catch, out@params)
              # kk <- NULL
              # for ( i in 1:nrow(dat.it$imatBobs)){
              #   print(dat.it$imatBobs[i,])
              #   kk <- c(kk, kk2$indicesB[[dat.it$imatBobs[i,1]]][,dat.it$imatBobs[i,2],,,,])
              # }
              # print(cbind(optres["Bpred" ==rownames(optres),"Estimate"],kk, optres["Bpred" ==rownames(optres),"Estimate"]-kk))
              # kk <- NULL
              # for ( i in 1:nrow(dat.it$imatPobs)){
              #   print(dat.it$imatPobs[i,])
              #   kk <- c(kk, kk2$indicesP[[dat.it$imatPobs[i,1]]][,dat.it$imatPobs[i,2],,,,])
              # }
              # print(cbind(optres["Ppred" ==rownames(optres),"Estimate"],kk, optres["Ppred" ==rownames(optres),"Estimate"]-kk))
              
            } # end of loop for iterations
            
            out@stock.bio <- calcPop(control@g, f, catch, out@params)$stock.bio
            
            aux <- calcIndices(control@g, findicesB, findicesP, catch, out@params)
            out@indicesB <- aux$indicesB
            out@indicesP <- aux$indicesP
            
            return(out)
            
          }
)
# }}}

#---------------------
#---------------------

bbmChecks <- function( object, indicesB=NULL, indicesP=NULL, findicesB, findicesP, control=NULL, inits=NULL) {
  
  
  # Checking that argument's class is correct
  
  if (!class(object) %in% c('FLQuant'))
    stop( "'object' argument must be of class 'FLQuant' or 'FLStock'.")
  
  if (!is.null(indicesB) & !class(indicesB) %in% c('FLQuants'))
    stop( "'indicesB' argument must be an 'FLIndices' or 'FLQuants' object.")
  
  if (!is.null(indicesP) & !class(indicesP) %in% c('FLQuants'))
    stop( "'indicesP' argument must be an 'FLIndices' or 'FLQuants' object.")
  
  if (!is.null(control) & class(control)!='bbmControl')
    stop( "'control' argument must be of class 'bbmControl'.")
  
  if (!is.null(inits) & class(inits)!='FLPar')
    stop( "'inits' argument must be of class 'FLPar'.")
  
  
  # Checking catch
  
  # no NA's in catch
  if (any(is.na(object))){
    stop("NA's not allowed in catch")
  }
  
  # Age dimension
  
  if (dim(object)[1]!=2)
    stop("'object' argument must have only 2 age classes (one for recruits and the other for adults).")
  
  # Unit and area dimensions
  
  if (sum(dim(object)[c(3,5)]>1)>0)
    warning("'object' should have length 1 for dimensions unit and area. Therefore the 1st value will be used.")
  
  # iter dimensions
  
  if (!is.null(inits)) 
    if (!AllEqual(c(dims(inits)$iter, dims(object)$iter)))
      stop("The number of iterations of inits and catch should be equal")
  
  # Season dimension
  
  aux <- periods(findicesB=findicesB, findicesP=findicesP)
  nper <- aux$nper
  
  if (dim(object)[4]!=1 & dim(object)[4]!=nper)
    stop(paste("'object' argument must have 1 or ", nper, "seasons (being this last value determined by the period of the year of the different indices)."))
  
  
  
  # Checking indicesB and indicesP (interactions and consistency with object)
  
  if ( !is.null(indicesB) & !is.null(indicesP) ) {
    
    # Interactions
    
    # - required names
    if (is.null(names(indicesB)) | is.null(names(indicesP)) )
      stop("Please assign survey names to 'indicesB' and 'indicesP'.")
    if (is.null(names(findicesB)) | is.null(names(findicesP)) )
      stop("Please assign survey names to 'findicesB' and 'findicesP'.")
    
    # - f with correct length
    if (length(indicesP)!=length(findicesP)) {
      stop("'findicesP' must have same length as 'indicesP")
    } else if (sum(names(indicesP)!=names(findicesP))>0) {
      stop("'findicesP' must have same names as 'indicesP")
    }
    
    # - at least one common index, with data in indicesP, and same f values
    ok.nam <- ok.dat <- FALSE
    for (i in names(indicesP))
      if (i %in% names(indicesB)) {
        ok.nam <- TRUE
        if (!all(is.na(indicesP[[i]]))) ok.dat <- TRUE
        if (findicesB[[i]]!=findicesP[[i]]) stop(paste("Different f values for '",i,"' in 'findicesB' and 'findicesP'.",sep=""))
      }
    if (ok.dat==FALSE) stop("Required not all NA values in 'indicesP'.")
    if (ok.nam==FALSE) stop("At least one common survey name is required among 'index' and 'indicesP'.")
    
    
    # Dimensions
    
    # Age dimension
    
    if ( sum(unlist(lapply( indicesB, function(x) dim(x)[1])) !=1)>0 | 
         (sum(unlist(lapply( indicesB, function(x) dim(x)[1])) !=1)==0 & sum(unlist(lapply( indicesB, function(x) dimnames(x)[1])) !="all")>0) )
      stop("Each index in 'indicesB' argument must have only 1 age class ('all'), as should give total biomass estimated by the survey.")
    
    if ( sum(unlist(lapply( indicesP, function(x) dim(x)[1])) !=1)>0 | 
         (sum(unlist(lapply( indicesP, function(x) dim(x)[1])) !=1)==0 & sum(unlist(lapply( indicesP, function(x) dimnames(x)[1])) !="all")>0) )
      stop("Each index in 'indicesP' argument must have only 1 age class ('all'), as should give proportion of recruits biomass estimated by the survey.")
    
    # Unit and area dimensions
    
    if (sum(unlist(lapply( indicesB, function(x) dim(x)[c(3,5)]))>1)>0)
      warning("'indicesB' should have length 1 for dimensions unit and area. Therefore the 1st value will be used.")
    if (sum(unlist(lapply( indicesP, function(x) dim(x)[c(3,5)]))>1)>0)
      warning("'indicesP' should have length 1 for dimensions unit and area. Therefore the 1st value will be used.")
    
    # iter dimensions
    
    if (!AllEqual(c(unlist(lapply( indicesB, function(x) dims(x)$iter)), 
                    unlist(lapply( indicesP, function(x) dims(x)$iter)),
                    dims(object)$iter)))
      stop("The number of iterations of indicesB, indicesP and catch should be equal")
    
    # Season dimension
    
    if (sum(unlist(lapply( indicesB, function(x) dim(x)[4]))!=1)>0 | sum(unlist(lapply( indicesP, function(x) dim(x)[4]))!=1)>0 )
      stop("'indicesB' and 'indicesP' must have only one season.")
    
    # Year dimension
    
    nyear <- dim(object)[2]
    # (indicesB, indicesP and object should have same range of years, otherwise add NA's in the indices objects)
    if ( !AllEqual(unlist(lapply( indicesB, function(x) dims(x)$minyear))) | !AllEqual(unlist(lapply( indicesP, function(x) dims(x)$minyear))) |
         !AllEqual(unlist(lapply( indicesB, function(x) dims(x)$maxyear))) | !AllEqual(unlist(lapply( indicesP, function(x) dims(x)$maxyear))) |
         dims(object)$minyear != min(unlist(lapply( indicesB, function(x) dims(x)$minyear))) | dims(object)$minyear != min(unlist(lapply( indicesP, function(x) dims(x)$minyear))) |
         dims(object)$maxyear < max(unlist(lapply( indicesB, function(x) dims(x)$maxyear))) | dims(object)$maxyear < max(unlist(lapply( indicesP, function(x) dims(x)$maxyear))) )
      stop("Check year ranges in the input objects, as years must be equal in 'object', 'indicesB' and 'indicesP'.")
    
    
    # check units
    
    ui <- lapply( indicesB, function(x) units(x) )
    #! necessary to compare also with indicesP?
    if (any(units(object)!=ui))
      stop("Units from 'object' and 'indicesB' must be equal.")
    
  }
  
  
  # check initis and param.fix: same number of parameters and with the same names
  
  if ( !is.null(inits) & !is.null(control)) {
    
    if (dim(inits)[1] != dim(control@param.fix)[1])
      stop("Different number of parameters in 'inits' and 'control@param.fix'.")
    
    if ( sum(dimnames(inits)$params != dimnames(control@param.fix)$params)>0 )
      stop("Parameter names in 'inits' and 'control@param.fix' are not unconsistent.")
    
    # check initis and param.fix: consistence in the number of indices and their names
    
    inits2    <- prepareFLPar(iter(inits,1))
    
    if (!is.null(control)) {
      param.fix <- prepareFLPar( control@param.fix, param.fix=T)
      nindicesB <- length(indicesB)
      if (length(inits2$logq) != nindicesB | length(param.fix$logq) != nindicesB )
        stop("'q' and 'psi' in 'inits' and in 'control@param.fix' must have one value for each index in indicesB.")
      nindicesP <- length(indicesP)
      if (length(inits2$xi) != nindicesP | length(param.fix$xi) != nindicesP) 
        stop("'xi' in 'inits' and in 'control@param.fix' must have one value for each index in indicesP.")
      if (length(inits2$logR) != nyear | length(param.fix$logR) != nyear)
        stop(paste("'R' in 'inits' and in 'control@param.fix' must be of length ", nyear, ", with one value for each year.", sep=""))
      if ( any(!unlist(as.list(control@param.fix)) %in% c(0,1)))
        stop("Only 0 and 1 values are accepted in control@param.fix")
    } else {
      nindicesB <- length(indicesB)
      if (length(inits2$logq) != nindicesB)
        stop("'q' and 'psi' in 'inits' must have one value for each index in indicesB.")
      nindicesP <- length(indicesP)
      if (length(inits2$xi) != nindicesP) 
        stop("'xi' in 'inits' must have one value for each index in indicesP.")
      if (length(inits2$logR) != nyear)
        stop(paste("'R' in 'inits' must be of length ", nyear, ", with one value for each year.", sep=""))
    }
    
  }
  
}

#---------------------
#---------------------

prepareInput <- function(catch, indicesB, indicesP, perindicesB, perindicesP, f, control, nper) {
  
  yrs   <- dimnames(catch)$year  # years
  nyear <- dim(catch)[2]         # number of years
  nindicesB <- length(indicesB)   # number of indices in Biomass
  nindicesP <- length(indicesP)   # number of indices of proportion of recruits
  
  
  # DATA transformation for input
  
  Crec <- catch[1,1:nyear,1,1:nper,1, 1, drop=T]  # recruits catch
  Cadu <- catch[2,1:nyear,1,1:nper,1, 1, drop=T]  # adults catch
  
  Bobs    <- NULL
  imatBobs <- NULL
  Pobs    <- NULL
  imatPobs <- NULL
  
  for (i in 1:nindicesB){
    aux <- indicesB[[i]][,,,,,,drop=T]
    aux <- aux[!is.na(aux)]
    Bobs <- c(Bobs, aux)
    imatBobs <- rbind( imatBobs, cbind( rep(i, length(aux)), match(names(aux), yrs)))
  }
  for (i in 1:nindicesP){
    aux <- indicesP[[i]][,,,,,,drop=T]
    aux <- aux[!is.na(aux)]
    Pobs <- c(Pobs, aux)
    imatPobs <- rbind( imatPobs, cbind( rep(i, length(aux)), match(names(aux), yrs)))
  }
  
  
  # create data list
  
  dat <- list(grec=control@g['rec'], gadu=control@g['adult'], 
              Crec=Crec, Cadu=Cadu,
              Bobs=Bobs, Pobs=Pobs,
              imatBobs=imatBobs, imatPobs=imatPobs,
              perindicesB=perindicesB, perindicesP=perindicesP,
              f = as.double(f))
  
  return(dat)
  
}

#---------------------
#---------------------

prepareFLPar <- function(input, param.fix=FALSE){
  
  if (!class(input) %in% c('FLPar'))
    stop( "'input' argument must be of class 'FLPar'.")
  
  if (dims(input)$iter >1){
    warning( "Only one iteration allowed in 'input'. Only the first iteration will be used")
    input <- iter(input, 1)    
  }
  
  #! include checks for the names

  parnames <- sapply(dimnames(input)$params, FUN=function(x) unlist(strsplit(x,split="_"))[1])
  names(parnames) <- NULL
  if(! any(parnames %in% c("q","psi","xi","B0","R","mur","psir"))){
    stop("Not appropriate parameters in 'FLPar'.")  
  }

  if (!param.fix) {
    obj <- list(logq=log(as.vector(iter(input, it)[parnames %in% "q"])),
                     logpsi=log(as.vector(iter(input, it)[parnames %in% "psi"])),
                     xi=as.vector(iter(input, it)[parnames %in% "xi"]),
                     logB0=log(as.vector(iter(input, it)[parnames %in% "B0"])),
                     logR=log(as.vector(iter(input, it)[parnames %in% "R"])),
                     mur=as.vector(iter(input, it)[parnames %in% "mur"]),
                     logpsir=log(as.vector(iter(input, it)[parnames %in% "psir"])))
    
  }else{
    obj <- list(logq=as.vector(input[parnames %in% "q"]),
                      logpsi=as.vector(input[parnames %in% "psi"]),
                      xi=as.vector(input[parnames %in% "xi"]),
                      logB0=as.vector(input[parnames %in% "B0"]),
                      logR=as.vector(input[parnames %in% "R"]),
                      mur=as.vector(input[parnames %in% "mur"]),
                      logpsir=as.vector(input[parnames %in% "psir"]))
  }
  return(obj)  
}

#---------------------
#---------------------

#' @rdname bbm
#' @aliases bbm,FLStock,ANY,ANY-method
#' @examples
#'
#' # Case:  object='FLStock'; indicesB=indicesP='FLIndices'; control='bbmControl'; inits='FLPar'
#' 
#' stock <- FLStock(catch.n=catch.ane, catch.wt=catch.ane*0+1)
#' units(stock@catch.wt) <- ''
#' stock@catch <- quantSums(stock@catch.n*stock@catch.wt)
#' 
#' run2 <- bbm(stock, indicesB=indicesB.ane, indicesP=indicesP.ane, control=control.ane, inits=inits.ane)
#' 

# bbm(object='FLStock',...)
setMethod('bbm', signature(object='FLStock', indicesB='ANY', indicesP='ANY'),
          function(object, indicesB, indicesP, findicesB=NULL, findicesP=NULL, control, inits)
          {
            
            # check the FLStock has only two age classess
            if(dim(object)[1]!=2)
              stop("'object' must have only two age classes: recruits and adults")
            
            # calculate total catches at age (i.e. for recruits and adults)
            catch <- object@catch.n*object@catch.wt
            
            # bbm(object='FLQuant', indicesB='FLQuant', indicesP='FLQuant', findicesB='vector', findicesP='vector', control='bbmControl', inits='FLPar')
            return(bbm(catch, indicesB=indicesB, indicesP=indicesP, findicesB=findicesB, findicesP=findicesP, control=control, inits=inits))
            
          }
)

#---------------------
#---------------------

#' @rdname bbm
#' @aliases bbm,FLQuant,FLIndices,FLIndices-method
#' @examples
#'
#' # Case:  object='FLQuant'; indicesB=indicesP='FLIndices'; control='bbmControl'; inits='FLPar'
#' 
#' run3 <- bbm(catch.ane, indicesB=indicesB.ane, indicesP=indicesP.ane, control=control.ane, inits=inits.ane) 
#' 

# bbm(indicesB='FLIndices',indicesP='FLIndices',...)
setMethod('bbm', signature(object='FLQuant', indicesB='FLIndices', indicesP='FLIndices'),
          function(object, indicesB, indicesP, findicesB=NULL, findicesP=NULL, control, inits)
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
            
            # bbm(object='FLQuant', indicesB='FLQuant', indicesP='FLQuant', findicesB='vector', findicesP='vector', control='bbmControl', inits='FLPar')
            return(bbm(object, indicesB=indicesB, indicesP=indicesP, findicesB=findicesB, findicesP=findicesP, control=control, inits=inits))
            
          }
)

#---------------------
#---------------------

#' @rdname bbm
#' @aliases bbm,FLQuant,FLQuant,FLQuant-method
#' @examples
#'
#' # Case:  object='FLQuant'; indicesB=indicesP='FLQuant'; control='bbmControl'; inits='FLPar'
#' 
#' namdel <- c("q_acoustic","psi_acoustic","xi_acoustic") # we will take only one of the indices --> need to delete the parameters related to other indices
#' control <- control.ane
#' control@param.fix <- control@param.fix[dimnames(control@param.fix)$params[!dimnames(control@param.fix)$params %in% namdel],]
#' inits   <- inits.ane[dimnames(inits.ane)$params[!dimnames(inits.ane)$params %in% namdel],]
#' run4 <- bbm( catch.ane, indicesB=indicesB.ane[[1]]@index, indicesP=indicesP.ane[[1]]@index, 
#'              findicesB=c( depm=(indicesB.ane[[1]]@range[['startf']]+indicesB.ane[[1]]@range[['endf']])/2), 
#'              findicesP=c( depm=(indicesP.ane[[1]]@range[['startf']]+indicesP.ane[[1]]@range[['endf']])/2), 
#'              control=control, inits=inits)
#' 
#' # Plot assessed populations
#' biomass <- FLQuants()
#' runs <- c('run','run0','run2','run3','run4')
#' names(runs) <- c('bc','fixed_qdepm','run2','run3','only_depm')
#' for (i in 1:length(runs)) biomass[[i]] <- quantSums(stock.bio(get(runs[i])))
#' names(biomass) <- names(runs)
#' plot( biomass)
#' 

# bbm(indicesB='FQuant',indicesP='FQuant',...)
setMethod('bbm', signature(object='FLQuant', indicesB='FLQuant', indicesP='FLQuant'),
          function(object, indicesB, indicesP, findicesB, findicesP, control, inits)
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
              
            # bbm(object='FLQuant', indicesB='FLQuant', indicesP='FLQuant', findicesB='vector', findicesP='vector', control='bbmControl', inits='FLPar')
            return(bbm(object, indicesB=indicesB, indicesP=indicesP, findicesB=findicesB, findicesP=findicesP, control=control, inits=inits))

          }
)

