#-------------------------------------------------------------------------------  
# intsBBM function.
# Created: Sonia Sanchez - 2018-05-24 08:10:31
# Changed: 2018-06-01 10:24:36 (ssanchez)
#------------------------------------------------------------------------------- 

# BBMpar.r - DESC
# bbm/R/BBMpar.r

# Copyright: European Union & AZTI, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


# BBMpar {{{

#' BBMpar class with information on parameters required by the BBM function.
#'
#' @name BBMpar
#' @rdname BBMpar
#' @aliases BBMpar BBMpar-methods BBMpar-class
#'
#' @param object An \code{FLQuant} with catch information (for recruits and adults) or an \code{FLStock}.
#' @param indices Abundance indices in biomass for same age classes (class: FLQuants).
#' @param control List with other extra variables:
#'  \itemize{
#'      \item{logq   }{: Logarithm of the catchability parameter for each survey.}
#'      \item{logpsi }{: Logarithm of the precisions of the observation equations of the total biomass index for each survey.}
#'      \item{xi     }{: Parameter related to the variance of the observation equations of recruits biomass proportion for each survey.}
#'      \item{logB0  }{: Logarithm of the total biomass at the begining of the first year.}
#'      \item{logR   }{: Logarithm of annual recruitments at the beggining of the year.}
#'      \item{mur    }{: Mean of the normal process for log(recruitment).}
#'      \item{logpsir}{: Logarithm of the precision of the normal process for log(recruitment.}
#'  }
#' 
#' @return An object of class BBMpar.
#' @section Validity: 
#'  \describe{
#'     \item{XXX}{XXX}
#' }
#' You can inspect the class validity function by using
#'    \code{getValidity(getClassDef('BBMpar'))}
#'
#'
#' @section Accessors:
#' All slots in the class have accessor and replacement methods defined that
#' allow retrieving and substituting individual slots.
#'
#' The values passed for replacement need to be of the class of that slot.
#' A numeric vector can also be used when replacing FLQuant slots, and the
#' vector will be used to substitute the values in the slot, but not its other
#' attributes.
#' @section Constructor:
#' A construction method exists for this class that can take named arguments for
#' any of its slots. All slots are then created to match the requirements of the
#' class validity. If \code{nyear} or \code{nidex} are provided, this is used
#' for sizing the logR slot and the slots related to the indices (logq, logpsi and xi), respectively.
#' @section Methods:
#' Methods exist for various calculations based on values stored in the class:
#'
#' \describe{
#'     \item{as.list}{: Converts the \code{BBMpar} object into a list.}
#' }
#' @author Leire Ibaibarriaga & Sonia Sanchez.
#' @seealso \link{BBM}
#' @keywords BBMpar
#' @examples 
#' # Load required libraries
#' library(bbm)
#' 
#' # Load data
#' data(ane)
#' 
#' # Generate an object of BBMpar class (different alternatives)
#' new("BBMpar")                # empty object
#' slotNames(BBMpar())          # slots
#' 
#' BBMpar( nyear=20, nindex=3)  # setting dimensions
#' BBMpar( logq=0)              # providing values
#' # BUG BBMpar( logq=0, xi=rep(5,3)) # but with some restrictions
#' # BUG BBMpar( logB0=rep(10,3))
#' # BUG BBMpar( logS=10)
#' 
#' # Available methods for the class
#' # - as.list()
#' 
#' l1 <- as.list(new("BBMpar"))
#' class(l1)
#' 
#' l2 <- base:::as.list(new("BBMpar"))
#' class(l2)
#' 
#' # Generate intitial values (output class is BBMpar)
#' inits <- initsBBM(catch.ane, indices=indices.ane, g=control.ane@g)
#' class(inits)
#' 
#' # Run assessment (inits must be of class BBMpar)
#' run <- BBM(catch.ane, indices=indices.ane, control=control.ane, inits=inits.ane)
#' run
#' 


setClass("BBMpar",
         representation(
           "list",
           logq    = "numeric", # [nindex]
           logpsi  = "numeric", # [nindex]
           xi      = "numeric", # [nindex]
           logB0   = "numeric",
           logR    = "numeric", # [nyear]
           mur     = "numeric",
           logpsir = "numeric"),
         prototype=prototype(
           logq    = as.numeric(NA), # [nindex]
           logpsi  = as.numeric(NA), # [nindex]
           xi      = as.numeric(NA), # [nindex]
           logB0   = as.numeric(NA),
           logR    = as.numeric(NA), # [nyear]
           mur     = as.numeric(NA),
           logpsir = as.numeric(NA)
         ),
         validity=function(object){
           
           # check dimensions
           
           if( length(object@logB0)!=1 | length(object@mur)!=1 | length(object@logpsir)!=1 )
             stop("Slots 'logB0', 'mur' and 'logpsir' must have length 1.")
           
           nyr <- length(object@logq)
           if( length(object@logpsi)!=nyr | length(object@xi)!=nyr )
             stop("All the slots related to indices (these are: logq, logpsi and xi) must have the same length.")
           
           # Everything is fine
           return(TRUE)
         }
)
#! NECESARY TO GENERALISE IT TO MORE THAN 1 ITERATION!


# BBMpar
setGeneric('BBMpar', function(object, ...) standardGeneric('BBMpar'))

# BBMpar()   {{{
setMethod('BBMpar', signature(object='missing'),
          function(object, nyear='missing', nindex='missing',...)
          {
            args <- list(...)
            
            res <- new("BBMpar")
            
            if (!missing(nindex)) {
              res@logq <- res@logpsi <- res@xi <- numeric(nindex)*NA
            }
            
            if (!missing(nyear)) {
              res@logR <- numeric(nyear)*NA
            }
            
            # Load given slots
            if(length(args) > 0)
              for(i in names(args))
                slot(res, i) <- args[[i]]
            
            # check object validity
            validObject(res)
            
            return(res)
          }
) # }}}


# as.list(BBM.par)   {{{
setMethod("as.list", "BBMpar",
          function(x){
            slotnames <- slotNames(x)[!slotNames(x) %in% ".Data"]
            slotlist <- list()
            for(sl in slotnames) {
              slotlist[[sl]] <- slot(x,  sl)
            }
            return(slotlist)
          }
) # }}}
          
