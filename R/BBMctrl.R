#-------------------------------------------------------------------------------  
# BBM.control class.
# Created: Sonia Sanchez - 2018-05-24 09:39:40
# Changed:  ()
#------------------------------------------------------------------------------- 

# BBMctrl.R - output of BBM function
# bbm/R/BBMctrl.R

# Copyright: European Union & AZTI, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# BBM.control

#' BBM.control class for setting the control parameters for the BBM function.
#'
#' @name BBM.control
#' @rdname BBM.control
#' @docType class
#' @aliases BBM.control BBM.control-methods BBM.control-class
#'
#' @section Slots:
#'  \itemize{
#'      \item{g        }{: \code{c(rec=numeric(1), adult=numeric(1))} - instantaneous rate of biomass decrease, g = M - G, for recruits (rec) and adults (adult).}
#'      \item{param.fix}{: \code{BBMpar}                              - with value 1 if the parameter is fixed and 0 otherwise.}
#'  }
#'
#' @section Validity: 
#'  \describe{
#'     \item{Dimensions}{: slot g must be a vector of length 2 and with the names 'rec' and 'adult'.}
#'     \item{Class     }{: slot param.fix must be of class \code{BBMpar}.}
#' }
#' You can inspect the class validity function by using
#'    \code{getValidity(getClassDef('BBM.control'))}
#'
# @section Accessors:
# All slots in the class have accessor and replacement methods defined that
# allow retrieving and substituting individual slots.
#
# The values passed for replacement need to be of the class of that slot.
# A numeric vector can also be used when replacing FLQuant slots, and the
# vector will be used to substitute the values in the slot, but not its other
# attributes.
#
#' @section Constructor:
#' A construction method exists for this class that can take named arguments for
#' any of its slots. All slots are then created to match the requirements of the
#' class validity. If \code{nyears} or \code{nyears} object is provided, this is used
#' for sizing but not stored in any slot.
#'
#'
#' @section Methods:
#' Methods exist for various calculations based on values stored in the class:
#'
#' \describe{
#'     \item{xxx}{: XXX.}
#' }
#'
#' @author Leire Ibaibarriaga & Sonia Sanchez
#' @seealso \link{BBM}
#' @keywords classes
#' 
#' @examples
#' 
#' # Load required libraries
#' library(bbm)
#' 
#' # Load data
#' data(ane)
#' 
#' # Generate an object of BBMpar class (different alternatives)
#' BBM.control()                # empty object
#' slotNames(BBM.control())          # slots
#' 
#' BBM.control( g=c(rec=0.68, adult=0.68), param.fix=BBMpar(nyear=20, nindex=3)) # setting values for the slots
#' # BBM.control( g=c(adult=0.68))                                                 # but with some restrictions
#' BBM.control( g=c(rec=0.68, adult=0.68))
#' 
#' # Run assessment (control must be of class BBM.control)
#' class(control.ane)
#' run <- BBM(catch.ane, indices=indices.ane, control=control.ane, inits=inits.ane)
#' run
#'
#' 


setClass("BBM.control",
  representation(
    "list",  # - Is it required different parameter by iteration?
           g         = "vector", 
           param.fix = "BBMpar"
         ),
         prototype=prototype(
           g         = c(rec=NA, adult=NA),
           param.fix = BBMpar()
         ),
         validity=function(object){
           
           # check names of slot g
           if( length(object@g) != 2 | is.null(names(object@g)) | sum(!names(object@g) %in% c("rec", "adult"))>0 )
             stop("Check slot g.")
           
           # check param.fix validity
           validObject(object@param.fix)
           
           # Everything is fine
           return(TRUE)
          
          }
)


# BBM.control
setGeneric('BBM.control', function(object, ...) standardGeneric('BBM.control'))

# BBM.control()   {{{
#' @rdname BBM.control
#' @aliases BBM.control,missing-method BBM.control,FLQuant-method BBM.controlcpp-class
setMethod('BBM.control', signature(object='missing'),
          # function(object, ...)
          function(object, ...)
          {
            args <- list(...)

            res <- new("BBM.control")
            
            # Load given slots
            if(length(args) > 0)
              for(i in names(args))
                slot(res, i) <- args[[i]]
            
            # check object validity
            validObject(res)
            
            return(res)
          }
) # }}}



