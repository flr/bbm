#-------------------------------------------------------------------------------  
# bbmControl class.
# Created: Sonia Sanchez - 2018-05-24 09:39:40
# Changed:  ()
#------------------------------------------------------------------------------- 

# bbmControl_class.R - output of bbm function
# bbm/R/bbmControl_class.R

# Copyright: European Union, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# bbmControl

#' @title bbmControl class
#' 
#' @description \code{bbmControl} class is used for setting the control parameters for the \code{bbm} function.
#'
#' @name bbmControl
#' @rdname bbmControl
#' @docType class
#' @aliases bbmControl bbmControl-class g,bbmControl-method param.fix,bbmControl-method
#'
#' @section Slots:
#'  \describe{
#'      \item{g        }{ Instantaneous rate of biomass decrease, g = M - G, for recruits (rec) and adults (adult), 
#'                         \code{c(rec=numeric(1), adult=numeric(1))}.}
#'      \item{param.fix}{ An \code{FLPar} with value 1 if the parameter is fixed at its initial value and 0 otherwise.}
#'  }
#'
#' @section Validity: 
#'  \describe{
#'     \item{Dimensions}{\code{g} must be a vector of length 2 and with the names 'rec' and 'adult'}
#'     \item{Class     }{\code{param.fix} must be of class \code{FLPar}}
#' }
#' You can inspect the class validity function by using \code{getValidity(getClassDef('bbmControl'))}
#'
#' @section Accessors:
#' All slots in the class have accessor methods defined that allow retrieving individual slots.
#' 
#' @section Constructor:
#' A construction method exists for this class that can take named arguments for
#' any of its slots. All slots are then created to match the requirements of the
#' class validity.
#'
#' @author Leire Ibaibarriaga & Sonia Sanchez
#' @seealso \link{bbm}
#' @keywords classes
#' 
#' @examples
#' 
#' # Load data
#' data(ane)
#' 
#' # Generate an object of FLPar class (different alternatives)
#' bbmControl()                     # empty object
#' slotNames(bbmControl())          # slots
#' 
#' bbmControl( g=c(rec=0.68, adult=0.68), param.fix=FLPar(nyear=20, nindex=3)) # setting values for the slots
#' 
#' # Run assessment (control must be of class bbmControl)
#' class(control.ane)
#' run <- bbm(catch.ane, indicesB=indicesB.ane, indicesP=indicesP.ane, control=control.ane, inits=inits.ane)
#' run
#' 


setClass("bbmControl",
         representation(
           g         = "vector", 
           param.fix = "FLPar"
         ),
         prototype=prototype(
           g         = c(rec=NA, adult=NA),
           param.fix = FLPar() # ADAPT TO FLPar
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

#---------------------
#---------------------

# -- CREATOR

# bbmControl()   {{{
#' @rdname bbmControl
#' @aliases bbmControl,missing-method
setMethod('bbmControl', signature(object='missing'),
          # function(object, ...)
          function(object, ...)
          {
            args <- list(...)

            res <- new("bbmControl")
            
            # Load given slots
            if(length(args) > 0)
              for(i in names(args))
                slot(res, i) <- args[[i]]
            
            # check object validity
            validObject(res)
            
            return(res)
          }
) # }}}


#---------------------
#---------------------


# -- ACCESSORS {{{

#' @rdname bbmControl-class
setMethod("g", "bbmControl", function(object) object@inputs)

#' @rdname bbmControl-class
setMethod("param.fix", "bbmControl", function(object)  object@convergence)

# }}}



