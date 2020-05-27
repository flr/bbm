#-------------------------------------------------------------------------------  
# genericMethods.
# Changed: 
#------------------------------------------------------------------------------- 

# generic_ethods - S4 generics
# FLCore/R/generic_methods.R

# Copyright: European Union, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.
# 
# Maintainer: Sonia Sanchez, AZTI



# globalVariables(c("qname"))
# 
# -- OVERLOADED methods/functions {{{

setGeneric("AIC", useAsDefault = stats::AIC)
setGeneric("BIC", useAsDefault = stats::BIC)
setGeneric("plot", useAsDefault = plot)
setGeneric("qqmath", useAsDefault = qqmath)

 # }}}



# -- CONSTRUCTORS, documented with each class {{{


# bbmControl
#' @rdname bbmControl
#' @aliases bbmControl bbmControl-methods
setGeneric("bbmControl", function(object, ...) standardGeneric("bbmControl"))


# bbmFit
#' @rdname bbmFit
#' @aliases bbmFit bbmFit-methods
setGeneric("bbmFit", function(object, ...) standardGeneric("bbmFit"))

# }}}


# bbmFitresiduals
#' @rdname bbmFitresiduals
#' @aliases bbmFitresiduals bbmFitresiduals-methods
setGeneric("bbmFitresiduals", function(object, ...) standardGeneric("bbmFitresiduals"))

# }}}



# -- ACCESSORS {{{

#' @rdname bbmControl-class
#' @aliases g g-methods
setGeneric("g",function(object, ...) standardGeneric("g"))

#' @rdname bbmControl-class
#' @aliases param.fix param.fix-methods
setGeneric("param.fix", function(object, ...) standardGeneric("param.fix"))

#' @rdname bbmFit-class
#' @aliases inputs inputs-methods
setGeneric("inputs", function(object, ...) standardGeneric("inputs"))

#' @rdname bbmFit-class
#' @aliases convergence convergence-methods
setGeneric("convergence", function(object, ...) standardGeneric("convergence"))

#' @rdname bbmFit-class
#' @aliases message message-methods
setGeneric("message", function(object, ...) standardGeneric("message"))

#' @rdname bbmFit-class
#' @aliases fitSumm fitSumm-methods
setGeneric("fitSumm", function(object, ...) standardGeneric("fitSumm"))

#' @rdname bbmFit-class
#' @aliases params.se params.se-methods
setGeneric("params.se", function(object, ...) standardGeneric("params.se"))

#' @rdname bbmFit-class
#' @aliases stock.bio stock.bio-methods
setGeneric("stock.bio", function(object, ...) standardGeneric("stock.bio"))

#' @rdname bbmFit-class
#' @aliases indicesB indicesB-methods
setGeneric("indicesB", function(object, ...) standardGeneric("indicesB"))

#' @rdname bbmFit-class
#' @aliases indicesP indicesP-methods
setGeneric("indicesP", function(object, ...) standardGeneric("indicesP"))

#' @rdname bbmFit-class
#' @aliases residuals.B residuals.B-methods
setGeneric("residuals.B", function(object, ...) standardGeneric("residuals.B"))

#' @rdname bbmFit-class
#' @aliases residuals.P residuals.P-methods
setGeneric("residuals.P", function(object, ...) standardGeneric("residuals.P"))

# }}}



# -- METHODS


# bbm {{{
#' @rdname bbm
#' @aliases bbm bbm-methods
setGeneric("bbm", function(object, indicesB, indicesP, ...)
  standardGeneric("bbm")) 

# }}}


# createInits {{{
#' @rdname createInits
#' @aliases createInits createInits-methods
setGeneric("createInits", function(object, indicesB, indicesP, ...)
  standardGeneric("createInits")) 

# }}}


# simIndices {{{
#' @rdname simIndices
#' @aliases simIndices simIndices-methods
setGeneric("simIndices", function(object, ...)
  standardGeneric("simIndices"))

# }}}



