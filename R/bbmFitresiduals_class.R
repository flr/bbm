#-------------------------------------------------------------------------------  
# bbmFitresiduals class.
# Created: Sonia Sanchez - 2018-06-27 13:26:29
# Changed: 
#------------------------------------------------------------------------------- 

# bbmFitresiduals_class.R - output of bbm function
# bbm/R/bbmFitresiduals_class.R

# Copyright: European Union, 2018
# Author: Leire Ibaibarriaga & Sonia Sanchez (AZTI) (<libaibarriaga@azti.es>, <ssanchez@azti.es>)
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# bbmFitresiduals

#' @title bbmFitresiduals class
#'
#' @description \code{bbmFitresiduals} class is used for storing the residuals for the fitted indices in the \code{bbm} function. 
#'              This includes residuals for indices in total biomass (indicesB) and as proportions of recruits in biomass (indicesP).
#' 
#' @name bbmFitresiduals
#' @rdname bbmFitresiduals
#' @docType class
#' @aliases bbmFitresiduals bbmFitresiduals-class residualsB,bbmFitresiduals-method residualsP,bbmFitresiduals-method
#'          plot,bbmFitresiduals-method qqmath,bbmFitresiduals-method
#'
#' @section Slots:
#'  \describe{
#'     \item{residuals.B}{Residuals for indices in total biomass (indicesB), \code{FLQuants}.}
#'     \item{residuals.P}{Residuals for indices as proportions of recruits in biomass (indicesP), \code{FLQuants}.}
#' }
#'
#'
#' @section Accessors:
#' All slots in the class have accessor methods defined that allow retrieving individual slots.
#'
#' @section Methods:
#' Methods exist for various calculations based on values stored in the class. these are:
#'
#' \describe{
#'     \item{plot  }{Method to plot the Pearson residuals.}
#'     \item{qqmath}{Method to do a qqplot of the Pearson residuals of survey indices.}
#' }
#' 
#' 
#' @author Leire Ibaibarriaga & Sonia Sanchez
#' @seealso \link{bbmFit}, \link{plot}, \link{qqmath}
#' @keywords classes
#' 
#' @examples
#' 
#' # Load data
#' data(ane)
#' 
#' # New element of class 'bbmFitresiduals' (empty)
#' bbmFitresiduals()
#' 
#' obj <- bbm( catch.ane, indicesB=indicesB.ane, indicesP=indicesP.ane, control=control.ane, inits=inits.ane)
#' 
#' res <- residuals(obj)
#' class(res)
#' slotNames(res)
#' 
#' # Starndardized residuals plots
#' plot(res)
#' qqmath(res)
#'

setClass("bbmFitresiduals",
         representation(
           residuals.B = "FLQuants",
           residuals.P = "FLQuants"
          ),
         prototype=prototype(
           residuals.B = new('FLQuants'),
           residuals.P = new('FLQuants')
          ),
         validity=function(object){
           # Everything is fine
           return(TRUE)
          }
)


#---------------------
#---------------------

# -- CREATOR

# bbmFitresiduals()   {{{
setMethod('bbmFitresiduals', signature(object='missing'),
          function(object, ...)
          {
            
            args <- list(...)
            
            res <- new("bbmFitresiduals")
            
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

#' @rdname bbmFitresiduals-class
setMethod("residuals.B", "bbmFitresiduals", function(object)  object@residuals.B)

#' @rdname bbmFitresiduals-class
setMethod("residuals.P", "bbmFitresiduals", function(object)  object@residuals.P)

# }}}

#---------------------
#---------------------

# -- METHODS

#' @title Plot of Pearson residuals of survey indices
#' @name plot
#' @docType methods
#' @rdname plot-methods
##' @aliases plot plot,bbmFitresiduals,missing-method
#' @description   Method to produce scatterplots of Pearson residuals of survey indices, if applied to a \code{bbmFitresiduals} object.
#' @param x       An \code{bbmFitresiduals} object with the Pearson residuals.
#' @param y       Ignored.
#' @param auxline A string defining the type of line to be added, by default uses 'smooth', 
#'                a common alternative is to use 'r', a regression, or leave it empty ''.
#' @param ...     Additional argument list.
#' @return If \code{class(x)=='bbmFitresiduals'}, a \code{plot} with stardardized residuals.
#' 
#' @examples
#' 
#' # Residuals
#' res <- residuals(obj)
#' plot(res)
#' 

setMethod("plot", c("bbmFitresiduals", "missing"), function(x, y=missing, auxline="smooth",...){
  
  # Residuals data
  
  residuals.B <- x@residuals.B
  residuals.P <- x@residuals.P
  
  df0 <- rbind(cbind(as.data.frame(residuals.B), idx = "total.biomass (log residuals)"), 
               cbind(as.data.frame(residuals.P), idx = "prop.recruits (residuals)"))
  df0$idx <- as.factor(df0$idx)
  
  # Plot settings
  
  args <- list()
  args$data   <- as.data.frame(df0)
  args$x      <- as.formula("data~factor(year)|factor(idx)*qname")
  args$type   <- c("p", auxline)
  args$groups <- quote(qname)
  args$cex    <- 0.3
  args$lwd    <- 2
  args$ylab   <- "Pearson residuals"
  args$xlab   <- ""
  args$panel  <- function(x,y,...){
    panel.abline(h=0, col.line="gray80")
    panel.xyplot(x,y,...)
  }
  args$par.settings <- list( superpose.symbol=list(col=1, pch=19, cex=0.2), 
                             superpose.line=list(col="gray75", lty=1, lwd=2),
                             strip.background=list(col="gray90"),
                             strip.border=list(col="black"),
                             box.rectangle=list(col="gray90"))
  args$main <- "residuals of survey indices"
  
  if(is(latticeExtra::useOuterStrips, "function")) latticeExtra::useOuterStrips(do.call("xyplot", args)) else do.call("xyplot", args)
  
})

#' @title Normal Q-Q Plot of surveys' indices
#' @name qqplot
#' @docType methods
#' @rdname qqmath-methods
#' @aliases qqmath qqmath,bbmFitresiduals,missing-method
#' @description Method to produce qqplots of Pearson residuals against theoretical distribution of surveys' indices.
#' @param x    An \code{bbmFitresiduals} object with the Pearson residuals.
#' @param data Ignored.
#' @param ...  Additional argument list that might never be used.
#' @return A qqplot with stardardized log residuals against a theoretical distribution.
#' 
#' @examples
#' 
#' # Load data
#' data(ane)
#' 
#' # Run the assessment
#' obj <- bbm( catch.ane, indicesB=indicesB.ane, indicesP=indicesP.ane, control=control.ane, inits=inits.ane)
#' 
#' # Residuals
#' res <- residuals(obj)
#' qqmath(res)
#' 

setMethod("qqmath", c("bbmFitresiduals", "missing"), function(x, data=missing, ...){
  
  # Residuals data
  
  residuals.B <- x@residuals.B
  residuals.P <- x@residuals.P
  
  df0 <- rbind(cbind(as.data.frame(residuals.B), idx = "total.biomass (log residuals)"), 
               cbind(as.data.frame(residuals.P), idx = "prop.recruits (residuals)"))
  df0$idx <- as.factor(df0$idx)
  
  # Plot settings
  
  args      <- list()
  args$data <- as.data.frame(df0)
  args$x    <- as.formula("data~factor(year)|factor(idx)*qname")
  args$ylab <- "Theoretical Quantiles"
  args$xlab <- "Sample Quantiles"
  args$prepanel <- prepanel.qqmathline
  args$panel <- function(x, ...){
    panel.qqmathline(x, col="gray50")
    panel.qqmath(x, ...)
  }
  args$par.settings <- list( strip.background=list(col="gray90")
                             # superpose.symbol=list(col="gray50", pch=19, cex=0.2),
                             # superpose.line=list(col=1, lty=1, lwd=2)
                            )
  args$pch  <- 19
  args$col  <- 1
  args$cex  <- 0.2
  args$main <- "Normal Q-Q Plot of survey indices' residuals"
  
  if(is(latticeExtra::useOuterStrips, "function")) latticeExtra::useOuterStrips(do.call("qqmath", args)) else do.call("qqmath", args)
  
})

