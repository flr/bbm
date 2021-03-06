% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcPop_fun.R
\name{calcPop}
\alias{calcPop}
\title{Function to estimate abundances}
\usage{
calcPop(g, f, catch, inits)
}
\arguments{
\item{g}{A \code{vector} with information on instantaneous rate of biomass decrease, g = M - G, for recruits (rec) and adults (adult).}

\item{f}{A \code{vector} with fraction of the year corresponding to each of the periods (determined by the period of the year of the different indices).}

\item{catch}{An \code{FLQuant} with catch information (for recruits and adults).}

\item{inits}{An \code{FLPar} with parameter values for the parameters required by the \code{bbm} function.}
}
\value{
A \code{list} with two elements: \code{stock.bio}, with information on estimated stock (\code{FLQuant});
and \code{ok}, an indicator on whether estimated parameters are valid (i.e. positive, \code{ok==TRUE}) or not (\code{ok==FALSE}).
}
\description{
This function estimates abundances in mass at age (for recruits and adults) by period, given information on growth, periods duration, catches
and some additional values.
}
\examples{

# Load required libraries
library(bbm)

# Load data
data(ane)

# Generate population estimates, given some estimated parameters
findicesB.ane <- unlist(lapply( indicesB.ane, function(x) mean(range(x)[c('startf','endf')])))
findicesP.ane <- unlist(lapply( indicesP.ane, function(x) mean(range(x)[c('startf','endf')])))
bioAge <- calcPop(g=control.ane@g, 
                  f=periods( findicesB=findicesB.ane, findicesP=findicesP.ane)$f, 
                  catch=catch.ane, inits=inits.ane)
class(bioAge)

# Check if valid ouput (i.e. positive biomass values)
bioAge$ok

# Estimates
bioAge$stock
 

}
\seealso{
\link{FLPar}, \link{bbmFLPar}
}
\author{
Leire Ibaibarriaga & Sonia Sanchez.
}
\keyword{calcPop}
