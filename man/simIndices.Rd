% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic_methods.R, R/simIndices_method.R
\name{simIndices}
\alias{simIndices}
\alias{simIndices-methods}
\alias{simIndices,FLQuant-method}
\alias{simIndices,FLStock-method}
\title{Method simIndices}
\usage{
simIndices(object, ...)

\S4method{simIndices}{FLQuant}(object, g, inits, findicesB = NULL, findicesP = NULL)

\S4method{simIndices}{FLStock}(object, g, inits, findicesB = NULL, findicesP = NULL)
}
\arguments{
\item{object}{An \code{FLQuant} with catch information (for recruits and adults).}

\item{g}{A \code{vector} with information on instantaneous rate of biomass decrease, g = M - G, for recruits (rec) and adults (adult).}

\item{inits}{An \code{FLPar} with initial values.}

\item{findicesB}{A \code{vector} with period of the year of the different indices in total biomass.
Optional parameter, if not set, then the survey is assumed to occur at the beggining of the year.}

\item{findicesP}{A \code{vector} with period of the year of the different indices of proportion of recruits in mass.
Optional parameter, if not set, then the survey is assumed to occur at the beggining of the year.}
}
\value{
A list with indices in biomass (Btot) and indices in proportion of recruits (Prec),
both elements of the list are FLIndices.
}
\description{
Function to generate indices of total biomass and of proportion of recruits in biomass,
given information on catches at age, instantaneous rate of biomass decrease (g = M - G ) and values for the bbm parameters.
}
\section{Methods}{

Methods exist for various calculations based on the output class (\code{FLPar}). For details: \code{?FLPar}.
}

\examples{

# Load data
data(ane)

# Case:  object='FLQuant'
indices1 <- simIndices( catch.ane, g=control.ane@g, inits=inits.ane, 
                        findicesB=unlist(lapply( indicesB.ane, function(x) mean(range(x)[c('startf','endf')]))),
                        findicesP=unlist(lapply( indicesP.ane, function(x) mean(range(x)[c('startf','endf')]))) )
class(indices1)
slotNames(indices1)


# Case:  object='FLStock'
stock <- FLStock(catch.n=catch.ane, catch.wt=catch.ane*0+1)
units(stock@catch.wt) <- ''
stock@catch <- quantSums(stock@catch.n*stock@catch.wt)

indices2 <- simIndices( stock, g=control.ane@g, inits=inits.ane, 
                        findicesB=unlist(lapply( indicesB.ane, function(x) mean(range(x)[c('startf','endf')]))),
                        findicesP=unlist(lapply( indicesP.ane, function(x) mean(range(x)[c('startf','endf')]))) )
class(indices2)

# Run assessment with the alternative indices
run  <- bbm(catch.ane, indicesB=indicesB.ane, indicesP=indicesP.ane, control=control.ane, inits=inits.ane)
run1 <- bbm(catch.ane, indicesB=indices1$Btot, indicesP=indices1$Prec, control=control.ane, inits=inits.ane)
run2 <- bbm(catch.ane, indicesB=indices2$Btot, indicesP=indices2$Prec, control=control.ane, inits=inits.ane)

# Plot assessed populations
plot( FLQuants( bc=quantSums(run@stock.bio)[,,,1,], alt1=quantSums(run1@stock.bio)[,,,1,], alt2=quantSums(run2@stock.bio)[,,,1,]))

}
\seealso{
\link{bbm}, \linkS4class{FLQuant}, \linkS4class{FLQuants}, \link{FLIndices}, \linkS4class{bbmControl}, \linkS4class{FLPar}, \link{bbmFLPar}
}
\author{
Leire Ibaibarriaga & Sonia Sanchez.
}
\keyword{FLPar}
\keyword{methods}
