# generics.R - DESC
# /generics.R

# Copyright European Union, 2018
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

# bbm
setGeneric('bbm', function(object, ...) standardGeneric('bbm'))

# bbmfit
setGeneric('bbmFit', function(object, ...) standardGeneric('bbmFit'))

# params.se
setGeneric('params.se', function(object, ...) standardGeneric('params.se'))
setGeneric('params.se<-', function(object, ..., value) standardGeneric('params.se<-'))
