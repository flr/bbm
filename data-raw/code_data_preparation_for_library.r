########################################################################################

# preparation for the bbm package of the Anchovy data in Ibaibarriaga et al. (2008)

########################################################################################

# load libraries

library(FLCore)
library(bbm)

########################################################################################

# read real data 1987-2006

bio.dat <- read.table("realdata_87-06_BBM.txt", header=F)
names(bio.dat) <- c("year","h1","h2","c1.per1","ctot.per1","c1.per2","ctot.per2","b1.depm","btot.depm","b1.ac","btot.ac")

# compute age 1 biomass proportion for depm and acoustics

bio.dat$p.depm <- bio.dat$b1.depm/bio.dat$btot.depm
bio.dat$p.ac <- bio.dat$b1.ac/bio.dat$btot.ac

# number of years

n <- dim(bio.dat)[1]

# remove the NA in the catch series to avoid problems with some functions 

bio.dat$c1.per2[n] <- 0.001/2
bio.dat$ctot.per2[n] <- 0.001

# surveys' timing

f <- 0.375

# create index variables for non missing data

years.p.depm <- (1:n)[!is.na(bio.dat$p.depm)]
years.btot.depm <- (1:n)[!is.na(bio.dat$btot.depm)]
years.b1.ac <- (1:n)[!is.na(bio.dat$p.ac)]
years.btot.ac <- (1:n)[!is.na(bio.dat$btot.ac)]

# years vector

years <- 1987:(1987+n-1)

# number of periods

nper <- 2

########################################################################################

# transformation to FLR objects

library(FLCore)

#-----------------------------
# catch FLQuant
#-----------------------------

catch <- FLQuant(dim=c(2, n, 1, nper, 1, 1), dimnames=list(age=1:2, year=years))
catch[1,,,1,,] <- bio.dat$c1.per1
catch[2,,,1,,] <- bio.dat$ctot.per1 - bio.dat$c1.per1
catch[1,,,2,,] <- bio.dat$c1.per2
catch[2,,,2,,] <- bio.dat$ctot.per2 - bio.dat$c1.per2

any(catch<=0)
any(is.na(catch))
catch

#-----------------------------
# indices: 
# one FLIndices for total biomass (indicesB) and 
# another FLIndices for percentage of recruits (indicesP)
#-----------------------------

flq <- FLQuant(dim=c(1, n, 1, 1, 1, 1), dimnames=list(age='all', year=years))

indicesB <- FLIndices(depm=FLIndex(index=flq), acoustic=FLIndex(index=flq))

indicesB$depm@index[] <- bio.dat$btot.depm
indicesB$depm@range["startf"] <- f
indicesB$depm@range["endf"] <- f

indicesB$acoustic@index[] <- bio.dat$btot.ac
indicesB$acoustic@range["startf"] <- f
indicesB$acoustic@range["endf"] <- f


indicesP <- FLIndices(depm=FLIndex(index=flq), acoustic=FLIndex(index=flq))

indicesP$depm@index[] <- bio.dat$p.depm
indicesP$depm@range["startf"] <- f
indicesP$depm@range["endf"] <- f

indicesP$acoustic@index[] <- bio.dat$p.ac
indicesP$acoustic@range["startf"] <- f
indicesP$acoustic@range["endf"] <- f

# Indices as FLQuants (instead of FLIndices)

indicesB_flqs <- FLQuants(depm=flq, acoustic=flq)
indicesB_flqs$depm[]     <- bio.dat$btot.depm
indicesB_flqs$acoustic[] <- bio.dat$btot.ac

indicesP_flqs <- FLQuants(depm=flq, acoustic=flq)
indicesP_flqs$depm[]     <- bio.dat$p.depm
indicesP_flqs$acoustic[] <- bio.dat$p.ac

# timing of the indices:vector(n) with the names of each survey 
#               (if FLIndex available --> f=mean(range(index)[c('startf','endf')]) )

findicesB <- c(depm=f, acoustic=f)
findicesP <- c(depm=f, acoustic=f)

#-----------------------------
# control: 
#-----------------------------

# dummy FLPar indicating which parameters are fixed (0 estimated and 1 fixed)

param.fix <- FLPar(c(rep(0,2), rep(0,2), rep(0,2), 0, rep(0, n), 0, 0),
      params=c(paste("q",names(findicesB),sep="_"),
               paste("psi",names(findicesB),sep="_"),
               paste("xi",names(findicesP),sep="_"),
               "B0",
               paste("R", dimnames(catch)$year,sep="_"),
               "mur","psir"), iter=1)

# named vector with annual biomass change rates

g <- c(rec=0.68, adult=0.68)

# bbmControl object

control <- new( "bbmControl", g=g, param.fix=param.fix)

rm( g, param.fix)

#------------------------
#  Initial values: FLPar
#------------------------

inits <- FLPar(c(rep(1,2), rep(50,2), rep(3,2), 50000, rep(80000, n), 10, 2),
               params=c(paste("q",names(findicesB),sep="_"),
                        paste("psi",names(findicesB),sep="_"),
                        paste("xi",names(findicesP),sep="_"),
                        "B0",
                        paste("R", dimnames(catch)$year,sep="_"),
                        "mur","psir"), iter=1)

########################################################################################

# save the RData in an external file in the data folder to be built-in in the package

catch.ane <- catch
indicesB.ane <- indicesB
indicesB_flqs.ane <- indicesB_flqs
findicesB.ane <- findicesB
indicesP.ane <- indicesP
indicesP_flqs.ane <- indicesP_flqs
findicesP.ane <- findicesP
control.ane <- control
inits.ane <- inits

# save( catch.ane, indicesB.ane, indicesB_flqs.ane, indicesP.ane, indicesP_flqs.ane, findicesB.ane, findicesP.ane, control.ane, inits.ane, 
#       file='../data/ane.RData')

save( catch.ane, indicesB.ane, indicesP.ane, control.ane, inits.ane, 
      file='../data/ane.RData')

########################################################################################

