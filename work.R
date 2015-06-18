
library(nimble)
library(coda)
source('defs.R')
load('dipperData.RData')
##ind <- c(1:20, 250:255);   nind<-length(ind);   first<-first[ind];   y<-y[ind,,drop=FALSE];   x_init<-x_init[ind,,drop=FALSE]
yDHMM <- 2 - y    ## changes y=1 (alive) to state=1, and y=0 (dead) to state=2


## regular dipper model
code <- nimbleCode({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for (i in 1:nind) {
        x[i, first[i]] <- 1
        for (t in (first[i] + 1):k) {
            mu_x[i, t] <- phi * x[i, t-1]
            mu_y[i, t] <- p * x[i, t]
            x[i, t] ~ dbin(mu_x[i, t], 1)
            y[i, t] ~ dbin(mu_y[i, t], 1)
        }
    }
})
constants <- list(k=k, nind=nind, first=first)
data      <- list(y=y)
inits     <- list(phi=0.6, p=0.9, x=x_init)


## DHMM dipper model
code <- nimbleCode({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    T[1, 1] <- phi
    T[2, 1] <- 1 - phi
    T[1, 2] <- 0
    T[2, 2] <- 1
    Z[1, 1] <- p
    Z[2, 1] <- 1 - p
    Z[1, 2] <- 0
    Z[2, 2] <- 1
    for (i in 1:nind) {
        y[i, first[i]:k] ~ dDHMM(length=k-first[i]+1, prior=prior[1:2], Z=Z[1:2,1:2], T=T[1:2,1:2], condition=condition[1:2])
    }
})
constants <- list(k=k, nind=nind, first=first, prior=c(1,0), condition=c(1,0))
data      <- list(y=yDHMM)
inits     <- list(phi=0.6, p=0.9)


## MCMC Suite
out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = 'nimble',
    niter = 100000,
    summaryStats = c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp', 'effectiveSize'),
    makePlot = FALSE
)

## output
out$timing
out$summary[, c('mean','sd','CI95_low','CI95_upp'), ]
out$summary[, 'effectiveSize', ]
out$summary[, 'effectiveSize', ] / out$timing['nimble']

## regular dipper MCMC
## > out$timing
##         nimble nimble_compile 
##         2.5737         0.5617 
## > out$summary[, c('mean','sd','CI95_low','CI95_upp'), ]
##                 phi          p
## mean     0.56124176 0.89709195
## sd       0.02499421 0.02882952
## CI95_low 0.51215809 0.83408893
## CI95_upp 0.61001026 0.94655026
## > out$summary[, 'effectiveSize', ]
##       phi         p 
## 14637.710  6632.726 
## > out$summary[, 'effectiveSize', ] / out$timing['nimble']
##      phi        p 
## 5687.419 2577.117 

## using dDHMM
## > out$timing
##         nimble nimble_compile 
##       4.376767       0.173100 
## > out$summary[, c('mean','sd','CI95_low','CI95_upp'), ]
##                phi          p
## mean     0.5612347 0.89571264
## sd       0.0249737 0.02895505
## CI95_low 0.5125123 0.83314978
## CI95_upp 0.6098984 0.94551182
## > out$summary[, 'effectiveSize', ]
##      phi        p 
## 19435.81 18234.94 
## > out$summary[, 'effectiveSize', ] / out$timing['nimble']
##      phi        p 
## 4440.677 4166.305 











