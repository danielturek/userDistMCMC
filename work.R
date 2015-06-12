
library(nimble)
source('defs.R')

## load dipper data
load('dipperData.RData')
##ind <- 1:2;   nind<-length(ind);   first<-first[ind];   y<-y[ind,,drop=FALSE];   x_init<-x_init[ind,,drop=FALSE]
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

out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = 'nimble',
    niter = 100000,
    makePlot = FALSE
)
out$summary
out$timing
## , , phi
##            mean    median         sd  CI95_low  CI95_upp
## nimble 0.560819 0.5608689 0.02529621 0.5113113 0.6103282
## , , p
##             mean    median         sd  CI95_low  CI95_upp
## nimble 0.8973973 0.8996027 0.02835741 0.8356764 0.9463729
## > out$timing
##         nimble nimble_compile 
##      2.5795667      0.5496667 





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
        y[i, first[i]:k] ~ dDHMM(length=k-first[i]+1, pi=pi[1:2], Z=Z[1:2,1:2], T=T[1:2,1:2], condition=condition[1:2])
    }
})
constants <- list(k=k, nind=nind, first=first, pi=c(1,0), condition=c(1,0))
data      <- list(y=yDHMM)
inits     <- list(phi=0.6, p=0.9)

## testing, delete after here
Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(10000)
samples <- as.matrix(Cmcmc$mvSamples)
apply(samples, 2, mean)
## delete up to here

out <- MCMCsuite(
    code, constants, data, inits,
    MCMCs = 'nimble',
    niter = 100000,
    makePlot = FALSE
)
out$summary
out$timing












