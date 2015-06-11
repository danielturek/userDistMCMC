
library(nimble)

## load dipper data
load('dipperData.RData')
ind <- c(1,135);   nind<-length(ind);   first<-first[ind];   y<-y[ind,,drop=FALSE];   x_init<-x_init[ind,,drop=FALSE]
yUDD <- 2 - y    ## changes y=1 (alive) to state=1, and y=0 (dead) to state=2

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
##             mean    median         sd  CI95_low  CI95_upp
## nimble 0.5606874 0.5608371 0.02507685 0.5113296 0.6096188
## , , p
##             mean    median         sd  CI95_low  CI95_upp
## nimble 0.8974232 0.8995736 0.02812045 0.8359778 0.9459544
## > out$timing
##         nimble nimble_compile 
##       2.575767       0.572700 





dCR <- nimbleFunction(
    run = function(x = double(1), T = double(2), Z = double(2), length = double(), log.p = double()) {
        pi <- asCol(T[,1])   ## cheap way to get discrete uniform prior
        pi <- pi * 0
        pi <- pi + 1
        pi <- pi/sum(pi)
        logL <- 0
        for(t in 1:length) {
            Zslice <- asCol(Z[x[t], ])
            if(t > 1) 
                logL <- logL + log(sum(pi * Zslice))
            piStar <- pi * Zslice / sum(pi * Zslice)   ## update conditional distribution
            pi <- T %*% piStar   ## propagate state distribution to time (t+1)
        }
        returnType(double())
        if(log.p == 0) return(exp(logL))
        return(logL)
    }
)

rCR <- nimbleFunction(
    run = function(n = integer(), T = double(2), Z = double(2), length = double()) {
        if(n != 1) print('should only specify n=1 in rCR() distribution')
        declare(crHist, double(1, length))
        for(i in 1:length)
            crHist[i] <- 1
        returnType(double(1))
        return(crHist)
    }
)

registerDistributions(list(
    dCR = list(
        BUGSdist = 'dCR(T, Z, length)',
        types    = c('value = double(1)', 'T = double(2)', 'Z = double(2)'),
        discrete = TRUE
    )
))


## UDD dipper model
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
        y[i, first[i]:k] ~ dCR(T[1:2, 1:2], Z[1:2, 1:2], length = k-first[i]+1)
    }
})
constants <- list(k=k, nind=nind, first=first)
data      <- list(y=yUDD)
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





## CdCR <- compileNimble(dCR)
## checkD <- function(ind) {
##     print(y[ind,])
##     x <- y[ind, first[ind]:k]
##     length <- length(x)
##     phi <- 0.6
##     p <- 0.9
##     T <- array(c(phi, 1-phi, 0, 1), c(2,2))
##     Z <- array(c(p,   1-p,   0, 1), c(2,2))
##     print(dCR(x, T, Z, length, 0))
##     print(CdCR(x, T, Z, length, 0))
##     print(exp(dCR(x, T, Z, length, 1)))
##     print(exp(CdCR(x, T, Z, length, 1)))
## }
## checkD(225)   ## 0.2484
## checkD(202)   ## 0.017496
## checkD(167)   ## 0.1246882

## CrCR <- compileNimble(rCR)
## checkR <- function(length) {
##     phi <- 0.6
##     p <- 0.9
##     T <- array(c(phi, 1-phi, 0, 1), c(2,2))
##     Z <- array(c(p,   1-p,   0, 1), c(2,2))
##     print(rCR(1, T, Z, length))
##     print(CrCR(1, T, Z, length))
## }
## checkR(1)
## checkR(2)
## checkR(3)
## checkR(10)







    nind  <- dim(y)[1]   ## number of individual MR histories
    k     <- dim(y)[2]   ## number of sighting occasions
    first <- apply(y, 1, function(MRhistory) min(which(MRhistory == 1)))
    nStates <- 2
    pi0 <- matrix(c(1, 0))   ## prior distribution of state variable at first encounter


run = function() {
    declare(first, double(1, nind))    ## delcare() is necessary here, to ensure 'first' has dimension = 1
    
    
    declare(Tmat, double(2, c(nStates, nStates)))   ## state transition matrix
    Tmat[1, 1] <- model[[phiNode]]
    Tmat[2, 1] <- 1 - model[[phiNode]]
    Tmat[1, 2] <- 0
    Tmat[2, 2] <- 1
    
    declare(Zmat, double(2, c(nStates, nStates)))   ## observation process matrix
    Zmat[1,1] <- model[[pNode]]
    Zmat[2,1] <- 1 - model[[pNode]]
    Zmat[1,2] <- 0
    Zmat[2,2] <- 1
    
    declare(Lvec, double(1, nind))   ## vector of conditional likelihood values
    
    for(ind in 1:nind) {
        pi <- pi0
        for(t in first[ind]:k) {
            declare(Zslice, double(2, c(2,1)))
            for(iState in 1:nStates)     { Zslice[iState, 1] <- Zmat[y[ind, t], iState] }    ## 'slice' of Z matrix relevant to y[ind, t]
            ## initialize or update cumulative likelihood
            if(t == first[ind])     { Lvec[ind] <- 1    ## first encounter (conditioned on y=1), initialize L <- 1
                                  } else                  { Lvec[ind] <- Lvec[ind] * sum(pi * Zslice) }
            piStar <- pi * Zslice / sum(pi * Zslice)   ## update distribution of state vector, conditional on observation y[ind, t]
            pi <- Tmat %*% piStar   ## propagate state distribution to time (t+1)
        }
        ##print('Lvec[ind] = ', Lvec[ind])
    }
    L <- prod(Lvec)
    returnType(double())
    return(log(L))
}




Rnf <- nimbleFunction(
    run = function(a = double(2)) {
        b <- asCol(a[2,])
        returnType(double(2))
        return(b)
    }
)

Cnf <- compileNimble(Rnf)

a <- array(1:6, c(3,2))
a
a[2,]
Rnf(a)
Cnf(a)



