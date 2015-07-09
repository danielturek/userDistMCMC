
############################
## dipper
############################

load('~/GitHub/legacy/dipper/dipperData.RData')
## optionally truncate data:
last <- apply(y, 1, function(hist) max(which(hist==1)))
yDHMM <- 2 - y
if(trunc) { ind <- 1:3;   nind<-length(ind);   first<-first[ind];   last<-last[ind];   y<-y[ind,,drop=FALSE];   yDHMM<-yDHMM[ind,,drop=FALSE];   x_init<-x_init[ind,,drop=FALSE] }

## dipper (with latent states, suitable for jags)
code <- quote({
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
dipper <- list(code=code, constants=constants, data=data, inits=inits)


## dipeprDHMM (y[i] ~ dDHMM(...) for nimble only)
code <- quote({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    T[1,1,1] <- phi
    T[2,1,1] <- 1 - phi
    T[1,2,1] <- 0
    T[2,2,1] <- 1
    T[1:2,1:2,2] <- nimArray(0, 2, 2)
    Z[1,1,1] <- p
    Z[2,1,1] <- 1 - p
    Z[1,2,1] <- 0
    Z[2,2,1] <- 1
    Z[1:2,1:2,2] <- nimArray(0, 2, 2)
    for (i in 1:nind) {
        y[i, first[i]:k] ~ dDHMM(length = k-first[i]+1,
                                 prior = prior[1:2],
                                 condition = condition[1:2],
                                 Z = Z[1:2,1:2,1:2], useZt = 0,
                                 T = T[1:2,1:2,1:2], useTt = 0,
                                 mult = 1)
    }
})
constants <- list(k=k, nind=nind, first=first, prior=c(1,0), condition=c(1,0))
data      <- list(y=yDHMM)
inits     <- list(phi=0.6, p=0.9)
dipperDHMM <- list(code=code, constants=constants, data=data, inits=inits)


## dipperCJS (y[i] ~ dCJS(...) for nimble only)
code <- quote({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for (i in 1:nind) {
        y[i, first[i]:k] ~ dCJS(length=k-first[i]+1, last=last[i]-first[i]+1, phi=phi, p=p)
    }
})
constants <- list(k=k, nind=nind, first=first, last=last)
data      <- list(y=y)
inits     <- list(phi=0.6, p=0.9)
dipperCJS <- list(code=code, constants=constants, data=data, inits=inits)




############################
## goose
############################




## goose data (multistate CR from "handbook of CR analyses", with large multiplicities)
y <- matrix(c(
    ## 0,0,0,1,
    ## 0,0,0,2,
    ## 0,0,0,3,
    0,0,1,0,
    0,0,1,1,
    0,0,1,2,
    0,0,2,0,
    0,0,2,1,
    0,0,2,2,
    0,0,2,3,
    0,0,3,0,
    0,0,3,1,
    0,0,3,2,
    0,0,3,3,
    0,1,0,0,
    0,1,0,1,
    0,1,0,2,
    0,1,1,0,
    0,1,1,1,
    0,1,1,2,
    0,1,2,0,
    0,1,2,1,
    0,1,2,2,
    0,1,2,3,
    0,1,3,0,
    0,1,3,2,
    0,2,0,0,
    0,2,0,1,
    0,2,0,2,
    0,2,0,3,
    0,2,1,0,
    0,2,1,1,
    0,2,1,2,
    0,2,2,0,
    0,2,2,1,
    0,2,2,2,
    0,2,2,3,
    0,2,3,0,
    0,2,3,1,
    0,2,3,2,
    0,2,3,3,
    0,3,0,0,
    0,3,0,1,
    0,3,0,2,
    0,3,0,3,
    0,3,1,0,
    0,3,1,1,
    0,3,1,3,
    0,3,2,0,
    0,3,2,1,
    0,3,2,2,
    0,3,2,3,
    0,3,3,0,
    0,3,3,1,
    0,3,3,2,
    0,3,3,3,
    1,0,0,0,
    1,0,0,1,
    1,0,0,2,
    1,0,0,3,
    1,0,1,0,
    1,0,1,1,
    1,0,1,2,
    1,0,2,0,
    1,0,2,1,
    1,0,2,2,
    1,0,3,1,
    1,0,3,2,
    1,0,3,3,
    1,1,0,0,
    1,1,0,1,
    1,1,0,2,
    1,1,1,0,
    1,1,1,1,
    1,1,1,2,
    1,1,2,0,
    1,1,2,1,
    1,1,2,2,
    1,2,0,0,
    1,2,0,1,
    1,2,0,2,
    1,2,1,0,
    1,2,1,1,
    1,2,1,2,
    1,2,2,0,
    1,2,2,1,
    1,2,2,2,
    1,2,3,0,
    1,3,0,0,
    1,3,1,0,
    2,0,0,0,
    2,0,0,1,
    2,0,0,2,
    2,0,0,3,
    2,0,1,0,
    2,0,1,1,
    2,0,1,2,
    2,0,2,0,
    2,0,2,1,
    2,0,2,2,
    2,0,2,3,
    2,0,3,0,
    2,0,3,2,
    2,1,0,0,
    2,1,0,1,
    2,1,0,2,
    2,1,1,0,
    2,1,1,1,
    2,1,1,2,
    2,1,2,0,
    2,1,2,1,
    2,1,2,2,
    2,1,3,1,
    2,2,0,0,
    2,2,0,1,
    2,2,0,2,
    2,2,1,0,
    2,2,1,1,
    2,2,1,2,
    2,2,2,0,
    2,2,2,1,
    2,2,2,2,
    2,2,3,0,
    2,2,3,2,
    2,3,0,0,
    2,3,2,0,
    2,3,2,1,
    2,3,3,0,
    3,0,0,0,
    3,0,0,1,
    3,0,0,2,
    3,0,0,3,
    3,0,1,0,
    3,0,1,2,
    3,0,2,0,
    3,0,2,2,
    3,0,2,3,
    3,0,3,0,
    3,0,3,2,
    3,0,3,3,
    3,1,0,0,
    3,2,0,0,
    3,2,0,2,
    3,2,2,0,
    3,2,2,1,
    3,2,2,2,
    3,2,3,0,
    3,2,3,3,
    3,3,0,0,
    3,3,0,1,
    3,3,0,2,
    3,3,0,3,
    3,3,2,0,
    3,3,3,0,
    3,3,3,2,
    3,3,3,3),
            ncol = 4, byrow = TRUE)
mult <- c(
    ## 158,
    ## 352,
    ## 271,
    317,
    62,
    31,
    748,
    56,
    234,
    7,
    504,
    14,
    71,
    72,
    643,
    57,
    36,
    98,
    42,
    10,
    70,
    13,
    25,
    1,
    1,
    2,
    1040,
    27,
    117,
    11,
    59,
    2,
    12,
    285,
    11,
    100,
    2,
    10,
    1,
    5,
    2,
    694,
    7,
    66,
    60,
    13,
    1,
    2,
    59,
    2,
    13,
    2,
    86,
    1,
    12,
    31,
    1135,
    40,
    39,
    2,
    63,
    22,
    8,
    68,
    8,
    13,
    1,
    2,
    1,
    201,
    25,
    11,
    56,
    38,
    1,
    21,
    4,
    1,
    81,
    9,
    12,
    11,
    5,
    4,
    22,
    6,
    4,
    1,
    3,
    1,
    1559,
    18,
    85,
    4,
    33,
    6,
    3,
    163,
    3,
    51,
    1,
    7,
    1,
    55,
    6,
    2,
    9,
    1,
    4,
    7,
    1,
    4,
    1,
    423,
    10,
    45,
    5,
    5,
    4,
    121,
    2,
    48,
    6,
    1,
    16,
    4,
    1,
    1,
    473,
    1,
    17,
    21,
    9,
    1,
    19,
    2,
    1,
    20,
    1,
    7,
    6,
    22,
    3,
    5,
    1,
    1,
    2,
    1,
    57,
    2,
    2,
    15,
    3,
    10,
    3,
    13)
## optional truncate:
if(trunc) { ind <- 1:3;     y<-y[ind,];     mult<-mult[ind] }
y[which(y == 0)] <- 4
first <- apply(y, 1, function(hist) which(hist!=4)[1])
nind <- dim(y)[1]
k <- dim(y)[2]
##x_init <- array(as.numeric(NA), dim(y))
##for(i in 1:dim(x_init)[1]) {
##    x_init[i, first[i]:k] <- y[i, first[i]]
##}
yExp <- array(0, c(sum(mult), k))
cur <- 0
for(i in seq_along(mult)) {
    for(j in 1:mult[i])
        yExp[cur+j,] <- y[i,]
    cur <- cur + mult[i]
}
firstExp <- apply(yExp, 1, function(hist) which(hist!=4)[1])
nindExp <- dim(yExp)[1]
kExp <- dim(yExp)[2]
x_knownExp <- yExp
x_knownExp[x_knownExp==4] <- NA
for(i in 1:dim(x_knownExp)[1]) {
    x_knownExp[i, firstExp[i]] <- NA
}
### below are the correct initial values for x for jags
### this took a long time
### they must result in a valid model logProb
x_initExp <- yExp
for(i in 1:dim(x_initExp)[1]) {
    x_initExp[i, 1:firstExp[i]] <- NA
    lastState <- yExp[i,firstExp[i]]
    for(j in (firstExp[i]+1):k) {
        if(x_initExp[i,j] == 4) {
            x_initExp[i,j] <- lastState
        } else {
            lastState <- x_initExp[i,j]
            x_initExp[i,j] <- NA
        }
    }
}






## gooesDHMM (multistate, y[i] ~ dDHMM(...) for nimble only, plus is has large multiplicities)
code <- quote({
    for(i in 1:6)
        p[i] ~ dunif(0, 1)
    for(t in 1:3) {
        Z[1,1,t] <- p[1]
        Z[2,1,t] <- 0
        Z[3,1,t] <- 0
        Z[4,1,t] <- 1 - p[1]
        Z[1,2,t] <- 0
        Z[2,2,t] <- p[2]
        Z[3,2,t] <- 0
        Z[4,2,t] <- 1 - p[2]
        Z[1,3,t] <- 0
        Z[2,3,t] <- 0
        Z[3,3,t] <- p[3]
        Z[4,3,t] <- 1 - p[3]
        Z[1,4,t] <- 0
        Z[2,4,t] <- 0
        Z[3,4,t] <- 0
        Z[4,4,t] <- 1
    }
    Z[1,1,4] <- p[4]
    Z[2,1,4] <- 0
    Z[3,1,4] <- 0
    Z[4,1,4] <- 1 - p[4]
    Z[1,2,4] <- 0
    Z[2,2,4] <- p[5]
    Z[3,2,4] <- 0
    Z[4,2,4] <- 1 - p[5]
    Z[1,3,4] <- 0
    Z[2,3,4] <- 0
    Z[3,3,4] <- p[6]
    Z[4,3,4] <- 1 - p[6]
    Z[1,4,4] <- 0
    Z[2,4,4] <- 0
    Z[3,4,4] <- 0
    Z[4,4,4] <- 1
    for(i in 1:3)
        phi[i] ~ dunif(0, 1)
    for(i in 1:2)
        for(j in 1:3)
            for(l in 1:2) {
                alpha[i,j,l] ~ dnorm(0, sd = 100)
                exal[i,j,l] <- exp(alpha[i,j,l])
            }
    for(j in 1:3)
        for(l in 1:2) {
            psi[1,j,l] <- exal[1,j,l] / (1 + exal[1,j,l] + exal[2,j,l])
            psi[2,j,l] <- exal[2,j,l] / (1 + exal[1,j,l] + exal[2,j,l])
            psi[3,j,l] <-           1 / (1 + exal[1,j,l] + exal[2,j,l])
        }
    for(t in 1:2)
        for(j in 1:3) {
            for(i in 1:3)
                T[i,j,t] <- phi[j] * psi[i,j,1]
            T[4,j,t] <- 1 - phi[j]
        }
    for(t in 3:4)
        for(j in 1:3) {
            for(i in 1:3)
                T[i,j,t] <- phi[j] * psi[i,j,2]
            T[4,j,t] <- 1 - phi[j]
        }
    for(t in 1:4) {
        T[1,4,t] <- 0
        T[2,4,t] <- 0
        T[3,4,t] <- 0
        T[4,4,t] <- 1
    }
    for (i in 1:nind)
        y[i, first[i]:k] ~ dDHMM(length = k-first[i]+1,
                                 prior = prior[1:4],
                                 condition = condition[1:4],
                                 Z = Z[1:k,1:k,first[i]:k], useZt = 1,
                                 T = T[1:k,1:k,first[i]:k], useTt = 1,
                                 mult = mult[i])
})
constants <- list(nind=nind, k=k, first=first, mult=mult,
                  prior=c(1/3, 1/3, 1/3, 0),
                  condition = c(1, 1, 1, 0))
data <- list(y=y)
inits <- list(p=rep(1/2,6), phi=rep(1/2,3), alpha=array(0,c(2,3,2)))
gooseDHMM <- list(code=code, constants=constants, data=data, inits=inits)






## gooesExp (multistate, x[i] latent states, large multiplicities expanded for jags run)
gooseExpJAGSfunction <- function(niter = 100000) {
    sink('goose.jags')
    cat('
model {
    for(i in 1:6) {
        p[i] ~ dunif(0, 1)
    }
    for(t in 1:3) {
        Z[1,1,t] <- p[1]
        Z[2,1,t] <- 0
        Z[3,1,t] <- 0
        Z[4,1,t] <- 1 - p[1]
        Z[1,2,t] <- 0
        Z[2,2,t] <- p[2]
        Z[3,2,t] <- 0
        Z[4,2,t] <- 1 - p[2]
        Z[1,3,t] <- 0
        Z[2,3,t] <- 0
        Z[3,3,t] <- p[3]
        Z[4,3,t] <- 1 - p[3]
        Z[1,4,t] <- 0
        Z[2,4,t] <- 0
        Z[3,4,t] <- 0
        Z[4,4,t] <- 1
    }
    Z[1,1,4] <- p[4]
    Z[2,1,4] <- 0
    Z[3,1,4] <- 0
    Z[4,1,4] <- 1 - p[4]
    Z[1,2,4] <- 0
    Z[2,2,4] <- p[5]
    Z[3,2,4] <- 0
    Z[4,2,4] <- 1 - p[5]
    Z[1,3,4] <- 0
    Z[2,3,4] <- 0
    Z[3,3,4] <- p[6]
    Z[4,3,4] <- 1 - p[6]
    Z[1,4,4] <- 0
    Z[2,4,4] <- 0
    Z[3,4,4] <- 0
    Z[4,4,4] <- 1
    for(i in 1:3) {
        phi[i] ~ dunif(0, 1)
    }
    for(i in 1:2) {
        for(j in 1:3) {
            for(l in 1:2) {
                alpha[i,j,l] ~ dnorm(0, 0.0001)
                exal[i,j,l] <- exp(alpha[i,j,l])
            }
        }
    }
    for(j in 1:3) {
        for(l in 1:2) {
            psi[1,j,l] <- exal[1,j,l] / (1 + exal[1,j,l] + exal[2,j,l])
            psi[2,j,l] <- exal[2,j,l] / (1 + exal[1,j,l] + exal[2,j,l])
            psi[3,j,l] <-           1 / (1 + exal[1,j,l] + exal[2,j,l])
        }
    }
    for(t in 1:2) {
        for(j in 1:3) {
            for(i in 1:3) {
                T[i,j,t] <- phi[j] * psi[i,j,1]
            }
            T[4,j,t] <- 1 - phi[j]
        }
    }
    for(t in 3:4) {
        for(j in 1:3) {
            for(i in 1:3) {
                T[i,j,t] <- phi[j] * psi[i,j,2]
            }
            T[4,j,t] <- 1 - phi[j]
        }
    }
    for(t in 1:4) {
        T[1,4,t] <- 0
        T[2,4,t] <- 0
        T[3,4,t] <- 0
        T[4,4,t] <- 1
    }
    for(i in 1:nind) {
        x[i,first[i]] <- y[i,first[i]]
        for(t in (first[i]+1):k) {
            x[i,t] ~ dcat(T[1:4, x[i,t-1], t-1])
            y[i,t] ~ dcat(Z[1:4, x[i,t], t])
        }
    }
}
',fill=TRUE) 
    sink() 
    jags_data <- list(y=yExp, first=firstExp, k=kExp, nind=nindExp, x=x_knownExp)
    initFunction <- function() list(p=rep(1/2,6), phi=rep(1/2,3), alpha=array(0,c(2,3,2)), x=x_initExp)
    ##parameters <- c('p','phi','psi')
    parameters <- c('p','phi','alpha')
    library(R2WinBUGS) 
    library(R2jags) 
    t <- system.time({ls <- jags(jags_data, initFunction, parameters, "goose.jags", n.chains = 1, n.thin = 1, n.iter = niter, n.burnin = 2000, working.directory = getwd())})
    ar <- ls$BUGSoutput$sims.array    ## using sims.array is consistent with MCMCsuite
    ar <- ar[,1,]     ## drop middle index
    ar <- ar[, -1]    ## drop deviance
    summaryStats <- c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp', 'effectiveSize')
    library(coda)
    CI95_low <- function(x) quantile(x, probs = 0.025)
    CI95_upp <- function(x) quantile(x, probs = 0.975)
    nSummaryStats <- length(summaryStats)
    nMonitorNodes <- dim(ar)[2]
    summaryArray <- array(NA, c(nSummaryStats, nMonitorNodes))
    dimnames(summaryArray) <- list(summaryStats, dimnames(ar)[[2]])
    summaryStatFunctions <- lapply(summaryStats, function(txt) eval(parse(text=txt)[[1]]))
    for(iStat in seq_along(summaryStats)) {
        summaryArray[iStat, ] <- apply(ar, 2, summaryStatFunctions[[iStat]])
    }
    theTime <- as.numeric(t[3] / 60)    ## this is consistent with MCMCsuite
    list(theTime=theTime, summary=summaryArray)
}






############################
## orchid
############################



## orchid (multistate, from Kery & Schaub, uses stochastic indexing => jags only!)
## (9.7. Real data example: the showy lady's slipper)
orchidJAGSfunction <- function(niter = 100000) {
    CH <- as.matrix(read.table("~/GitHub/userDistMCMC/orchids.txt", sep=" ", header = FALSE))
    ## optionally truncate data:
    if(trunc) {     ind <- 1:5;     CH <- CH[ind,]     }
    n_occasions <- dim(CH)[2]
    ## Compute vector with occasion of first capture 
    f <- numeric() 
    for (i in 1:dim(CH)[1]){f[i] <- min(which(CH[i,]!=0))}
    ## Recode CH matrix: note, a 0 is not allowed by WinBUGS!
    ## 1 = seen vegetative, 2 = seen flowering, 3 = not seen 
    rCH <- CH  # Recoded CH 
    rCH[rCH==0] <- 3 
    ## Specify model in BUGS language 
    sink("ladyslipper.jags") 
    cat("
model {
    ## -------------------------------------------------
    ## Parameters:
    ## s: survival probability 
    ## psiV: transitions from vegetative 
    ## psiF: transitions from flowering 
    ## psiD: transitions from dormant 
    ## -------------------------------------------------
    ## States (S):
    ## 1 vegetative 
    ## 2 flowering 
    ## 3 dormant 
    ## 4 dead 
    ## Observations (O):  
    ## 1 seen vegetative 
    ## 2 seen flowering 
    ## 3 not seen 
    ## -------------------------------------------------
    ## Priors and constraints 
    ## Survival: uniform 
    for (t in 1:(n_occasions-1)){  
        s[t] ~ dunif(0, 1) 
    }
    ## Transitions: gamma priors 
    for (i in 1:3){
        a[i] ~ dgamma(1, 1) 
        psiD[i] <- a[i]/sum(a[1:3]) 
        b[i] ~ dgamma(1, 1) 
        psiV[i] <- b[i]/sum(b[1:3]) 
        c[i] ~ dgamma(1, 1) 
        psiF[i] <- c[i]/sum(c[1:3]) 
    }
    ## Define state-transition and observation matrices 	
    for (i in 1:nind){
        ## Define probabilities of state S(t+1) given S(t) 
        for (t in 1:(n_occasions-1)){
            ps[1,i,t,1] <- s[t] * psiV[1]
            ps[1,i,t,2] <- s[t] * psiV[2]
            ps[1,i,t,3] <- s[t] * psiV[3]
            ps[1,i,t,4] <- 1-s[t]
            ps[2,i,t,1] <- s[t] * psiF[1]
            ps[2,i,t,2] <- s[t] * psiF[2]
            ps[2,i,t,3] <- s[t] * psiF[3]
            ps[2,i,t,4] <- 1-s[t]
            ps[3,i,t,1] <- s[t] * psiD[1]
            ps[3,i,t,2] <- s[t] * psiD[2]
            ps[3,i,t,3] <- s[t] * psiD[3]
            ps[3,i,t,4] <- 1-s[t]
            ps[4,i,t,1] <- 0 
            ps[4,i,t,2] <- 0 
            ps[4,i,t,3] <- 0 
            ps[4,i,t,4] <- 1 
            ## Define probabilities of O(t) given S(t) 
            po[1,i,t,1] <- 1 
            po[1,i,t,2] <- 0 
            po[1,i,t,3] <- 0 
            po[2,i,t,1] <- 0 
            po[2,i,t,2] <- 1 
            po[2,i,t,3] <- 0 
            po[3,i,t,1] <- 0 
            po[3,i,t,2] <- 0 
            po[3,i,t,3] <- 1 
            po[4,i,t,1] <- 0 
            po[4,i,t,2] <- 0 
            po[4,i,t,3] <- 1 
        } #t 
    } #i 
    ## Likelihood 
    for (i in 1:nind){
        ## Define latent state at first capture 
        z[i,f[i]] <- y[i,f[i]]
        for (t in (f[i]+1):n_occasions){
            ## State process: draw S(t) given S(t-1) 
            z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:4]) 
            ## Observation process: draw O(t) given S(t) 
            y[i,t] ~ dcat(po[z[i,t], i, t-1, 1:3]) 
        } #t 
    } #i
}
",fill=TRUE) 
    sink() 
    ## Function to create known latent states z 
    known_state_ms <- function(ms, notseen){
        ## notseen: label for not seen 
        state <- ms 
        state[state==notseen] <- NA 
        for (i in 1:dim(ms)[1]){
            m <- min(which(!is.na(state[i,]))) 
            state[i,m] <- NA 
        }
        return(state) 
    }
    ## Bundle data 
    jags_data <- list(y = rCH, f = f, n_occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known_state_ms(rCH, 3)) 
    ## Note that all initial values of the unknown z must all be 3 (not 1 or 2). Therefore the function ms_init_z need some slight changes. 
    ms_init_z <- function(ch, f){
        for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
        states <- 3
        v <- which(ch==states) 
        ch[-v] <- NA 
        ch[v] <- states 
        return(ch) 
    }
    ## Initial values 
    initFunction <- function(){list(s = runif((dim(rCH)[2]-1), 0, 1), z = ms_init_z(rCH, f))}
    ## Parameters monitored 
    ##parameters <- c("s", "psiV", "psiF", "psiD")
    parameters <- c('a', 'b', 'c', 's')
    library(R2WinBUGS) 
    library(R2jags) 
    ## Call JAGS from R (BRT 3 min) 
    t <- system.time({ls <- jags(jags_data, initFunction, parameters, "ladyslipper.jags", n.chains = 1, n.thin = 1, n.iter = niter, n.burnin = 2000, working.directory = getwd())})
    ar <- ls$BUGSoutput$sims.array    ## using sims.array is consistent with MCMCsuite
    ar <- ar[,1,]     ## drop middle index
    ar <- ar[, -1]    ## drop deviance
    summaryStats <- c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp', 'effectiveSize')
    library(coda)
    CI95_low <- function(x) quantile(x, probs = 0.025)
    CI95_upp <- function(x) quantile(x, probs = 0.975)
    nSummaryStats <- length(summaryStats)
    nMonitorNodes <- dim(ar)[2]
    summaryArray <- array(NA, c(nSummaryStats, nMonitorNodes))
    dimnames(summaryArray) <- list(summaryStats, dimnames(ar)[[2]])
    summaryStatFunctions <- lapply(summaryStats, function(txt) eval(parse(text=txt)[[1]]))
    for(iStat in seq_along(summaryStats)) {
        summaryArray[iStat, ] <- apply(ar, 2, summaryStatFunctions[[iStat]])
    }
    theTime <- as.numeric(t[3] / 60)    ## this is consistent with MCMCsuite
    list(theTime=theTime, summary=summaryArray)
}




## orchidDHMM (multistate, from Kery & Schaub, y[i] ~ dDHMM(...) for nimble only)
## (9.7. Real data example: the showy lady's slipper)
code <- quote({
    ## -------------------------------------------------
    ## Parameters:
    ## s: survival probability 
    ## psiV: transitions from vegetative 
    ## psiF: transitions from flowering 
    ## psiD: transitions from dormant 
    ## -------------------------------------------------
    ## States (S):
    ## 1 vegetative 
    ## 2 flowering 
    ## 3 dormant 
    ## 4 dead 
    ## Observations (O):  
    ## 1 seen vegetative 
    ## 2 seen flowering 
    ## 3 not seen 
    ## -------------------------------------------------
    ## Priors and constraints 
    ## Survival: uniform 
    for (t in 1:(k-1)){  
        s[t] ~ dunif(0, 1)
    }
    ## Transitions: gamma priors 
    for (i in 1:3){
        a[i] ~ dgamma(1, 1) 
        psiD[i] <- a[i]/sum(a[1:3]) 
        b[i] ~ dgamma(1, 1) 
        psiV[i] <- b[i]/sum(b[1:3]) 
        c[i] ~ dgamma(1, 1) 
        psiF[i] <- c[i]/sum(c[1:3]) 
    }
    ## Define state-transition and observation matrices 	
    for (t in 1:(k-1)) {
        T[1,1,t] <- s[t] * psiV[1]
        T[2,1,t] <- s[t] * psiV[2]
        T[3,1,t] <- s[t] * psiV[3]
        T[4,1,t] <- 1-s[t]
        T[1,2,t] <- s[t] * psiF[1]
        T[2,2,t] <- s[t] * psiF[2]
        T[3,2,t] <- s[t] * psiF[3]
        T[4,2,t] <- 1-s[t]
        T[1,3,t] <- s[t] * psiD[1]
        T[2,3,t] <- s[t] * psiD[2]
        T[3,3,t] <- s[t] * psiD[3]
        T[4,3,t] <- 1-s[t]
        T[1,4,t] <- 0
        T[2,4,t] <- 0
        T[3,4,t] <- 0
        T[4,4,t] <- 1
    }
    T[1,1,k] <- 1
    T[2,1,k] <- 0
    T[3,1,k] <- 0
    T[4,1,k] <- 0
    T[1,2,k] <- 0
    T[2,2,k] <- 1
    T[3,2,k] <- 0
    T[4,2,k] <- 0
    T[1,3,k] <- 0
    T[2,3,k] <- 0
    T[3,3,k] <- 1
    T[4,3,k] <- 0
    T[1,4,k] <- 0
    T[2,4,k] <- 0
    T[3,4,k] <- 0
    T[4,4,k] <- 1
    for (i in 1:nind) {
        y[i, f[i]:k] ~ dDHMM(length = k-f[i]+1,
                             prior = prior[1:4],
                             condition = condition[1:3],
                             Z = Z[1:3,1:4,1:2],    useZt = 0,
                             T = T[1:4,1:4,f[i]:k], useTt = 1,
                             mult = mult[i])
    }
})

CH <- as.matrix(read.table("~/GitHub/userDistMCMC/orchids.txt", sep=" ", header = FALSE))
## first, remove all individuals not seen until the last occasion:
f <- numeric()
for (i in 1:dim(CH)[1])     f[i] <- min(which(CH[i,]!=0))
CH <- CH[which(f!=11), ]  ## remove all individuals not seen until the last occasion:
## optionally truncate data:
if(trunc) {     ind <- 1:5;     CH <- CH[ind,]     }
rCH <- CH  # Recoded CH 
rCH[rCH==0] <- 3
nind <- dim(rCH)[1]
k <- dim(rCH)[2]
f <- numeric()   ## vector with occasion of first capture 
for (i in 1:dim(CH)[1])     f[i] <- min(which(CH[i,]!=0))
Z <- array(c(1,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0), c(3,4,2))
constants <- list(f=f, k=k, nind=nind, prior=c(1/3,1/3,1/3,0), condition=c(1,1,0), mult=rep(1,nind), Z=Z)
data <- list(y = rCH)
inits <- list(s = rep(1/2,k-1), a = rep(1,3), b = rep(1,3), c = rep(1,3))

orchidDHMM <- list(code=code, constants=constants, data=data, inits=inits)







save(dipper,
     dipperDHMM,
     dipperCJS,
     ##dipperSeasonalDHMM,
     gooseDHMM,
     gooseExpJAGSfunction,
     orchidJAGSfunction,
     orchidDHMM,
     file = '~/GitHub/userDistMCMC/models.RData')







## (not using seasonal dipper model)
## dipperSeasonalDHMM (with seasonal survival, y[i] ~ dDHMM(...) for nimble only)
## code <- quote({
##     phi_flood ~ dunif(0,1)
##     phi_non   ~ dunif(0,1)
##     p         ~ dunif(0,1)
##     T[1,1,1] <- phi_non
##     T[2,1,1] <- 1 - phi_non
##     for(i in 2:3) {
##         T[1,1,i] <- phi_flood
##         T[2,1,i] <- 1 - phi_flood
##     }
##     for(i in 4:6) {
##         T[1,1,i] <- phi_non
##         T[2,1,i] <- 1 - phi_non
##     }
##     T[1,1,7] <- 1
##     T[2,1,7] <- 0
##     for(i in 1:7) {
##         T[1,2,i] <- 0
##         T[2,2,i] <- 1
##     }
##     Z[1,1,1] <- p
##     Z[2,1,1] <- 1 - p
##     Z[1,2,1] <- 0
##     Z[2,2,1] <- 1
##     Z[1:2,1:2,2] <- nimArray(0, 2, 2)
##     for (i in 1:nind)
##         y[i, first[i]:k] ~ dDHMM(length = k-first[i]+1,
##                                  prior = prior[1:2],
##                                  condition = condition[1:2],
##                                  Z = Z[1:2,1:2,1:2],        useZt = 0,
##                                  T = T[1:2,1:2,first[i]:k], useTt = 1,
##                                  mult = 1)
## })
## constants <- list(k=k, nind=nind, first=first, prior=c(1,0), condition=c(1,0))
## data      <- list(y=yDHMM)
## inits     <- list(phi_flood=0.6, phi_non=0.6, p=0.9)
## dipperSeasonalDHMM <- list(code=code, constants=constants, data=data, inits=inits)
