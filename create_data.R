
############################
## dipper
############################

load('~/GitHub/legacy/dipper/dipperData.RData')
## optionally truncate data:
#ind <- c(1:3);   nind<-length(ind);   first<-first[ind];   y<-y[ind,,drop=FALSE];   yDHMM<-yDHMM[ind,,drop=FALSE];   x_init<-x_init[ind,,drop=FALSE]
yDHMM <- 2 - y

## regular dipper model
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

## DHMM dipper model
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
    for (i in 1:nind)
        y[i, first[i]:k] ~ dDHMM(length=k-first[i]+1, prior=prior[1:2], condition=condition[1:2],
                                 Z=Z[1:2,1:2,1:2], useZt=0,
                                 T=T[1:2,1:2,1:2], useTt=0,
                                 mult=1)
})
constants <- list(k=k, nind=nind, first=first, prior=c(1,0), condition=c(1,0))
data      <- list(y=yDHMM)
inits     <- list(phi=0.6, p=0.9)
dipperDHMM <- list(code=code, constants=constants, data=data, inits=inits)

## dipper with seasonal variation in phi (flood, non-flood)
code <- quote({
    phi_flood ~ dunif(0,1)
    phi_non   ~ dunif(0,1)
    p         ~ dunif(0,1)
    T[1,1,1] <- phi_non
    T[2,1,1] <- 1 - phi_non
    for(i in 2:3) {
        T[1,1,i] <- phi_flood
        T[2,1,i] <- 1 - phi_flood
    }
    for(i in 4:6) {
        T[1,1,i] <- phi_non
        T[2,1,i] <- 1 - phi_non
    }
    T[1,1,7] <- 1
    T[2,1,7] <- 0
    for(i in 1:7) {
        T[1,2,i] <- 0
        T[2,2,i] <- 1
    }
    Z[1,1,1] <- p
    Z[2,1,1] <- 1 - p
    Z[1,2,1] <- 0
    Z[2,2,1] <- 1
    Z[1:2,1:2,2] <- nimArray(0, 2, 2)
    for (i in 1:nind)
        y[i, first[i]:k] ~ dDHMM(length=k-first[i]+1, prior=prior[1:2], condition=condition[1:2],
                                 Z=Z[1:2,1:2,1:2],        useZt=0,
                                 T=T[1:2,1:2,first[i]:k], useTt=1,
                                 mult=1)
})
constants <- list(k=k, nind=nind, first=first, prior=c(1,0), condition=c(1,0))
data      <- list(y=yDHMM)
inits     <- list(phi_flood=0.6, phi_non=0.6, p=0.9)
dipperSeasonalDHMM <- list(code=code, constants=constants, data=data, inits=inits)



## canadian good CR data from "handbook of CR analyses"
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
y[which(y == 0)] <- 4
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
first <- apply(y, 1, function(hist) which(hist!=4)[1])
nind <- dim(y)[1]
k <- dim(y)[2]

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
        y[i, first[i]:k] ~ dDHMM(length=k-first[i]+1, prior=prior[1:k], condition=condition[1:k],
                                 Z=Z[1:k,1:k,first[i]:k], useZt=1,
                                 T=T[1:k,1:k,first[i]:k], useTt=1,
                                 mult=mult[i])
})
constants <- list(nind=nind, k=k, first=first, mult=mult,
                  prior=c(1/3, 1/3, 1/3, 0),
                  condition = c(1, 1, 1, 0))
data <- list(y=y)
inits <- list(p=rep(1/2,6), phi=rep(1/2,3), alpha=array(0,c(2,3,2)))
goose <- list(code=code, constants=constants, data=data, inits=inits)



save(dipper,
     dipperDHMM,
     dipperSeasonalDHMM,
     goose,
     file = '~/GitHub/userDistMCMC/models.RData')












