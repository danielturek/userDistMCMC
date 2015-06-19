

dDHMM <- nimbleFunction(
    run = function(x = double(1),
        length = double(), prior = double(1), condition = double(1),
        Z = double(2), Zt = double(3), useZt = double(),
        T = double(2), Tt = double(3), useTt = double(),
        log.p = double()) {
        verb <- 0
        pi <- prior / sum(prior)    ## normalize prior, assume no negative values
        if(verb) { print('initial pi:'); print(pi) }
        logL <- 0
        for(t in 1:length) {
            if(verb) { print('**************************************')
                       print('t: ', t)
                       print('pi entering into this iteration:'); print(pi) }
            if(useZt) Zcurrent <- Zt[,,t] else Zcurrent <- Z
            if(t == 1) {    ## condition on first observation
                for(i in 1:dim(Zcurrent)[1])
                    Zcurrent[i, ] <- Zcurrent[i, ] * condition[i]
                for(j in 1:dim(Zcurrent)[2]) {
                    s <- sum(Zcurrent[ ,j])
                    if(s != 0) Zcurrent[ ,j] <- Zcurrent[ ,j] / s
                }
            }
            Zpi <- Zcurrent[x[t], ] * pi
            sumZpi <- sum(Zpi)
            logL <- logL + log(sumZpi)
            if(useTt) Tcurrent <- Tt[,,t] else Tcurrent <- T
            pi <- (Tcurrent %*% asCol(Zpi) / sumZpi)[, 1]
            if(verb) { print('Zcurrent:'); print(Zcurrent)
                       print('Zpi:'); print(Zpi)
                       print('likelihood contribution: ', sumZpi)
                       print('updated pi exiting this iteration:'); print(pi)
                       print('**************************************') }
        }
        returnType(double())
        if(log.p) return(logL) else return(exp(logL))
    }
)


rDHMM <- nimbleFunction(
    run = function(n = integer(),
        length = double(), prior = double(1), condition = double(1),
        Z = double(2), Zt = double(3), useZt = double(),
        T = double(2), Tt = double(3), useTt = double()) {
        if(n != 1) print('should only specify n=1 in rDHMM() distribution')
        print('STILL NEED TO WRITE THE rDHMM() METHOD!')
        returnType(double(1))
        return(nimVector(1, length))
    }
)


registerDistributions(list(
    dDHMM = list(
        BUGSdist = 'dDHMM(length, prior, condition, Z, Zt, useZt, T, Tt, useTt)',
        types    = c('value = double(1)',
            'length = double()', 'prior = double(1)', 'condition = double(1)',
            'Z = double(2)', 'Zt = double(3)', 'useZt = double()',
            'T = double(2)', 'Tt = double(3)', 'useTt = double()'),
        discrete = TRUE
    )
))





