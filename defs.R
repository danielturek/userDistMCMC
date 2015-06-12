

dDHMM <- nimbleFunction(
    run = function(x = double(1), length = double(), pi = double(1), Z = double(2), T = double(2), condition = double(1), log.p = double()) {
        pi <- pi / sum(pi)    ## normalize prior, assume no negative values
        ##print('initial pi:');  print(pi)
        logL <- 0
        for(t in 1:length) {
            ##print('**************************************')
            ##print('t: ', t)
            ##print('pi entering into this iteration:');  print(pi)
            if(t == 1) {
                ## conditioning on first observation
                Zcond <- Z
                for(i in 1:dim(Zcond)[1])
                    for(j in 1:dim(Zcond)[2])
                        Zcond[i,j] <- Zcond[i,j] * condition[i]
                for(j in 1:dim(Zcond)[2]) {
                    s <- sum(Zcond[ ,j])
                    if(s != 0)
                        for(i in 1:dim(Zcond)[1])
                            Zcond[i,j] <- Zcond[i,j] / s
                }
                ##print('Zcond:');  print(Zcond)
                Zpi <- Zcond[x[t], ] * pi
            } else
                Zpi <- Z[x[t], ] * pi
            ##print('Zpi:');   print(Zpi)
            sumZpi <- sum(Zpi)
            logL <- logL + log(sumZpi)
            pi <- (T %*% asCol(Zpi) / sumZpi)[, 1]
            ##print('likelihood contribution: ', sumZpi)
            ##print('updated pi exiting this iteration:');   print(pi)
            ##print('**************************************')
        }
        returnType(double())
        if(log.p != 0) return(logL) else return(exp(logL))
    }
)


rDHMM <- nimbleFunction(
    run = function(n = integer(), length = double(), pi = double(1), Z = double(2), T = double(2), condition = double(1)) {
        if(n != 1) print('should only specify n=1 in rDHMM() distribution')
        print('STILL NEED TO WRITE THE rDHMM() METHOD!')
        returnType(double(1))
        return(nimVector(1, length))
    }
)


registerDistributions(list(
    dDHMM = list(
        BUGSdist = 'dDHMM(length, pi, Z, T, condition)',
        types    = c('value = double(1)', 'length = double()', 'pi = double(1)', 'Z = double(2)', 'T = double(2)', 'condition = double(1)'),
        discrete = TRUE
    )
))





