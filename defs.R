

library(igraph)
library(nimble)
library(coda)
library(R6)
library(abind)
library(ggplot2)



dDHMM <- nimbleFunction(
    run = function(x = double(1),
        length = double(), prior = double(1), condition = double(1),
        Z = double(3), useZt = double(),
        T = double(3), useTt = double(),
        mult = double(),
        log.p = double()) {
        verb <- 0
        pi <- prior / sum(prior)    ## normalize prior, assume no negative values
        if(verb) { print('initial pi:'); print(pi) }
        logL <- 0
        for(t in 1:length) {
            if(verb) { print('**************************************')
                       print('t: ', t)
                       print('pi entering into this iteration:'); print(pi) }
            if(useZt) Zcurrent <- Z[,,t] else Zcurrent <- Z[,,1]
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
            logL <- logL + log(sumZpi) * mult
            if(useTt) Tcurrent <- T[,,t] else Tcurrent <- T[,,1]
            pi <- (Tcurrent %*% asCol(Zpi) / sumZpi)[, 1]
            if(verb) { print('Zcurrent:'); print(Zcurrent)
                       print('Zpi:'); print(Zpi)
                       print('likelihood contribution: ', sumZpi)
                       print('contribution exponent: ', mult)
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
        Z = double(3), useZt = double(),
        T = double(3), useTt = double(),
        mult = double()) {
        if(n != 1) print('should only specify n=1 in rDHMM() distribution')
        print('STILL NEED TO WRITE THE rDHMM() METHOD!')
        returnType(double(1))
        return(nimVector(1, length))
    }
)




dCJS <- nimbleFunction(
    run = function(x = double(1),
        length = double(), last = double(), phi = double(), p = double(),
        log.p = double()) {
        L <- 1
        if(length > last)
            for(i in 1:(length-last))
                L <- 1-phi + phi*(1-p)*L
        logL <- log(L) + (last-1)*(log(phi*p))
        returnType(double())
        if(log.p) return(logL) else return(exp(logL))
    }
)

rCJS <- nimbleFunction(
    run = function(n = integer(),
        length = double(), last = double(), phi = double(), p = double()) {
        if(n != 1) print('should only specify n=1 in rCJS() distribution')
        print('STILL NEED TO WRITE THE rCJS() METHOD!')
        returnType(double(1))
        return(nimVector(1, length))
    }
)



registerDistributions(list(
    dDHMM = list(
        BUGSdist = 'dDHMM(length, prior, condition, Z, useZt, T, useTt, mult)',
        types = c('value = double(1)',
            'length = double()', 'prior = double(1)', 'condition = double(1)',
            'Z = double(3)', 'useZt = double()',
            'T = double(3)', 'useTt = double()',
            'mult = double()'),
        discrete = TRUE
    ),
    dCJS = list(
        BUGSdist = 'dCJS(length, last, phi, p)',
        types = c('value = double(1)',
            'length = double()', 'last = double()', 'phi = double()', 'p = double()'),
        discrete = TRUE
    )
))






resultsObjectDef <- function(...) resultsObjectClassDef$new(...)

resultsObjectClassDef <- R6Class(
    'resultsObjectClassDef',
    public = list(
        out = list(),
        statNames  = c('mean', 'sd', 'CI95_low', 'CI95_upp', 'effectiveSize'),
        statNames0 = c('mean', 'sd', 'CI95_low', 'CI95_upp'),
        niter = NULL,
        initialize = function(niter) {
            self$niter <- niter
        },
        run = function(name, mod, MCMCs, monitors=character(), MCMCnames) {
            if(is.list(mod)) {
                suiteOut <- MCMCsuite(
                    mod$code, mod$constants, mod$data, mod$inits,
                    MCMCs = MCMCs,
                    niter = self$niter,
                    monitors = monitors,
                    summaryStats = self$statNames,
                    makePlot = FALSE)
                if(is.null(self$out[[name]])) self$out[[name]] <- self$newListEntry(dimnames(suiteOut$summary)[[3]])
                MCMCnamesToUse <- if(!missing(MCMCnames)) MCMCnames else MCMCs
                newTimings <- suiteOut$timing[MCMCs]
                names(newTimings) <- MCMCnamesToUse
                self$out[[name]]$timing <- c(self$out[[name]]$timing, newTimings)
                newSummary <- suiteOut$summary[, self$statNames0, , drop=FALSE]
                dimnames(newSummary)[[1]] <- MCMCnamesToUse
                self$out[[name]]$summary <- abind(self$out[[name]]$summary, newSummary, along=1)
                ar <- self$out[[name]]$ESS
                arE <- self$out[[name]]$Efficiency
                newAr <- ar
                newArE <- arE
                for(mcmc in MCMCs) {
                    newAr <- rbind(newAr, suiteOut$summary[mcmc, 'effectiveSize',])
                    newArE <- rbind(newArE, suiteOut$summary[mcmc, 'effectiveSize',]/suiteOut$timing[mcmc])
                }
                dimnames(newAr) <- list(c(dimnames(ar)[[1]],MCMCnamesToUse), dimnames(ar)[[2]])
                dimnames(newArE) <- list(c(dimnames(arE)[[1]],MCMCnamesToUse), dimnames(arE)[[2]])
                self$out[[name]]$ESS <- newAr
                self$out[[name]]$Efficiency <- newArE
                return(invisible(NULL))
            } else if(is.function(mod)) {
                ret <- mod(niter=self$niter)
                MCMCnamesToUse <- if(!missing(MCMCnames)) MCMCnames else stop('must provide MCMCnames')
                self$out[[name]]$timing <- c(self$out[[name]]$timing, ret$theTime)
                names(self$out[[name]]$timing)[length(self$out[[name]]$timing)] <- MCMCnamesToUse
                alreadyParamNames <- dimnames(self$out[[name]]$summary)[[3]]
                newSummary <- ret$summary[self$statNames0,alreadyParamNames]  ## reorder
                newSummary3D <- array(as.numeric(newSummary), c(1,dim(newSummary)[1],dim(newSummary)[2]))
                dimnames(newSummary3D) <- list(MCMCnamesToUse, self$statNames0, alreadyParamNames)
                self$out[[name]]$summary <- abind(self$out[[name]]$summary, newSummary3D, along=1)
                newESS <- ret$summary['effectiveSize',]
                newESS <- newESS[alreadyParamNames]    ## reorder
                oldMCMCnames <- dimnames(self$out[[name]]$ESS)[[1]]
                self$out[[name]]$ESS <- rbind(self$out[[name]]$ESS, newESS)
                dimnames(self$out[[name]]$ESS)[[1]] <- c(oldMCMCnames, MCMCnamesToUse)
                newE <- newESS / ret$theTime   ## newESS already reordered
                self$out[[name]]$Efficiency <- rbind(self$out[[name]]$Efficiency, newE)
                dimnames(self$out[[name]]$Efficiency)[[1]] <- c(oldMCMCnames, MCMCnamesToUse)
                return(invisible(NULL))
            } else stop('bad first argument')
        },
        newListEntry = function(paramNames) {
            l <- list()
            l$timing <- numeric()
            l$summary <- array(as.numeric(NA), c(0,length(self$statNames0),length(paramNames)))
            dimnames(l$summary) <- list(NULL, self$statNames0, paramNames)
            l$ESS <- array(as.numeric(NA), c(0,length(paramNames)))
            dimnames(l$ESS) <- list(NULL, paramNames)
            l$Efficiency <- array(as.numeric(NA), c(0,length(paramNames)))
            dimnames(l$Efficiency) <- list(NULL, paramNames)
            return(l)
        }
    )
)




## Multiple plot function
##
## ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
## - cols:   Number of columns in layout
## - layout: A matrix specifying the layout. If present, 'cols' is ignored.
##
## If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
## then plot 1 will go in the upper left, 2 will go in the upper right, and
## 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    ## Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    ## If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        ## Make the panel
        ## ncol: Number of columns of plots
        ## nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots==1) {
        print(plots[[1]])
    } else {
        ## Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        ## Make each plot, in the correct location
        for (i in 1:numPlots) {
            ## Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                  layout.pos.col = matchidx$col))
        }
    }
}


results_plots <- function(out) {
    for(name in names(out)) {
        E <- out[[name]]$Efficiency
        df <- data.frame(mcmc=rep(dimnames(E)[[1]],each=dim(E)[2]), param=rep(dimnames(E)[[2]],dim(E)[1]), E=as.numeric(t(E)))
        dev.new(width=10, height=5)
        ymax <- max(df$E)
        p1 <- ggplot(df, aes(mcmc, E, fill=mcmc)) + stat_summary(fun.y='mean', geom='bar') + ggtitle(paste0(name, '\nMean')) + ylim(c(0,ymax)) + theme(legend.position='none')
        p2 <- ggplot(df, aes(mcmc, E, fill=mcmc)) + stat_summary(fun.y='min', geom='bar') + ggtitle(paste0(name, '\nMin')) + ylim(c(0,ymax)) + theme(legend.position='none')
        p3 <- ggplot(df, aes(mcmc, E, colour=mcmc)) + geom_point(size=3) + ylim(c(0,ymax)) + ggtitle(paste0(name, '\npoints')) + theme(legend.position='none')
        p4 <- ggplot(df, aes(mcmc, E, fill=mcmc)) + stat_summary(fun.y='mean', geom='bar') + ylim(c(0,ymax))
        multiplot(p1, p2, p3, p4, cols=4)
        dev.copy2pdf(file=paste0('~/GitHub/userDistMCMC/plot_', name, '.pdf'))
    }
}









## run_suite <- function(lst, monitors=character(), niter=100000, MCMCs='nimble') {
##     out <- MCMCsuite(
##         lst$code, lst$constants, lst$data, lst$inits,
##         MCMCs = MCMCs,
##         niter = niter,
##         monitors = monitors,
##         summaryStats = c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp', 'effectiveSize'),
##         makePlot = FALSE)
##     message('> timing:')
##     print(out$timing)
##     message('> summary statistics:')
##     print(out$summary[, c('mean','sd','CI95_low','CI95_upp'), ])
##     message('> effectiveSize:')
##     print(out$summary[, 'effectiveSize', ])
##     message('> effectiveSize / timing:')
##     print(out$summary[, 'effectiveSize', ] / out$timing['nimble'])
##     return(invisible(out))
## }
