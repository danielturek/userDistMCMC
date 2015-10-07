

library(ggplot2)
library(dplyr)
load('~/GitHub/userDistMCMC/results.RData')
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
makePaperPlot <- function(model_arg, mcmcs, df = results$df) {
    ##browser()
    df <- filter(df, model == model_arg)
    df$mcmc <- as.character(df$mcmc)
    df <- df[df$mcmc != 'nimbleCJS',]               ## remove the nimble CJS alg,
    df$mcmc[df$mcmc=='nimbleCJS2'] <- 'nimbleCJS'   ## then rename CJS2 to CJS
    df$mcmc[df$mcmc=='jagsExp'] <- 'jags'
    df$mcmc <- as.factor(df$mcmc)
    df <- filter(df, mcmc %in% mcmcs)
    newDF <- data.frame(mcmc = character(), Mean = numeric(), Minimum = numeric(), stringsAsFactors=FALSE)
    for(mcmc in unique(df$mcmc)) {
        thismean <- mean(df[df$mcmc==mcmc, 'Efficiency'])
        thismin <- min(df[df$mcmc==mcmc, 'Efficiency'])
        newDF <- rbind(newDF, data.frame(mcmc = mcmc, Mean = thismean, Minimum = thismin, stringsAsFactors=FALSE))
    }
    df <- newDF
    df$mcmc <- factor(as.character(df$mcmc), levels=mcmcs)   ## re-order, according to mcmcs argument
    ##ymax <- max(c(df$Mean, df$Minimum))
    p1 <- ggplot(df, aes(mcmc, Minimum, fill=mcmc)) + geom_bar(stat='identity') + theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank()) + ylab('Minimum sampling efficiency')
    p2 <- ggplot(df, aes(mcmc, Mean,    fill=mcmc)) + geom_bar(stat='identity') + theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank()) + ylab('Mean sampling efficiency')
    p4 <- ggplot(df, aes(mcmc, Mean, fill=mcmc)) + geom_bar(stat='identity') + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'left') + labs(fill='MCMC Algorithm      ')
    dev.new(width=6, height=3)
    multiplot(p1, p2,    p4, cols=3)
    localFile <- paste0('~/GitHub/userDistMCMC/plot_', model_arg, '.pdf')
    paperFile <- paste0('~/GitHub/nimble/nimblePapers/CR-MCMC/plot_', model_arg, '.pdf')
    dev.copy2pdf(file=localFile)
    system(paste0('cp ', localFile, ' ', paperFile))
}

makePaperPlot(model_arg='dipper', mcmcs = c('jags', 'jagsPoisson', 'nimble', 'nimbleCJS'))
makePaperPlot(model_arg='orchid', mcmcs = c('jags', 'nimbleDHMM2', 'autoBlockDHMM2'))
makePaperPlot(model_arg='goose', mcmcs = c('jags', 'nimbleDHMM', 'autoBlockDHMM'))



