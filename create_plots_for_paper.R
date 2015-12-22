

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
makePaperPlot <- function(model_arg, mcmcs, df) {
    df <- filter(df, model == model_arg)
    df$mcmc <- as.factor(df$mcmc)
    df <- filter(df, mcmc %in% mcmcs)
    df$mcmc <- factor(as.character(df$mcmc), levels=mcmcs)   ## re-order, according to mcmcs argument
    p1 <- ggplot(df, aes(mcmc, Minimum, fill=mcmc)) + geom_bar(stat='identity') + theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank()) + ylab('Sampling efficiency (ESPS)') + ggtitle('Minimum')
    p2 <- ggplot(df, aes(mcmc, Mean,    fill=mcmc)) + geom_bar(stat='identity') + theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank()) + ylab('Sampling efficiency (ESPS)') + ggtitle('Mean')
    p4 <- ggplot(df, aes(mcmc, Mean, fill=mcmc)) + geom_bar(stat='identity') + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'left') + labs(fill='MCMC Algorithm      ')
    wid <- if(model_arg=='goose') 6.5 else 6
    dev.new(width=wid, height=3)
    multiplot(p1, p2,    p4, cols=3)
    localFile <- paste0('~/GitHub/userDistMCMC/plot_', model_arg, '.pdf')
    paperFile <- paste0('~/GitHub/nimble/nimblePapers/CR-MCMC/plot_', model_arg, '.pdf')
    dev.copy2pdf(file=localFile)
    system(paste0('cp ', localFile, ' ', paperFile))
    ## adding eps files, for Biometrics submission
    localFile <- paste0('~/GitHub/userDistMCMC/plot_', model_arg, '.eps')
    paperFile <- paste0('~/GitHub/nimble/nimblePapers/CR-MCMC/plot_', model_arg, '.eps')
    dev.copy2eps(file=localFile)
    system(paste0('cp ', localFile, ' ', paperFile))
}

df <- results$df
newDF <- data.frame(model = character(), mcmc = character(), Mean = numeric(), Minimum = numeric(), stringsAsFactors=FALSE)
for(mod in unique(df$model)) {
    for(mcmc in unique(df[df$model==mod,]$mcmc)) {
        thismean <- mean(df[df$model==mod & df$mcmc==mcmc, 'Efficiency'])
        thismin <- min(df[df$model==mod & df$mcmc==mcmc, 'Efficiency'])
        newDF <- rbind(newDF, data.frame(model = mod, mcmc = mcmc, Mean = thismean, Minimum = thismin, stringsAsFactors=FALSE))
    }
}
df <- newDF
df$mcmc <- as.character(df$mcmc)   ## make 'mcmc' column as character strings
df[df$model=='dipper' & df$mcmc=='nimble', ]$mcmc <- 'Latent State'
df[df$model=='dipper' & df$mcmc=='jags', ]$mcmc <- 'Latent State (JAGS)'
df[df$model=='dipper' & df$mcmc=='nimbleCJS2', ]$mcmc <- 'Filtering'
df[df$model=='dipper' & df$mcmc=='jagsPoisson', ]$mcmc <- 'Filtering (JAGS)'
df[df$model=='orchid' & df$mcmc=='jags', ]$mcmc <- 'Latent State (JAGS)'
df[df$model=='orchid' & df$mcmc=='nimbleDHMM2', ]$mcmc <- 'Filtering'
df[df$model=='orchid' & df$mcmc=='autoBlockDHMM2', ]$mcmc <- 'Filtering & Blocking'
df[df$model=='goose'  & df$mcmc=='jagsExp', ]$mcmc <- 'Latent State (JAGS)'
df[df$model=='goose'  & df$mcmc=='nimbleDHMM', ]$mcmc <- 'Filtering (RR)'
df[df$model=='goose'  & df$mcmc=='autoBlockDHMM', ]$mcmc <- 'Filtering & Blocking (RR)'

makePaperPlot('dipper', c('Latent State', 'Latent State (JAGS)', 'Filtering', 'Filtering (JAGS)'), df)
makePaperPlot('orchid', c('Latent State (JAGS)', 'Filtering', 'Filtering & Blocking'), df)
makePaperPlot('goose',  c('Latent State (JAGS)', 'Filtering (RR)', 'Filtering & Blocking (RR)'), df)



