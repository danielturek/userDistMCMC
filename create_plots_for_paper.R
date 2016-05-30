## NEW FIGRURES for revision of paper for EES journal
## this will also include scatterplot or boxplot of values
library(ggplot2)
library(dplyr)
load('~/github/userDistMCMC/results.RData')
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    if (is.null(layout)) layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
    if (numPlots==1) { print(plots[[1]])
    } else { grid.newpage()
             pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
             for (i in 1:numPlots) {
                 matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                 print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
             }
         }
}
df <- results$df
df$mcmc <- as.character(df$mcmc)   ## make 'mcmc' column as character strings
df[df$model=='dipper' & df$mcmc=='nimble', ]$mcmc <- 'Latent State'
df[df$model=='dipper' & df$mcmc=='jags', ]$mcmc <- 'Latent State (JAGS)'
df[df$model=='dipper' & df$mcmc=='nimbleCJS2', ]$mcmc <- 'Filter'
df[df$model=='dipper' & df$mcmc=='jagsPoisson', ]$mcmc <- 'Filter (JAGS)'
df[df$model=='orchid' & df$mcmc=='jags', ]$mcmc <- 'Latent State (JAGS)'
df[df$model=='orchid' & df$mcmc=='nimbleDHMM2', ]$mcmc <- 'Filter'
df[df$model=='orchid' & df$mcmc=='autoBlockDHMM2', ]$mcmc <- 'Filter & Block'
df[df$model=='goose'  & df$mcmc=='jagsExp', ]$mcmc <- 'Latent State (JAGS)'
df[df$model=='goose'  & df$mcmc=='nimbleDHMM', ]$mcmc <- 'Filter (RR)'
df[df$model=='goose'  & df$mcmc=='autoBlockDHMM', ]$mcmc <- 'Filter & Block (RR)'
dipperMCMCs <- c('Latent State', 'Latent State (JAGS)', 'Filter', 'Filter (JAGS)')
orchidMCMCs <- c('Latent State (JAGS)', 'Filter', 'Filter & Block')
gooseMCMCs <- c('Latent State (JAGS)', 'Filter (RR)', 'Filter & Block (RR)')
dipperWidth <- 6.7
orchidWidth <- 6.7
gooseWidth  <- 6.7
dipperBarWidth <- 0.9
orchidBarWidth <- 0.8
gooseBarWidth <- 0.8
makePaperPlot2 <- function(model_arg, mcmcs, thiswidth, thisBarWidth, thisdf) {
    thisdf <- filter(thisdf, model == model_arg)
    thisdf$mcmc <- as.factor(thisdf$mcmc)
    thisdf <- filter(thisdf, mcmc %in% mcmcs)
    thisdf$mcmc <- factor(as.character(thisdf$mcmc), levels=mcmcs)   ## re-order, according to mcmcs argument
    p1 <- ggplot(thisdf, aes(mcmc, Efficiency, fill=mcmc), xlab='') + stat_summary(fun.y='min', geom='bar', width=thisBarWidth) + ggtitle('Minimum') + ylab('Sampling efficiency (ESPS)') + theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank())
    p2 <- ggplot(thisdf, aes(mcmc, Efficiency, fill=mcmc), xlab='') + stat_summary(fun.y='mean', geom='bar', width=thisBarWidth) + ggtitle('Mean') + ylab('Sampling efficiency (ESPS)') + theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank())
    ##p3 <- ggplot(thisdf, aes(mcmc, Efficiency, colour=mcmc), xlab='') + geom_point(size=1) + ggtitle('Values') + ylab('Sampling efficiency (ESPS)') + theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank())
    p3 <- ggplot(thisdf, aes(mcmc, Efficiency, colour=mcmc), xlab='') + geom_boxplot(width=0.6) + ggtitle('Boxplot') + ylab('Sampling efficiency (ESPS)') + theme(legend.position='none', axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank())
    p4 <- ggplot(thisdf, aes(mcmc, Efficiency, fill=mcmc)) + geom_bar(stat='identity') + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = 'left') + labs(fill='MCMC Algorithm     ')
    dev.new(width=thiswidth, height=3)
    multiplot(p1, p2, p3, p4, cols=4)
    localFile <- paste0('~/github/userDistMCMC/plot_', model_arg, '.pdf')
    paperFile <- paste0('~/github/nimble/nimblePapers/CR-MCMC/EES_revision/plot_', model_arg, '.pdf')
    dev.copy2pdf(file=localFile)
    system(paste0('cp ', localFile, ' ', paperFile))
}
makePaperPlot2('dipper', dipperMCMCs, dipperWidth, dipperBarWidth, df)
makePaperPlot2('orchid', orchidMCMCs, orchidWidth, orchidBarWidth, df)
makePaperPlot2('goose',  gooseMCMCs,  gooseWidth,  gooseBarWidth, df)


## original paper plots for first submissions to Biometrics, EES, etc.
library(ggplot2)
library(dplyr)
load('~/github/userDistMCMC/results.RData')
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
    localFile <- paste0('~/github/userDistMCMC/plot_', model_arg, '.pdf')
    paperFile <- paste0('~/github/nimble/nimblePapers/CR-MCMC/plot_', model_arg, '.pdf')
    dev.copy2pdf(file=localFile)
    system(paste0('cp ', localFile, ' ', paperFile))
    ## adding eps files, for Biometrics submission
    localFile <- paste0('~/github/userDistMCMC/plot_', model_arg, '.eps')
    paperFile <- paste0('~/github/nimble/nimblePapers/CR-MCMC/plot_', model_arg, '.eps')
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


## NEW figures for presentations on CR-MCMC
## for Dipper model, only showing Latent State and Filtering MCMCs
## THIS WILL OVERWRITE 'df' object !!!!!
model_arg <- 'dipper'
mcmcs <- c('Latent State', 'Filtering')
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
dev.copy2pdf(file='~/Downloads/dipper_samplingEfficiency.pdf')
