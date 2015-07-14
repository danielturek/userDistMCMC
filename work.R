

setwd('~/GitHub/userDistMCMC')
source('defs.R')
load('~/GitHub/userDistMCMC/results.RData')
##load('~/GitHub/userDistMCMC/resultsGandalf.RData')
##load('~/GitHub/userDistMCMC/resultsNew.RData')


results_check(results$out)
results_plot(results$out)



## rm(list=ls())
## a <- 2
## source('defs.R')
## if(a==1) load('~/GitHub/userDistMCMC/results.RData')
## if(a==2) load('~/GitHub/userDistMCMC/resultsGandalf.RData')
## lapply(results$out, function(a) a$timing)
