

setwd('~/GitHub/userDistMCMC')
source('defs.R')
load('~/GitHub/userDistMCMC/results.RData')
##load('~/GitHub/userDistMCMC/resultsNew.RData')


results_check(results$out)
results_plot(results$out)

