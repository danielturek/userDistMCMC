

setwd('~/GitHub/userDistMCMC')
source('defs.R')
load('~/GitHub/userDistMCMC/results.RData')
##load('~/GitHub/userDistMCMC/resultsNew.RData')

results_check(out)   ### eventually can DELETE this line

results_plot(out)   ### eventually can DELETE this line


##results_check(results$out)   ### eventually use this line instead

##results_plot(results$out)   ### eventually use this line instead

