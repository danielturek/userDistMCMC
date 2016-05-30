

setwd('~/github/userDistMCMC')
source('defs.R')
load('~/github/userDistMCMC/results.RData')
##load('~/github/userDistMCMC/resultsNew.RData')


results$check()
results$quickplot()




## migrating to new results
setwd('~/github/userDistMCMC')
source('defs.R')
load('~/github/userDistMCMC/resultsSave.RData')
resultsOld <- results
results <- resultsObjectDef(niter=resultsOld$niter)
results$out <- resultsOld$out
results$processOutIntoDF()
save(results, file='~/github/userDistMCMC/results.RData')






