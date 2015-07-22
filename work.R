

setwd('~/GitHub/userDistMCMC')
source('defs.R')
load('~/GitHub/userDistMCMC/results.RData')
##load('~/GitHub/userDistMCMC/resultsNew.RData')


results$check()
results$plot()




## migrating to new results
setwd('~/GitHub/userDistMCMC')
source('defs.R')
load('~/GitHub/userDistMCMC/resultsSave.RData')
resultsOld <- results
niter <- resultsOld$niter
results <- resultsObjectDef(niter=niter)
results$out <- resultsOld$out
results$processOutIntoDF()
save(results, file='~/GitHub/userDistMCMC/results.RData')






