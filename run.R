
fast <- FALSE
trunc <- if(fast) TRUE else FALSE
niter <- if(fast) 5000 else 100000
setwd('~/GitHub/userDistMCMC')
source('defs.R')
source('create_data.R')
load('models.RData')

results <- resultsObjectDef(niter=niter)

##results$run('dipper', dipper,     MCMCs = c('nimble','autoBlock','jags'))
##save(results, file='~/GitHub/userDistMCMC/resultsNew.RData')

##results$run('dipper', dipperCJS,  MCMCs = c('nimble','autoBlock'), MCMCnames = c('nimbleCJS', 'autoBlockCJS'))
##save(results, file='~/GitHub/userDistMCMC/resultsNew.RData')

##results$run('orchid', orchidDHMM, MCMCs = c('nimble','autoBlock'), MCMCnames = c('nimbleDHMM','autoBlockDHMM'))
results$run('orchid', orchidDHMM, MCMCs = c('nimble','nimble_slice'), MCMCnames = c('nimbleDHMM','sliceDHMM'))
save(results, file='~/GitHub/userDistMCMC/resultsNew.RData')

results$run('orchid', orchidJAGSfunction,                          MCMCnames = c('jags'))
save(results, file='~/GitHub/userDistMCMC/resultsNew.RData')

##results$run('goose',  gooseDHMM,  MCMCs = c('nimble','autoBlock'), MCMCnames = c('nimbleDHMM','autoBlockDHMM'))
##save(results, file='~/GitHub/userDistMCMC/resultsNew.RData')

##results$run('goose',  gooseExpJAGSfunction,                        MCMCnames = c('jagsExp'))
##save(results, file='~/GitHub/userDistMCMC/resultsNew.RData')


##results$run('dipper', dipperDHMM, MCMCs = c('nimble','autoBlock'), MCMCnames = c('nimbleDHMM','autoBlockDHMM'))













