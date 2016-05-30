
fast <- FALSE
trunc <- if(fast) TRUE else FALSE
niter <- if(fast) 5000 else 100000
setwd('~/github/userDistMCMC')
source('defs.R')
source('create_data.R')
load('models.RData')

results <- resultsObjectDef(niter=niter)

results$run('dipper', dipper,     MCMCs = c('nimble','nimble_slice','autoBlock','jags'))
save(results, file='~/github/userDistMCMC/resultsNew.RData')

results$run('dipper', dipperCJS,  MCMCs = c('nimble','nimble_slice','autoBlock'), MCMCnames = c('nimbleCJS','sliceCJS','autoBlockCJS'))
save(results, file='~/github/userDistMCMC/resultsNew.RData')

results$run('dipper', dipperCJS2,  MCMCs = c('nimble','nimble_slice','autoBlock'), MCMCnames = c('nimbleCJS2','sliceCJS2','autoBlockCJS2'))
save(results, file='~/github/userDistMCMC/resultsNew.RData')

##results$run('dipper', dipperChiJAGSfunction,                       MCMCnames = c('jagsChi'))
##save(results, file='~/github/userDistMCMC/resultsNew.RData')

results$run('dipper', dipperPoissonJAGSfunction,                   MCMCnames = c('jagsPoisson'))
save(results, file='~/github/userDistMCMC/resultsNew.RData')

results$run('orchid', orchidDHMM, MCMCs = c('nimble','autoBlock'), MCMCnames = c('nimbleDHMM','autoBlockDHMM'))
save(results, file='~/github/userDistMCMC/resultsNew.RData')

results$run('orchid', orchidDHMM2, MCMCs = c('nimble','autoBlock'), MCMCnames = c('nimbleDHMM2','autoBlockDHMM2'))
save(results, file='~/github/userDistMCMC/resultsNew.RData')

results$run('orchid', orchidJAGSfunction,                          MCMCnames = c('jags'))
save(results, file='~/github/userDistMCMC/resultsNew.RData')

results$run('goose',  gooseDHMM,  MCMCs = c('nimble','autoBlock'), MCMCnames = c('nimbleDHMM','autoBlockDHMM'))
save(results, file='~/github/userDistMCMC/resultsNew.RData')

results$run('goose',  gooseExpJAGSfunction,                        MCMCnames = c('jagsExp'))
save(results, file='~/github/userDistMCMC/resultsNew.RData')

## out of ordering: slice sampling for 'orchid' model, which will take forever:
results$run('orchid', orchidDHMM, MCMCs = c('nimble_slice'),       MCMCnames = c('sliceDHMM'))
save(results, file='~/github/userDistMCMC/resultsNew.RData')

## out of ordering: slice sampling for 'orchid' model, which will take forever:
results$run('orchid', orchidDHMM2, MCMCs = c('nimble_slice'),       MCMCnames = c('sliceDHMM2'))
save(results, file='~/github/userDistMCMC/resultsNew.RData')



##results$run('orchid', orchidDHMM,  MCMCs = c('nimble'), MCMCnames = c('nimbleDHMM'))
##results$run('orchid', orchidJAGSfunction,               MCMCnames = c('jags'))
##results$run('orchid', orchidDHMM2, MCMCs = c('nimble'), MCMCnames = c('nimbleDHMM2'))

##results$check()
##results$quickplot()

