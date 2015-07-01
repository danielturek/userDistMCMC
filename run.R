
trunc <- FALSE
niter <- 100000
setwd('~/GitHub/userDistMCMC')
source('defs.R')
source('create_data.R')
load('models.RData')

results <- resultsObjectDef(niter=niter)

results$run('dipper', dipper,     MCMCs = c('nimble','nimble_slice','jags'))
results$run('dipper', dipperDHMM, MCMCs = c('nimble','nimble_slice'), MCMCnames = c('nimbleDHMM','nimble_sliceDHMM'))
results$run('dipper', dipperCJS,  MCMCs = c('nimble','nimble_slice'), MCMCnames = c('nimbleCJS','nimble_sliceCJS'))
results$run('orchid', orchidDHMM, MCMCs = c('nimble','nimble_slice'), MCMCnames = c('nimbleDHMM','nimble_sliceDHMM'), monitors = c('s','psiV','psiF','psiD'))
results$run('orchid', orchidJAGSfunction,                             MCMCnames = c('jags'))
results$run('goose',  gooseDHMM,  MCMCs = c('nimble','nimble_slice'), MCMCnames = c('nimbleDHMM','nimble_sliceDHMM'), monitors = c('p','phi','psi'))
##results$run(dipperSeasonalDHMM, MCMCs = 'nimble')   ## stopped using

out <- results$out
##print(out)


save(out, file = '~/GitHub/userDistMCMC/results.RData')

