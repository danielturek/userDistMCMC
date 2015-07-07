
fast <- TRUE
trunc <- if(fast) TRUE else FALSE
niter <- if(fast) 5000 else 100000
setwd('~/GitHub/userDistMCMC')
source('defs.R')
source('create_data.R')
load('models.RData')

results <- resultsObjectDef(niter=niter)

results$run('dipper', dipper,     MCMCs = c('nimble','nimble_slice','jags','autoBlock'))
results$run('dipper', dipperDHMM, MCMCs = c('nimble','nimble_slice','autoBlock'), MCMCnames = c('nimbleDHMM','nimble_sliceDHMM','autoBlockDHMM'))
results$run('dipper', dipperCJS,  MCMCs = c('nimble','nimble_slice'), MCMCnames = c('nimbleCJS','nimble_sliceCJS'))
results$run('orchid', orchidDHMM, MCMCs = c('nimble','nimble_slice','autoBlock'), MCMCnames = c('nimbleDHMM','nimble_sliceDHMM','autoBlockDHMM'), monitors = c('s','psiV','psiF','psiD'))
results$run('orchid', orchidJAGSfunction,                             MCMCnames = c('jags'))
results$run('goose',  gooseDHMM,  MCMCs = c('nimble','nimble_slice','autoBlock'), MCMCnames = c('nimbleDHMM','nimble_sliceDHMM','autoBlockDHMM'), monitors = c('p','phi','psi'))

out <- results$out
save(out, file = '~/GitHub/userDistMCMC/results.RData')

results$run('goose',  gooseExpJAGSfunction,                           MCMCnames = c('jags'))

out <- results$out
save(out, file = '~/GitHub/userDistMCMC/results.RData')


##results$run(dipperSeasonalDHMM, MCMCs = 'nimble')   ## stopped using


