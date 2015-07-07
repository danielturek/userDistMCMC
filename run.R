
fast <- TRUE
trunc <- if(fast) TRUE else FALSE
niter <- if(fast) 5000 else 100000
setwd('~/GitHub/userDistMCMC')
source('defs.R')
source('create_data.R')
load('models.RData')

results <- resultsObjectDef(niter=niter)

results$run('dipper', dipper,     MCMCs = c('nimble','autoBlock','jags'))
results$run('dipper', dipperCJS,  MCMCs = c('nimble','autoBlock'), MCMCnames = c('nimbleCJS', 'autoBlockCJS'))
results$run('orchid', orchidDHMM, MCMCs = c('nimble','autoBlock'), MCMCnames = c('nimbleDHMM','autoBlockDHMM'), monitors = c('s','psiV','psiF','psiD'))
results$run('orchid', orchidJAGSfunction,                          MCMCnames = c('jags'))
results$run('goose',  gooseDHMM,  MCMCs = c('nimble','autoBlock'), MCMCnames = c('nimbleDHMM','autoBlockDHMM'), monitors = c('p','phi','psi'))
##out<-results$out; save(out, file = '~/GitHub/userDistMCMC/results.RData')
results$run('goose',  gooseExpJAGSfunction,                        MCMCnames = c('jags'))

out<-results$out; save(out, file = '~/GitHub/userDistMCMC/results.RData')


##results$run('dipper', dipperDHMM, MCMCs = c('nimble','autoBlock'), MCMCnames = c('nimbleDHMM','autoBlockDHMM'))
