
rrr

CH <- as.matrix(read.table("~/github/userDistMCMC/orchids.txt", sep=" ", header = FALSE))

f <- numeric()
for (i in 1:dim(CH)[1])
    f[i] <- min(which(CH[i,]!=0))
CH <- CH[which(f!=11), ]
rCH <- CH
rCH[rCH==0] <- 3
nind <- dim(rCH)[1]
k <- dim(rCH)[2]
f <- numeric()
for (i in 1:dim(CH)[1])
    f[i] <- min(which(CH[i,]!=0))

y <- rCH

save(list = c('y', 'nind', 'k', 'f'), file='~/github/userDistMCMC/orchidData.RData')

