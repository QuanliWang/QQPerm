data(igm.data)
n.samples <- 2993
samples <-sample(n.samples)
temp <- igm.data

temp$data <- temp$data[,samples]
temp$is.case <- temp$is.case[samples,]
dim(temp$is.case) <- c(n.samples,1)

rownames(temp$is.case) <- paste('sample',1:n.samples,sep='')
colnames(temp$data) <- paste('sample',1:n.samples,sep='')


rownames(temp$data) <- paste('gene',1:18308,sep='')

example.data <- temp

n.permutations <- 1000 #too low for real analysis, default value is 1000.

#caclualte expected and observed distributions of P-values using igm data
Ps <- igm.get.pvalues(example.data$data,example.data$is.case,n.permutations)

lambda <-estlambda2(Ps$observed,Ps$perm, plot = TRUE, adjust.xy = TRUE)

