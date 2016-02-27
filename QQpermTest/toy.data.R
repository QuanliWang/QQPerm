data(igm.data)
n.samples <- 2803
samples <-sample.int(2993,n.samples)
temp <- igm.data
temp$data <- temp$data[,samples]
temp$is.case <- temp$is.case[samples,]
dim(temp$is.case) <- c(n.samples,1)
row.names(temp$is.case) <- 1:n.samples
colnames(temp$data) <- 1:n.samples


example.data <- temp
rownames(example.data$data) <- 1:18308

n.permutations <- 5 #too low for real analysis, default value is 1000.

#caclualte expected and observed distributions of P-values using igm data
Ps <- igm.get.pvalues(toy.data$data,toy.data$is.case,n.permutations)


