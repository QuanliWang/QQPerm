library(QQperm)
sample.file <- '/Users/qw2192/Desktop/w2015/QQPlot_byPermuation/SUDEP_PEDMAP_v5_2936ctrls.txt'
input.matrix.file <-'/Users/qw2192/Desktop/w2015/QQPlot_byPermuation/v5_CCDS_noben_0005_ExACEVS0005_exonprune_2936controls_gene.sample.matrix.txt'
n.permutations <- 100

## compute observed and permuted p-values
data <- igm.read.data(sample.file, matrix.file)

Ps <- igm.get.pvalues(data$data,data$is.case,n.permutations)

pdf("QQ_output.pdf")
qqplot(Ps$perm, Ps$observed)
lambda <-estlambda2(Ps$observed,Ps$perm, plot = TRUE, adjust.xy = TRUE)
dev.off()
