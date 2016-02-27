#' Read sample file and genotype collipsing matrix in IGM format.
#'
#' @author Slave Petrovski and Quanli Wang
#' @param sample.file The input sample file.
#' @param matrix.file The input matrix file.
#' @param filter.list A list of genes that will be excluded from the analysis.
#'
#' @return Returns A list contains data matrix (matrix) and phenotype indicators (is.case).
#'
#' @examples
#' #igm.data <- igm.read.data(sample.file, matrix.file)
#'
igm.read.data <- function(sample.file, matrix.file, filter.list = NULL) {
  #get sample and phenotype info
  ped <- read.table(sample.file,header = FALSE, sep = '\t', as.is = TRUE)
  samples <- ped[,2]
  is.case <- as.matrix(ped[,6] ==2)
  rownames(is.case) <- samples

  #get count matrix
  data <- read.table(matrix.file,header = TRUE, sep = '\t', as.is = TRUE)
  genes <- data[,1]
  data <- data[,samples]
  data[data > 0] <- 1

  #convert data to matrix
  matrix <- as.matrix(data[])
  rownames(matrix)<-genes

  if (!is.null(filter.list)) {
    genes.common <- intersect(genes,filter.list)
    matrix <- matrix[genes.common,]
  }
  out <- list()
  out$data <- matrix
  out$is.case <- is.case
  out
}

#' Generate NULL distribution of P-values through label switching and permutation and compute distribution of observed P-values.
#' @author Slave Petrovski and Quanli Wang
#' @param matrix The input genotype matrix, with rows for genes and columns for samples.
#' @param is.case The case/control indicator.
#' @param n.permutations Number of label permutaitons.
#'
#' @return Returns A list contains observed P-values (observed) and permuted P-Values for NULL distribution (perm).
#'
#' @examples
#' #Ps  <- igm.get.pvalues(matrix, is.case)
#'
igm.get.pvalues <-function(matrix, is.case, n.permutations = 1000) {
  #Number of cases and controls
  n.samples <- length(is.case)
  n.cases <- sum(is.case)
  n.controls <- n.samples - n.cases

  #pre-compute all possible contingency.table for perofrmance
  #this will create a look up table
  contingency.table = matrix(0,2,2)
  Fisher.precompute <- matrix(0, n.cases + 1, max(rowSums(matrix)) + 1)
  for (i in 1: dim(Fisher.precompute)[1]) {
    for (j in 1:dim(Fisher.precompute)[2]) {
      contingency.table[1,1] <- i-1
      contingency.table[1,2] <- n.cases - contingency.table[1,1]
      contingency.table[2,1] <- j -1
      contingency.table[2,2] <- n.controls - contingency.table[2,1]
      Fisher.precompute[i,j] <- fisher.test(contingency.table)$p.value
    }
  }

  #permutaiton, save all p-values just in case median will be needed later on
  P.Values <- matrix(1,dim(matrix)[1],n.permutations)
  total.1 <- rowSums(matrix)
  for (i in 1: n.permutations) {
    K <- sample.int(n.samples, size = n.cases, replace = FALSE)
    Labels.1.1 <- rowSums(matrix[,K])
    Labels.0.1 <- total.1 - Labels.1.1
    P.Values[,i] <- sort(Fisher.precompute[cbind(Labels.1.1+1,Labels.0.1+1)])
  }
  P.perm <- rowMeans(P.Values)


  #compute observed p values
  K <- which(is.case)
  Labels.1.1 <- rowSums(matrix[,K])
  Labels.0.1 <- total.1 - Labels.1.1
  P.observed <- sort(Fisher.precompute[cbind(Labels.1.1+1,Labels.0.1+1)])

  out <- list()
  out$perm <- P.perm
  out$observed <- P.observed
  out
}
