#' An example collapsing data matrx and associated case/control status
#'
#' @format A list contains both genotype data matrix and case/control status:
#'
#'  \describe{
#'     \item{data}{a named data matrix, with rows for genes and columns for samples. Each cell in the matrix represent the counts of variants for given (gene,sample) pair.}
#'    \item{is.case}{case/control indicator. The order of indicators must matches that of the samples defined in data matrix.}
#' }
#' @source \url{http://ask.slave.for.it/}
"example.data"


#' An example collapsing data matrx and associated case/control status
#'
#' @format A list contains both genotype data matrix and case/control status:
#'
#'  \describe{
#'     \item{data}{a named data matrix, with rows for genes and columns for samples. Each cell in the matrix represent the counts of variants for given (gene,sample) pair.}
#'    \item{is.case}{case/control indicator. The order of indicators must matches that of the samples defined in data matrix.}
#' }
#' @source \url{http://ask.slave.for.it/}
"toy.data"


#' Distributions of expected and observed P-values from igm.data data set.
#'
#' @format A list contains both expected and observed P-values.
#'
#'  \describe{
#'     \item{perm}{A vector for expected P-values.}
#'    \item{observed}{A vector for observed P-values.}
#' }
#' @source \url{http://ask.slave.for.it/}
"example.Ps"
