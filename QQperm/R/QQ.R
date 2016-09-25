#' @title Estimate the inflation factor for a distribution of P-values or 1df chi-square test.
#' @description Estimate the inflation factor for a distribution of P-values or 1df chi-square test using a NULL distribution of P-values provided by user.
#'  The implementation is similar to the default implementation in GenABEL but does not assume the NULL distribution of P-values to be uniform.
#' @author Slave Petrovski and Quanli Wang
#' @param p.o Observed P-values from the data.
#' @param p.e Expected P-values from the NULL distribution, usually obtained through label permutations of the data.
#' @param plot Indicate if a plot should be produced.
#' @param filter Indicate if the filter should be applied. This parameter behaves the same as in estlambda in GenABEL.
#' @param adjust.xy Indicate if the x-axis and y-axis should be adjusted to their own range.
#' @param ... arguments passed the plot function.
#'
#' @return Returns a list containing the estimated lambda value (estimate) and its standard deviation (se).
#'
#' @examples
#' #load pre-computed p-values for igm dataset
#' library(QQperm)
#' data("example.Ps")
#'
#' #print output to pdf file only if not running in interactive mode
#' if (!interactive()) {
#'  pdf("lambda.pdf")
#' }
#'
#' #estimate inflation factor and generate plot.
#' lambda <-estlambda2(example.Ps$observed,example.Ps$perm, plot = TRUE, adjust.xy = TRUE)
#'
#' if (!interactive()) {
#'  dev.off()
#' }

estlambda2 <-function (p.o, p.e, plot = FALSE, filter = TRUE, adjust.xy = FALSE, ...) {
  p.o <- p.o[which(!is.na(p.o))]
  p.e <- p.e[which(!is.na(p.e))]
  ntp <-length(p.o)
  if (ntp != length(p.e)) {
    stop("data does not match")
  }
  if (ntp == 1) {
    warning(paste("One measurement, lambda = 1 returned"))
    return(list(estimate = 1, se = 999.99))
  }
  if (ntp < 10)
    warning(paste("number of points is too small:", ntp))
   p.o[p.o>1] = 1
   p.e[p.e>1] = 1

   p.o <- qchisq(p.o, 1, lower.tail = FALSE)
   p.e <- qchisq(p.e, 1, lower.tail = FALSE)

  if (filter) {
    to.be.removed <- which((abs(p.o) < 1e-08) | (abs(p.e) < 1e-08))
    p.o[to.be.removed] <- NA
    p.e[to.be.removed] <- NA
  }

  p.o <- sort(p.o)
  p.e <- sort(p.e)

  out <- list()

  s <- summary(lm(p.o ~ 0 + p.e))$coeff
  out$estimate <- s[1, 1]
  out$se <- s[1, 2]

  if (plot) {
    lim <- c(0, max(p.o, p.e, na.rm = TRUE))
    #oldmargins <- par()$mar
    #par(mar = oldmargins + 0.2, cex=1.5)
    if (adjust.xy) {
      xlim <- c(0,round(max(p.e) + 1))
      ylim <- c(0,round(max(p.o) + 1))
    } else {
      max_lim <- max(max(p.e), max(p.o)) + 1
      xlim <- c(0, max_lim)
      ylim <- c(0, max_lim)
    }

    #default parameters
    def_args <- list(pch=16, xlim=xlim, ylim=ylim,
                     xlab = expression("Expected " ~ chi^2),
                     ylab = expression("Observed " ~ chi^2),...)


    ## Next, get a list of ... arguments
    dotargs <- list(...)

    ## And call the plot function passing NA, your ... arguments, and the default
    ## arguments that were not defined in the ... arguments.
    tryCatch(do.call("plot", c(list(x=p.e, y=p.o), def_args[!names(def_args) %in% names(dotargs)], dotargs)), warn=stop)

    abline(a = 0, b = 1)
    abline(a = 0, b = out$estimate, col = "red")
    rp = vector('expression',1)
    if (out$estimate >= 1.0) {
      dig<- 5
    } else {
      dig<- 4
    }
    rp[1] = substitute(expression(lambda == MYVALUE),
                       list(MYVALUE = format(out$estimate,dig=dig)))[2]
    legend('top', legend = rp, bty = 'n')
    #par(mar = oldmargins)
  }
  out
}

#' QQ plot of observed P-values vs expected P-values.
#'
#' @author Slave Petrovski and Quanli Wang
#' @param P.perm Expected P-values from NULL distribution, which is generated through permutation in our example.
#' @param P.observed Observed P-values from data.
#' @param adjust.xy Indicates if the x-axis and y-axis should be adjusted based on their own range in the plot.
#' @param ... arguments passed the plot function.
#'
#' @return None.
#'
#' @examples
#' library(QQperm)
#' data("example.Ps")
#'
#' #print output to pdf file only if not running in interactive mode
#' if (!interactive()) {
#'  pdf("qqplot.pdf")
#' }
#'
#' qqplot(example.Ps$perm, example.Ps$observed)
#'
#' if (!interactive()) {
#'  dev.off()
#' }
qqplot <- function(P.perm, P.observed,adjust.xy = TRUE,...) {
  #do QQ plot
  e <- -log10(P.perm)
  o <- -log10(P.observed)

  if (adjust.xy) {
    xlim <- c(0,round(max(e) + 1))
    ylim <- c(0,round(max(o) + 1))
  } else {
    max_lim <- max(max(e), max(o)) + 1
    xlim <- c(0, max_lim)
    ylim <- c(0, max_lim)
  }

  #default parameters
  def_args <- list(pch=16, xlim=xlim, ylim=ylim,
                   xlab = expression("Expected " ~ chi^2),
                   ylab = expression("Observed " ~ chi^2),...)

  ## Next, get a list of ... arguments
  dotargs <- list(...)
  ## And call the plot function passing NA, your ... arguments, and the default
  ## arguments that were not defined in the ... arguments.
  tryCatch(do.call("plot", c(list(x=e, y=o), def_args[!names(def_args) %in% names(dotargs)], dotargs)), warn=stop)
  # Add diagonal
  abline(0,1,col="red")
}
