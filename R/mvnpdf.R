#' Compute the density of a multivariate normal distribution
#'
#' mvnpdf() computes the value of the density of any multivariate normal distribution, at n points supplied by the user.
#'
#' @param x a matrix, with n columns (the observations) and  p rows
#' @param mean a vector of means
#' @param varcovM a variance-covariance matrix
#' @param Log logical, indicating whether the function should return the log of the density? Default is \code{TRUE}
#'
#' @return a list containing the matrix x, and a vector of length n of the multivariate normal distribution density values at those points
#' @export
#'
#' @importFrom mvtnorm dmvnorm
#'
#' @examples
#' matrix <- matrix(data = rnorm(100), nrow = 10, ncol = 10)
#' mean <- apply(matrix, 1, mean)
#' covmat <- diag(1, 10, 10)
#' mvnpdf(x = matrix, mean = mean, varcovM = covmat)
#'

mvnpdf <- function(x, mean = rep(0, nrow(x)), varcovM = diag(nrow(x)), Log = TRUE){

  n <- ncol(x)
  p <- nrow(x)
  x0 <- x - mean
  Rinv <- solve(varcovM)
  LogDetvarcovM <- log(det(varcovM))

  y <- NULL
  for (j in 1:n) {
    yj <- - p/2 * log(2*pi) - 0.5 * LogDetvarcovM -
      0.5 * t(x0[, j]) %*% Rinv %*% x0[, j]
    y <- c(y, yj)
  }

  if (!Log) {
    y <- exp(y)
  }

  res <- list(x = x, y = y)
  return(res)
}
