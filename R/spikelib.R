#' KN.test
#'
#' Estimate the number of spikes in a spiked covariance model using the procedure described in Kritchman and Nadler, (2009)
#'
#' @param S the sample covariance matrix
#' @param numObs  number of data observations
#' @param alpha significance level for the test
#'
#' @return
#' \enumerate{
#' \item numOfSpikes - The number of spikes K
#' \item sigma.sq - Estimated noise variance
#' }
#'
#' @details This method calculates the number of spikes in a spiked covariance model
#' by iteratively comparing a standardized eigenvalue to the quantile of Tracy-Wisdom
#' distribution for given significance level. The procedure is described in Kritchman and Nadler, (2009).
#'
#' @seealso \code{\link{sigmasq.est}}
#'
#' @examples
#' library(casp)
#' S<-diag(c(10,8,6,4,rep(1,6)))
#' KN.test.out<-KN.test(S,50)
#'
#' @references
#' Kritchman, S. and Nadler, B., 2009. Non-parametric detection of the number of signals: Hypothesis testing and random matrix theory.
#' IEEE Transactions on Signal Processing, 57(10), pp.3930-3941.
#'
#' @export

KN.test <- function(S, numObs, alpha = 0.05) {

  lambdas<- eigen(S)$values
  p <- length(lambdas); max.spikes <- min(numObs, p) - 1; sigma.sq <- 0
  s_alpha <- (-3/2 * log(4*sqrt(pi) * alpha/100 ))^(2/3)

  for(spikes in 1:max.spikes) {

    mu_np <- 1/numObs*(sqrt(numObs-1/2) + sqrt((p - spikes)-1/2))^2
    sigma_np <- sqrt(mu_np/numObs) * (1/sqrt(numObs-1/2) + 1 / sqrt((p-spikes)-1/2) )^(1/3)

    sigma.sq.temp <- sigmasq.est(lambdas,numObs,spikes)
    if(lambdas[spikes] <= sigma.sq.temp*(mu_np + s_alpha * sigma_np)) break
    sigma.sq <- sigma.sq.temp
  }

  est.sigma.sq <- if(spikes > 1) sigma.sq else  sum(lambdas)/p
  return(list(numOfSpikes = spikes - 1, sigma.sq = est.sigma.sq))
}

#' sigmasq.est
#'
#' Estimates the unknown noise variance based on the procedure described in Kritchman and Nadler, (2009)
#'
#' @param lambdas eigenvalues to which the distribution is fitted
#' @param numObs  number of data observations
#' @param numOfSpikes number of spikes in the spike covariance model
#'
#' @return
#' \enumerate{
#' \item sigma.sq - Estimated noise variance
#' }
#'
#' @details This method calculates an unbiased estimate of the noise variance
#' based on the procedure described in Kritchman and Nadler, (2009). It involves
#' solving a couple of nonlinear equations. This function is called by \code{\link{KN.test}}.
#'
#' @seealso \code{\link{KN.test}}
#'
#' @examples
#' library(casp)
#' lambdas<-c(10,8,6,4,rep(1,6))
#' sigmasq.est(lambdas,100,4)
#'
#' @references
#' Kritchman, S. and Nadler, B., 2009. Non-parametric detection of the number of signals: Hypothesis testing and random matrix theory.
#' IEEE Transactions on Signal Processing, 57(10), pp.3930-3941.
#'
#' @export

sigmasq.est <- function(lambdas , numObs, numOfSpikes) {

  p <- length(lambdas)
  sigma.sq <- 1/(p-numOfSpikes) * sum(lambdas[(numOfSpikes+1):p])
  counter = 0

  while(TRUE) {

    tmp.b <- lambdas[1:numOfSpikes] + sigma.sq - sigma.sq*(p-numOfSpikes)/numObs
    discriminant <- tmp.b^2 - 4 * lambdas[1:numOfSpikes]*sigma.sq
    if (any( discriminant < 0)) break

    rho <- 0.5*(tmp.b + sqrt(discriminant))
    sigma.sq.new <- 1/(p-numOfSpikes) * (sum(lambdas) - sum(rho))
    if (abs(sigma.sq.new - sigma.sq)/sigma.sq < 1e-8 || counter>500) break
    sigma.sq <- max(sigma.sq.new,1e-6)
    counter = counter+1
  }

  return(sigma.sq)
}
