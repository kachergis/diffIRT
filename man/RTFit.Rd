\name{RTFit}
\alias{RTFit}
\alias{print.RTFit}

\title{Item fit of responses and response times for diffIRT models}
\description{
  Chi-square tests of the fit of the joint distribution of responses and response times, per item and over
  all items, for a \code{diffIRT} object.}

\usage{
RTFit(object, bins=10, control=list())
\method{print}{RTFit}(x, digits=3, ...)
}
\arguments{
  \item{object}{ A \code{diffIRT} object.}
  \item{bins}{ Number of response time bins per item. The bin edges are the quantiles of the item's predicted
    marginal response time distribution (deciles for \code{bins=10}).}
  \item{control}{ A list of control values: \code{nq} (Gauss-Hermite nodes per latent variable used to integrate the
    cell probabilities over the person distribution, default 21), \code{panels} (Gauss-Legendre panels per bin, default 12)
    and \code{h} (step for the numerical derivatives of the cell probabilities, default 1e-4).}
  \item{x}{ An \code{RTFit} object.}
  \item{digits}{ Number of digits for the p-values.}
  \item{...}{ Not used.}
}

\details{
  For each item, the persons who answered it are classified into \code{2 x bins} cells: the response (0 or 1)
  crossed with the response time bin. The observed cell proportions \eqn{p} are compared with the cell
  probabilities \eqn{\pi(\hat\theta)} implied by the fitted model, integrated over the person distribution.

  The test accounts for the estimation of the item and population parameters. Because \eqn{\hat\theta} is the
  maximum likelihood estimate, \eqn{Cov(p, \hat\theta) = \Delta V}, with \eqn{\Delta} the derivatives of the cell
  probabilities with respect to the free parameters and \eqn{V} the covariance matrix of the estimates
  (\code{\link{vcov}}). The residuals \eqn{r = p - \pi(\hat\theta)} therefore have covariance
  \eqn{\Sigma = \Xi - \Delta V \Delta'}, where \eqn{\Xi} is the multinomial covariance of the cell proportions
  (for two items, scaled by the number of persons who answered both). The statistic \eqn{r' \Sigma^{-} r} is
  asymptotically chi-square with degrees of freedom equal to the rank of \eqn{\Sigma}: \code{2 x bins - 1} per
  item (occasionally one less, when a direction is numerically singular), and the sum over items for the overall
  test. Bins defined by estimated quantiles give the same limiting distribution as bins fixed at the true
  quantiles (Moore & Spruill, 1975). Missing responses are assumed to be missing completely at random.

  If \code{object} was fitted without \code{se=TRUE}, the covariance matrix of the estimates is computed first.

  In simulations (\code{sim/fit_statistics.R} in the source repository; N = 500, 10 items, 10\% missing), the
  item and overall tests rejected a correct model at close to the nominal 5\% rate, and detected an item whose
  response times were delayed by 0.3 s for 30\% of the persons.
}
\value{
  An object of class \code{RTFit}, a list with
  \item{items}{ A data frame with, per item, the number of persons \code{n}, the statistic \code{stat}, its
    degrees of freedom \code{df} and p-value \code{p}.}
  \item{overall}{ The statistic, df and p-value of the test over all items.}
  \item{cells}{ A list with, per item, the observed and expected cell proportions.}
  \item{edges}{ A list with, per item, the response time bin edges.}
  \item{bins, model}{ As used.}
}
\references{
Maydeu-Olivares, A., & Joe, H. (2005). Limited- and full-information estimation and goodness-of-fit testing in
2^n contingency tables: A unified framework. \emph{Journal of the American Statistical Association}, \bold{100}, 1009-1020.

Moore, D.S., & Spruill, M.C. (1975). Unified large-sample theory of general chi-squared statistics for tests of fit.
\emph{The Annals of Statistics}, \bold{3}, 599-616.
}
\seealso{
\code{\link{RespFit}} for the fit of the responses, \code{\link{QQdiff}} for QQ-plots of the response times,
\code{\link{PersonFit}} for person fit.
}
\examples{
\dontrun{
data=simdiff(N=500,nit=10,model="D")
res=diffIRT(data$rt,data$x,model="D",se=TRUE)
RTFit(res)
}
}
\keyword{models}
