\name{predict.diffIRT}
\alias{predict.diffIRT}
\alias{fitted.diffIRT}
\alias{simulate.diffIRT}

\title{Predictions and simulated data from a fitted diffIRT model}
\description{
  \code{predict} gives the expected probability of a correct response and the expected response time per item,
  for the population or for given person parameters. \code{fitted} gives these expectations for each person
  in the data, at their factor scores. \code{simulate} draws new datasets from the fitted model.}

\usage{
\method{predict}{diffIRT}(object, gamma = NULL, theta = NULL, ...)
\method{fitted}{diffIRT}(object, ...)
\method{simulate}{diffIRT}(object, nsim = 1, seed = NULL, keep.missing = TRUE, ...)
}
\arguments{
  \item{object}{ A \code{diffIRT} object.}
  \item{gamma, theta}{ Person boundary separation and person drift rate (vectors of equal length). If both are
        \code{NULL}, \code{predict} returns population-level expectations, integrated over the estimated
        person distribution.}
  \item{nsim}{ Number of datasets to simulate.}
  \item{seed}{ Optional random seed (passed to \code{\link{set.seed}}).}
  \item{keep.missing}{ Logical; give the simulated data the same pattern of missing responses as the data
        in \code{object}. Default \code{TRUE}.}
  \item{...}{ For \code{fitted}: passed to \code{\link{factest}}. For \code{simulate}: passed to
        \code{\link{simdiff}} (e.g. \code{max.iter}).}
}
\details{
  For a person with boundary separation \eqn{gamma} and drift \eqn{theta} on item \eqn{i}, the boundary is
  \eqn{a = gamma/a[i]} and the drift is \eqn{v = theta - v[i]} (D-diffusion) or \eqn{v = theta/v[i]}
  (Q-diffusion). With an unbiased starting point, P(correct) = 1/(1 + exp(-a v)) and the expected response
  time is \eqn{Ter[i] + a/(2v) tanh(a v/2)} (\eqn{Ter[i] + a^2/4} for \eqn{v = 0}).

  \code{simulate} draws new person parameters from the estimated population distribution for each dataset,
  so it can be used for parametric bootstraps and posterior predictive checks.
}
\value{
  \code{predict} without person parameters: a data frame with columns \code{prob} and \code{rt}, one row per
  item. With person parameters, and \code{fitted}: a list with matrices \code{prob} and \code{rt}
  (persons x items).

  \code{simulate}: a list of \code{nsim} datasets, each as returned by \code{\link{simdiff}}, with attribute
  \code{"seed"}.
}
\seealso{\code{\link{diffIRT}}, \code{\link{factest}}, \code{\link{simdiff}}}
\examples{
\dontrun{
data=simdiff(N=500,nit=10,model="D")
res=diffIRT(data$rt,data$x,model="D")
predict(res)                                  # population-level expectations per item
p=fitted(res)                                 # per person, at the factor scores
sims=simulate(res,nsim=10,seed=1)             # parametric bootstrap datasets
}
}
\keyword{models}
