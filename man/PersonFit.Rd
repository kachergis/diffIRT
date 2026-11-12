\name{PersonFit}
\alias{PersonFit}
\alias{print.PersonFit}

\title{Person fit for diffIRT models}
\description{
  Flags persons whose responses and response times are unlikely under a fitted \code{diffIRT} model.}

\usage{
PersonFit(object, nsim=1000, seed=NULL, control=list())
\method{print}{PersonFit}(x, alpha=.05, ...)
}
\arguments{
  \item{object}{ A \code{diffIRT} object.}
  \item{nsim}{ Number of persons simulated from the fitted model for the reference distribution.}
  \item{seed}{ Optional seed for the simulation.}
  \item{control}{ A list of control values: \code{nq} (Gauss-Hermite nodes per latent variable for the
    log-likelihood, default 11) and \code{threads} (see \code{\link{diffIRT}}).}
  \item{x}{ A \code{PersonFit} object.}
  \item{alpha}{ Significance level used in the printed summary.}
  \item{...}{ Not used.}
}

\details{
  The statistic is each person's marginal log-likelihood \eqn{l_i = \log p(x_i, rt_i)} (responses and response times,
  integrated over the person distribution). Its distribution under the model is obtained by Monte Carlo: \code{nsim}
  persons are simulated from the fitted model, and for each pattern of answered items their log-likelihood is
  computed over those items. The p-value is the proportion of simulated persons with a log-likelihood at most
  \eqn{l_i}; \code{lz} standardizes \eqn{l_i} by the mean and standard deviation of the reference distribution.
  Because the observed and simulated statistics are the same function of the data, the p-values are calibrated
  (apart from the estimation of the item parameters, which matters little for large samples).

  The statistic detects persons whose responses and response times are inconsistent across items, e.g.
  inattentive responding on part of the test (slow responses at chance accuracy) or random responding. Fast
  responses at chance accuracy on \emph{all} items are not misfit: they are what the model predicts for a person with
  a small boundary separation, and show up as a low \code{gamma[p]} in \code{\link{factest}}. Rapid guessing on
  part of the test is detected when the item parameters are known, but item parameters estimated with the guessers
  included are biased towards them (the guesses pull the non-decision times down), which masks them. In
  simulations (\code{sim/fit_statistics.R}; N = 500, 10 items, 5\% aberrant persons on half the items), 4-6\% of
  regular persons were flagged at \eqn{\alpha = .05}; see \code{NEWS.md} for the power against each kind of aberrance.
}
\value{
  An object of class \code{PersonFit}, a list with
  \item{persons}{ A data frame with, per person, the log-likelihood \code{ll}, the standardized \code{lz}, the
    p-value \code{p} and the number of answered items \code{items}.}
  \item{nsim, nq}{ As used.}
}
\seealso{
\code{\link{RTFit}} and \code{\link{RespFit}} for item fit, \code{\link{factest}} for person parameters.
}
\examples{
\dontrun{
data=simdiff(N=500,nit=10,model="D")
res=diffIRT(data$rt,data$x,model="D")
pf=PersonFit(res,seed=1)
pf
head(pf$persons[order(pf$persons$p),])
}
}
\keyword{models}
