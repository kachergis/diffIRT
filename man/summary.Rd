\name{summary}
\alias{summary.diffIRT}

\title{ Summary function for diffIRT objects}

\description{
  Summarizes the modeling results of a diffIRT object.
}

\usage{
\method{summary}{diffIRT}(object, digits=3, \dots)

}

\arguments{
  \item{object}{ a \code{diffIRT} object from which a summary of the modeling results is wanted.}
  \item{digits}{integer; number of decimal places the output is rounded to.}
  \item{\dots}{ additional parameters, currently not used.}
}

\details{
  Parameter estimates are displayed, including the standard errors if these have been estimated (see the \code{se}
  argument of \code{\link{diffIRT}}. In addition, -2 times the log-likelihood function, AIC, BIC, sample size adjusted
  BIC, and DIC, are provided.
}

\author{ Dylan Molenaar \email{d.molenaar@uva.nl}}

\seealso{
\code{\link{diffIRT}} for fitting diffusion IRT models.
}
\examples{
\dontrun{

# simulate data according to D-diffusion model
data=simdiff(N=100,nit=10,model="D")

# fit a D-diffusion model
res=diffIRT(data$rt,data$x,model="D")

# use the summary function to obtain modeling results
summary(res)
}
}
\keyword{models}
\keyword{multivariate}

