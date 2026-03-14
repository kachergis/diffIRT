\name{anova}
\alias{anova.diffIRT}

\title{ Conduct a likelihood ratio test on two nested \code{diffIRT} models}

\description{
  Using this function, a likelihood ratio test is conducted on two nested \code{diffIRT} model. Results are 
  printed to the screen togheter with the AIC, BIC, sample size adjusted BIC, and DIC comparative fit indices for both models.
}

\usage{
\method{anova}{diffIRT}(object, object2, \dots)

}

\arguments{
  \item{object}{ a \code{diffIRT} object that is nested within the model provided in \code{object2}.}
  \item{object2}{ a \code{diffIRT} object.}
  \item{\dots}{ additional parameters, currently not used.}
}

\author{ Dylan Molenaar \email{d.molenaar@uva.nl}}

\seealso{
\code{\link{diffIRT} for fitting diffusion IRT models.}
}
\examples{
\dontrun{

# simulate data according to a D-diffusion model 
# with equal a[i] parameters 
data=simdiff(100,10,
  ai=rep(.3,10),vi=seq(-1,1,length=10),ter=runif(10,2,3),
  gamma=rlnorm(100,0,.3),theta=rnorm(100,0,.5),
  model="D")

# fit a full D-diffusion model
res=diffIRT(data$rt,data$x,model="D")

# fit a D-diffusion model subject to an 
# equality constraint across all a[i] parameters
res2=diffIRT(data$rt,data$x,model="D",constrain=c(rep(1,10),2:21,22,23))

# use the anova function to conduct the likelihood ratio test
anova(res2,res)

}
}
\keyword{models}
\keyword{multivariate}

