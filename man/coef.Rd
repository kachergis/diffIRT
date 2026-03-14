\name{coef}
\alias{coef.diffIRT}

\title{
  Return estimated parameters}
\description{
  Returns estimated item and population parameters from a \code{diffIRT} object}

\usage{
\method{coef}{diffIRT}(object, \dots)
}

\arguments{
  \item{object}{ a \code{diffIRT} object from which parameter estimates need to be extracted.}
  \item{\dots}{ additional parameters, currently not used.}
  }

\value{
  Returns a list with two entries:
  \item{item}{ the estimated item parameters, item boundary (\code{ai}), item drift (\code{vi}), 
    and item non-desicion time (\code{ter}).}
  \item{pop}{ the estimated population parameters, omega_ap and omega_vp.}
}

\author{ Dylan Molenaar \email{d.molenaar@uva.nl}}

\seealso{  
\code{\link{diffIRT} for fitting diffusion IRT models.}
}
\examples{  
\dontrun{
# simulate data accroding to D-diffusion IRT model
data=simdiff(N=100,nit=10,model="D")          

# fit the D-diffusion IRT model
res1=diffIRT(data$rt,data$x)

# extract parameter estimates 
coef(res1)
}  
}
\keyword{models}
\keyword{multivariate}


