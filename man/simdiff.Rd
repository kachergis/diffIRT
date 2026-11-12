\name{simdiff}
\alias{simdiff}

\title{
  Simulate data according to the D-diffusion or Q-diffusion IRT model.}
\description{
  This function simulates responses and response time data according to the D-diffusion or Q-diffusion IRT model. 
  }

\usage{
simdiff(N,nit,ai=NULL,vi=NULL,gamma=NULL,theta=NULL,ter=NULL,
        model="D",max.iter=19999,eps=1e-15,omega=NULL,rho=0,miss=0)
}
\arguments{
  \item{N}{ number of subjects.}
  \item{nit}{ number of items.}
  \item{ai}{ a vector of length \code{nit} containing the true values for the item boundary separation, a[i].}
  \item{vi}{ a vector of length \code{nit} containing the true values for the item drift rate, v[i].}
  \item{gamma}{ a vector of length \code{N} containing the true values for the person boundary separation, gamma[p].}
  \item{theta}{ a vector of length \code{N} containing the true values for the person drift rate, theta[p].}
  \item{ter}{ a vector of length \code{nit} containing the true values for the item non-decision time, ter[i].}
  \item{model}{ string; Either "D" to fit the D-diffusion IRT model or "Q" to fit the Q-diffusion IRT model.}
  \item{max.iter}{ maximum number of iterations for the rejection algorithm. See \bold{Details}.}
  \item{eps}{ convergence criterion for the rejection algorithm. See \bold{Details}}
  \item{omega}{ when the true parameters are not given: the standard deviations \code{c(omega[gamma], omega[theta])}
        of log \code{gamma[p]} and of \code{theta[p]} (D model) or log \code{theta[p]} (Q model). Default
        \code{c(.3, .5)} for the D model and \code{c(.3, .3)} for the Q model.}
  \item{rho}{ when the true parameters are not given: the correlation between log \code{gamma[p]} and
        \code{theta[p]} (D) or log \code{theta[p]} (Q). Default 0. Note that \code{\link{diffIRT}} assumes the two are
        independent.}
  \item{miss}{ proportion of responses (and their response times) set to \code{NA}, completely at random. Default 0.}
  }

\details{
  Function \code{simdiff} is an extension of the rejection algorithm outlined in Tuerlinckx et al. (2001). In this algorithm,
  a proposal response time is sampled from an exponential distribution. This proposal is accepted as actual response 
  time when a specific condition is satisfied (see Eq. 16 in Tuerlinckx, 2001). As this condition requires the 
  approximation of an infinite sum, a convergence criterion needs to be specified (see the argument \code{eps}). When
  the condition is not satisfied, a new proposal response time is sampled. This is repeated until the proposal response
  time is accepted or when \code{max.iter} has been reached.
  }

\value{
  Returns a list with the following entries:
  \item{rt}{ the simulated matrix of response times}
  \item{x}{ the simulated matrix of responses}
  \item{ai}{ true values for \code{ai[i]}}
  \item{vi}{ true values for \code{vi[i]}}
  \item{gamma}{ true values for \code{gamma[p]}}
  \item{theta}{ true values for \code{theta[p]}}
  \item{ter}{ true values for \code{ter[i]}}  
}
\references{ 
Tuerlinckx, F., Maris, E., Ratcliff, R., & De Boeck, P. (2001). A comparison of four methods 
for simulating the diffusion process. \emph{Behavior Research Methods, Instruments & Computers}, \bold{33}, 443-456.
}
\author{ Dylan Molenaar \email{d.molenaar@uva.nl}}

\seealso{  
\code{\link{diffIRT}} for fitting diffusion IRT models.
}
\examples{
\dontrun{
# simulate data accroding to D-diffusion model
data=simdiff(N=100,nit=10,model="D")                   

}  
}
\keyword{models}
\keyword{multivariate}


