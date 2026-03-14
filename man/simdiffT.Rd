\name{simdiffT}
\alias{simdiffT}

\title{
  Simulate data according to the traditional diffusion model.}
\description{
  This function simulates responses and response time data according to the traditional diffusion model for a single subject on
  a given number of trails. The parameters of the traditional diffusion model include: boundary separation, mean drift rate, 
  standard deviation of drift rate, variance of the process, and ter.}

\usage{
simdiffT(N,a,mv,sv,ter,vp,max.iter=19999,eps=1e-15)
}
\arguments{
  \item{N}{ number of trails.}
   \item{a}{ boundary separation.}
  \item{mv}{ mean of the normally distributed drift rates across trails.}
  \item{sv}{ standard deviation of the normally distributed drift rate across trails.}
  \item{ter}{ non-decision time.}
  \item{vp}{ variance of the process, which is a scaling parameter. Default equals 1.}
  \item{max.iter}{ Maximum number of iterations for the rejection algorithm. See \bold{Details}.}
  \item{eps}{ Convergence criterion for the rejection algorithm. See \bold{Details}}
  }

\details{
  Function \code{simdiffT} is an application of the rejection algorithm outlined in Tuerlinckx et al. (2001) subject 
  to normally distributed inter-trail variability in drift. In this algorithm, a proposal response time is sampled 
  from an exponential distribution. This proposal is accepted as actual response time when a specific condition is 
  satisfied (see Eq. 16 in Tuerlinckx, 2001). As this condition requires the approximation of an infinite sum, 
  a convergence criterion needs to be specified (see the argument \code{eps}). When the condition is not satisfied, 
  a new proposal response time is sampled. This is repeated until the proposal response time is accepted or when 
  \code{max.iter} has been reached.
  }

\value{
  Returns a list with the following entries:
  \item{rt}{ the simulated matrix of response times}
  \item{x}{ the simulated matrix of responses}
}
\references{ 
	Molenaar, D., Tuerlinkcx, F., & van der Maas, H.L.J. (2015). Fitting Diffusion Item Response Theory Models for Responses and Response Times Using the R Package diffIRT. 
	\emph{Journal of Statistical Software}, \bold{66(4)}, 1-34. URL \url{http://www.jstatsoft.org/v66/i04/}.

Tuerlinckx, F., Maris, E., Ratcliff, R., & De Boeck, P. (2001). A comparison of four methods 
for simulating the diffusion process. \emph{Behavior Research Methods, Instruments & Computers}, \bold{33}, 443-456.
}
\author{ Dylan Molenaar \email{d.molenaar@uva.nl}}

\seealso{  
\code{\link{diffIRT}} for fitting diffusion IRT models.
}
\examples{
\dontrun{

# simulate data accroding to the traditional diffusion model
set.seed(1310)
a=2
v=1
ter=2
sdv=.3
N=10000

data=simdiffT(N,a,v,sdv,ter)
rt=data$rt
x=data$x

# fit the traditional diffusion model (i.e., a restricted D-diffusion model, 
# see application 3 of the paper by Molenaar et al., 2013) 

diffIRT(rt,x,model="D",constrain=c(1,2,3,0,4),start=c(rep(NA,3),0,NA)) 

# this constrained model is a traditional diffusion model
# please note that the estimated a[i] value = 1/a
# and that the estimated v[i] value = -v                  

}  
}
\keyword{models}
\keyword{multivariate}


