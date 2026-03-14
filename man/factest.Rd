\name{factest}
\alias{factest}

\title{
  Estimating factor scores for diffIRT models}
\description{
  This function estimates the person drift rate and person boundary separation for \code{diffIRT} objects.}

\usage{
factest(object,start=NULL,se=F, control=list())
}
\arguments{
  \item{object}{ A \code{diffIRT} object for which factor scores need to be estimated.}
  \item{start}{ If \code{NULL} starting values are automatically chosen. Otherwise, \code{start} 
        should be a vector of size 2 x \code{N}, where \code{N} denotes the number of subjects. The first \code{N} elements
        correspond to the starting values for person boundary separation (\code{ap}), the next \code{N} elements correspond 
        to the starting values for person drift rate (\code{vp}.
        \code{NA} are allowed.}
  \item{se}{ Logical; Denoting whether standard errors of the parameters should be estimated (can be time consuming, therefore, default is \code{F}).}
  \item{control}{ a list of control values for the optimisation} 
  \describe{
        \item{method}{ The optimisation method used by \code{\link{optim}}. Default "BFGS".}
        \item{eps}{ See \code{\link{optim}} for details and default.}
        \item{delta}{ See \code{\link{optim}} for details and default.}
        \item{trace}{ See \code{\link{optim}} for details and default.}
        \item{fnscale}{ See \code{\link{optim}} for details and default.}
        \item{parscale}{ See \code{\link{optim}} for details and default.}
        \item{maxit}{ See \code{\link{optim}} for details. Default is 1999. }
        \item{reltol}{ See \code{\link{optim}} for details and default.}
    }
  }

\details{
  \code{factest} returns empirical Bayes estimates of the person drift rate and 
    the person boundary separation. See \link{diffIRT} for more explanation concerning the parameters in the 
    D-diffusion and Q-diffusion IRT model.
}
\value{
  Function \code{factest} returns a matrix of parameter estimates and - if \code{se}=T - standard errors.
}
\references{ 
Navarro, D.J. & Fuss, I.G. (2009). Fast and accurate calculations for first-passagetimes in Wiener diffusion models.
\emph{Journal of mathematical psychology}, \bold{53}, 222-230.

Tuerlinckx, F., & De Boeck, P. (2005). Two interpretations of the discrimination parameter. \emph{Psychometrika}, \bold{70}, 629-650.

van der Maas, H.L.J., Molenaar, D., Maris, G., Kievit, R.A., & Borsboom, D. (2011). Cognitive	Psychology Meets Psychometric Theory: 
On the Relation Between Process Models for 	Decision Making and Latent Variable Models for Individual Differences. 
\emph{Psychological Review}, \bold{118}, 339-356.
}
\author{ Dylan Molenaar \email{d.molenaar@uva.nl}}

\seealso{  
\code{\link{diffIRT}} for fitting diffusion IRT models.
\code{\link{simdiff}} for simulating data according to the D-diffusion or Q-diffusion IRT model.
\code{\link{QQdiff}} and \link{RespFit} for model fit assesment.
}
\examples{
\dontrun{
 # simulate data accroding to D-diffusion model
data=simdiff(N=100,nit=10,model="D")                   

# fit an unconstrained model
res1=diffIRT(data$rt,data$x,model="D")          

# estimate factor scores
fs=factest(res1) 
}  
}
\keyword{models}
\keyword{multivariate}
