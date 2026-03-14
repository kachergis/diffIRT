\name{QQdiff}
\alias{QQdiff}

\title{
  Assessing diffIRT model fit for the response times using QQ-plots}
\description{
  This function plots the observed response times against the predicted response times for a \code{diffIRT} object.}

\usage{
QQdiff(object, items, plot=2, breaks=15, quant=NULL, maxRT=NULL)
}
\arguments{
  \item{object}{ A \code{diffIRT} object for which the QQ-plots need to be created.}
  \item{items}{ A vector denoting for which items the QQ-plots need to be created.}
  \item{plot}{ Integer; 1: only QQ-plots. 2: both a QQ-plot and a histogram containing the predicted and observed distribution.}
  \item{breaks}{ Number of breaks to be used in \link{hist} when \code{plot}=2.}
  \item{quant}{ The number of quantiles to be used. If \code{NULL}, the number of quantiles will equal the number of subjects.} 
  \item{maxRT}{ The maximum response time used in finding the quantiles of the theoretical distribution. If \code{NULL}, twice the maximum 
    observed response time is used for each item. Increasing \code{maxRT} will increase computation time and should 
    only be used when \link{uniroot} produces errors.} 
  }

\details{
  \code{QQdiff} calculates the predicted quantiles in the marginal response time distribution of the given 
    model (D-diffusion or Q-diffusion). 
}
\value{
  Function \code{QQdiff} returns a list with entries:
  \item{qexp}{ a vector with predicted quantiles.}
  \item{qobs}{ a vector with observed quantiles.}
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
\code{\link{factest}} for estimation of factor scores (person drift rate and person boundary separation).
\code{\link{simdiff}} for simulating data according to the D-diffusion or Q-diffusion IRT model.
}
\examples{
\dontrun{
# open rotation data
data(rotation)
x=rotation[,1:10]
rt=rotation[,11:20]
 
# fit an unconstrained Q-diffusion model
res1=diffIRT(rt,x,model="Q")           

# make QQ-plots and histograms for items 1 to 4.
QQdiff(res1, items=1:4, plot=2, maxRT=rep(50,4))
}  
}
\keyword{models}
\keyword{multivariate}
