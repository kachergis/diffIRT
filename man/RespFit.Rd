\name{RespFit}
\alias{RespFit}

\title{
  Assessing diffIRT model fit for the responses using limited information goodness-of-fit testing.}
\description{
  This function uses the procedure by Maydeu-Olivares & Joe (2005) to asses the goodness-of-fit of the responses 
  from a \code{diffIRT} object.}

\usage{
RespFit(object, order=2)
}
\arguments{
  \item{object}{ A \code{diffIRT} object for which the Maydeu-Olivares & Joe test needs to be conducted.}
  \item{order}{ Integer; The order of the moments to be compared see \bold{details}.}
  }

\details{
  \code{RespFit} is an implementation of the method outlined in Maydeu-Olivares & Joe (2005). The traditional Pearson 
  chi-square method are sub optimal in this case because in common IRT settings, contingency tables tend to be sparse. This 
  causes the asymptotic distribution of the traditional test statistic to depart from its theoretical distribution. In 
  the method proposed by Maydeu-Olivares & Joe, this problem is overcome by focussing on the first \emph{r} 
  moments (specified in \code{order}) of the observed and predicted response distributions. Choosing \code{order}
  to be equal to the number of items will result in the traditional chi-square test statistic. Commonly
  \code{order} is chosen to be small (e.g., 1 or 2).}

\note{
The degrees of freedom for the test statistic differ between the D-diffusion and Q-diffusion model as for the 
Q-diffusion model, \eqn{a[i]} and \eqn{v[i]} are not simultaneously identified in response data only. See Molenaar, Tuerlinckx, & van der Maas (2013) 
for more details.}

\section{Warning}{ 
For large numbers of items, this test becomes computationally infeasible.}

\value{
  Returns an object of class \code{RespFit} with entries:
  \item{Z}{ A matrix with predicted statistics, observed statistics, and Z-values}
  \item{Mr}{ The test statistic.}
  \item{df}{ Degrees of freedom.}
  \item{order}{ Order of the test statistic.}
  
}
\references{ 
Maydeu-Olivares, A., & Joe, H. (2005). Limited and full information estimation and testing in 2n contingency tables: 
A unified framework. \emph{Journal of the American Statistical Association}, \bold{100}, 1009-1020.

	Molenaar, D., Tuerlinkcx, F., & van der Maas, H.L.J. (2015). Fitting Diffusion Item Response Theory Models for Responses and Response Times Using the R Package diffIRT. 
	\emph{Journal of Statistical Software}, \bold{66(4)}, 1-34. URL \url{http://www.jstatsoft.org/v66/i04/}.
}

\author{ Dylan Molenaar \email{d.molenaar@uva.nl}}

\seealso{  
\code{\link{QQdiff}} for model fit assessment of the response times.
\code{\link{diffIRT}} for fitting diffusion IRT models.
\code{\link{factest}} for estimation of factor scores (person drift rate and person boundary separation).
\code{\link{simdiff}} for simulating data according to the D-diffusion or Q-diffusion IRT model.
}
\examples{
\dontrun{
 # open extraversion data
data(extraversion)
x=extraversion[,1:10]
rt=extraversion[,11:20]
 
# fit an unconstrained D-diffusion model
res1=diffIRT(rt,x,model="D")          

# Conduct the limited-information test
RespFit(res1, 2)
}  
}
\keyword{models}
\keyword{multivariate}

