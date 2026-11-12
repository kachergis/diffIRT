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
The degrees of freedom are the number of moments minus the number of parameter combinations that the responses
identify, i.e. the rank of the matrix of derivatives of the moments with respect to the free parameters. Response
data alone do not identify every parameter: in the D-diffusion model the scale of the drift is free (2 x items + 1
identified combinations without constraints); in the Q-diffusion model accuracy depends only on \eqn{a[i] v[i]} and on
\eqn{omega[gamma]^2 + omega[theta]^2} (items + 1). See Molenaar, Tuerlinckx, & van der Maas (2013) for more details.

The moments, their covariance matrix and their derivatives are computed directly from joint probabilities of at most
2 x \code{order} items, so the test is feasible for long tests (about 1 s for 20 items and 11 s for 40 items with
\code{order=2}). diffIRT 1.5 and earlier enumerated all 2^\code{nit} response patterns (infeasible beyond about 14
items). Their Q-model derivatives were wrong (the item drifts were left on the log scale), and their df
subtracted 2 more than the number of estimated parameters, so their Mr values and p-values differ from the current ones.}

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

