\name{factest}
\alias{factest}

\title{
  Estimating factor scores for diffIRT models}
\description{
  This function estimates the person drift rate and person boundary separation for \code{diffIRT} objects.}

\usage{
factest(object, start=NULL, se=F, control=list(), type=c("MAP","EAP","PV"), npv=5)
}
\arguments{
  \item{object}{ A \code{diffIRT} object for which factor scores need to be estimated.}
  \item{start}{ If \code{NULL} starting values are automatically chosen. Otherwise, \code{start} 
        should be a vector of size 2 x \code{N}, where \code{N} denotes the number of subjects. The first \code{N} elements
        correspond to the starting values for person boundary separation (\code{ap}), the next \code{N} elements correspond 
        to the starting values for person drift rate (\code{vp}.
        \code{NA} are allowed.}
  \item{se}{ Logical; Denoting whether standard errors of the parameters should be estimated (can be time consuming, therefore, default is \code{F}).
        For \code{type="EAP"} these are the posterior standard deviations. Ignored for \code{type="PV"}.}
  \item{type}{ \code{"MAP"} (default): posterior modes, as in earlier versions. \code{"EAP"}: posterior means.
        \code{"PV"}: plausible values (random draws from each person's posterior). See \bold{Details}.}
  \item{npv}{ Number of sets of plausible values for \code{type="PV"}. Default 5.}
  \item{control}{ a list of control values for the estimation. Each person's scores are found
        separately by a one-dimensional Brent search (in parallel when compiled with OpenMP).}
  \describe{
        \item{eps}{ Precision of the series for the diffusion density. Default: the value used to fit \code{object}.}
        \item{tol}{ Convergence tolerance of the per-person search. Default 1e-8.}
        \item{halfwidth}{ Half-width of the search interval around each person's starting value.
              Default \code{max(3, 10*sd)}, with sd the estimated population standard deviation of the
              latent variable. A warning is given when estimates end at the edge of the interval.}
        \item{adaptive}{ Logical; if the model was fitted with adaptive quadrature (the default in \code{\link{diffIRT}}),
              integrate over the other latent variable on each person's adaptive grid and start each search at the
              person's posterior mean. Default \code{TRUE}; \code{FALSE} uses the non-adaptive grid.}
        \item{threads}{ Number of OpenMP threads for the per-person searches. Default: the value used to fit
              \code{object} (see \code{\link{diffIRT}}); \code{NULL} uses the OpenMP default.}
        \item{nq.eap}{ Number of adaptive Gauss-Hermite nodes per latent variable for \code{type="EAP"}. Default 21.}
        \item{pv.draws}{ Number of proposal draws per person for each plausible value (importance resampling).
              Default 200. A warning is given if the effective sample size falls below 10.}
        \item{pv.items}{ Logical; for \code{type="PV"}, draw the item and population parameters for each set of plausible
              values from their asymptotic normal distribution (\code{\link{vcov}}), so that the plausible values reflect
              the uncertainty in these estimates. Needs a fit with \code{se=TRUE}. Default: \code{TRUE} if the fit
              has standard errors, \code{FALSE} otherwise.}
    }
  }

\details{
  \code{factest} returns empirical Bayes estimates of the person drift rate and 
    the person boundary separation. See \link{diffIRT} for more explanation concerning the parameters in the 
    D-diffusion and Q-diffusion IRT model.

  The prior is a normal density with the estimated population variance, on log boundary
    separation and on person drift (D-diffusion) or log person drift (Q-diffusion).
    diffIRT 1.5 and earlier used \code{exp(-.5*P/sd)} in place of \code{exp(-.5*P^2/sd^2)}, so
    their factor scores differ from the current ones (most visibly for person drift).

  \code{type="MAP"} finds the mode of each person's posterior for one latent variable at a time, integrating over
    the other. \code{type="EAP"} returns posterior means (and, with \code{se=TRUE}, posterior standard deviations) of
    \code{gamma[p]} and \code{theta[p]} on the reported scale, computed on a 21 x 21 adaptive Gauss-Hermite grid
    centred on each person's posterior (accurate to about 1e-7). EAP scores are less sensitive than MAP scores to the
    skew of the posterior of the boundary separation, which is lognormal, and are better behaved for short tests.

  Point estimates (MAP or EAP) are shrunk towards the population mean, so their variance underestimates the
    population variance and their use as observed variables biases secondary analyses (correlations, regressions).
    Plausible values avoid this: they are random draws from each person's posterior, and analyses are run on each
    set and combined with Rubin's rules. Draws are made by importance resampling: \code{pv.draws} candidates from a
    bivariate t (5 df) distribution centred on the person's posterior, weighted by posterior over proposal density.
    With \code{pv.items=TRUE} each set first draws the item and population parameters from their asymptotic
    distribution. The plausible values reproduce the population variances: see \code{sim/pv_variance.R} in the
    source repository. Plausible values use the random number generator, so call \code{\link{set.seed}} first for
    reproducible results.
}
\value{
  For \code{type="MAP"} or \code{"EAP"}, a matrix with columns \code{gamma[p]} and \code{theta[p]}, and, if
  \code{se}=T, \code{gamma[p].se} and \code{theta[p].se} (standard errors from the curvature for MAP, posterior standard
  deviations for EAP). For \code{type="PV"}, a list of \code{npv} such matrices (without standard errors), one per set
  of plausible values, with attribute \code{"items"} saying whether item parameters were drawn.
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
eap=factest(res1, type="EAP", se=TRUE)

# plausible values: refit with standard errors so that they reflect item parameter uncertainty
res2=diffIRT(data$rt,data$x,model="D",se=TRUE)
set.seed(1)
pv=factest(res2, type="PV", npv=10)
# e.g. correlation of person drift with an external variable z, averaged over the sets:
# mean(sapply(pv, function(p) cor(p[,"theta[p]"], z)))
}  
}
\keyword{models}
\keyword{multivariate}
