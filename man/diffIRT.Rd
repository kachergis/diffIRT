\name{diffIRT}
\alias{diffIRT}

\title{
  Fitting diffusion IRT models}
\description{
  This function fits the D-diffusion or the Q-diffusion IRT model to response and response time data using 
  marginal maximum likelihood. Item parameters that are estimated are item drift rate, v[i] item boundary separation, a[i], 
  and item specific non-decision time, Ter[i]. Population parameters that are estimated are omega[gamma] and omega[theta], which are 
  scale parameters for the person boundary separation and person drift rate respectively. Parameters can be submitted 
  to equality and/or fixed value constraints.}

\usage{
diffIRT(rt, x, model="D", constrain=NULL,
        start=NULL, se=F, control=list())
}
\arguments{
  \item{rt}{ A matrix of size \code{N} by \code{nit} containing the response times, where \code{N} is the number of subjects and 
    \code{nit} is the number of items. NA's are allowed. If a given element in rt is NA, the corresponding element in x is also treated as missing 
    and vice versa.}
  \item{x}{ A matrix of \code{N} by \code{nit} containing the responses. A's are allowed. If a given element in \code{rt} is NA, the corresponding 
    element in \code{x} is also treated as missing and vice versa.}
  \item{model}{ String; Either "D" to fit the D-diffusion IRT model or "Q" to fit the Q-diffusion IRT model.}
  \item{constrain}{ If \code{NULL}, the unconstrained model is fitted. Otherwise, one could use a manual setup or a preprogrammed setup. In the manual setup 
    a vector of length \code{3*nit+2} should be provided. Each element of this vector corresponds to a parameter in the model. Each parameter should be consecutively numbered 
    from '1' onwards. Equality constraints can be imposed by giving two or more parameters the same number. Fixed parameter constraints 
    can be imposed by putting \code{0}'s at the corresponding elements. The value to which those parameters should be fixed need then to be 
    be specified in the \code{start} argument. The first 1 to \code{nit} elements in \code{constrain} correspond to the item boundary parameters 
    (\code{ai}), the next \code{nit} elements correspond to the item drift parameters (\code{vi[i]}), the next \code{nit} elements 
    correspond to the item non-decision time parameters (\code{ter[i]}), and the final 2 elements contain the population scale parameters of the
    person boundary distribution (\code{omega[gamma]}) and the person drift distribution (\code{omega[theta]}). Note that \code{omega[gamma]} and \code{omega[theta]} 
    can be fixed to 0 to fit a model without random drift rates and random person boundaries respectively, see \bold{Examples}. In the preprogrammed setups argument \code{start} 
    should be either equal to \code{"ai.equal", "vi.equal", "ter.equal",} or \code{"scale.equal"} to equate respectively \code{a[i], v[i], or Ter[i]} 
    across items or to equate \code{omega[gamma]} and \code{omega[kappa]}. See \bold{Details} for more explanation on the parameters in the different models and see \bold{Examples} for examples on constraining parameters. }
  \item{start}{ If \code{NULL} the starting values are automatically chosen, see \bold{Details}. Otherwise, a vector of length \code{3 x nit + 2} should be 
    provided. Each element of this vector corresponds to a parameter from the model, similarly as in the \code{constrain} argument. For each 
    parameter a starting value can be provided. If \code{NA}, the starting value for that parameter is chosen automatically. If a parameter 
    is submitted to a fixed value constraint in the \code{constrain} argument, the value of that parameter should be provided (i.e., an \code{NA} 
    is not allowed). See \bold{Examples}.}
  \item{se}{ Logical; Denoting whether standard errors of the parameters should be estimated (this will increase computation time, default is \code{F}).}
  \item{control}{ a list of control values for the optimizations} 
  \describe{
        \item{nq}{ Vector of any length containing the number of quadrature points that should be used in successive optimizations. Optimisation will start 
          using the number of quadratures points specified in the first element. When converged, the resulting parameter estimates will be used as 
          starting values in a new optimisation using the number of quadrature points from the second element. This process continues until the
          end of vector \code{nq} is reached. Default is c(7).}
        \item{method}{ The optimisation method used by \code{\link{optim}}. Default "BFGS".}
        \item{eps}{ Precision with which the infinite integral in the likelihood function is approximated. Default is 0.01. See \bold{Details}.}
        \item{delta}{ Precision used in the finite difference approximation of the gradient. Default is 1e-7.}
        \item{trace}{ See \code{\link{optim}} for details and default.}
        \item{fnscale}{ See \code{\link{optim}} for details and default.}
        \item{parscale}{ See \code{\link{optim}} for details and default.}
        \item{maxit}{ See \code{\link{optim}} for details. Default is 1999.}
        \item{reltol}{ See \code{\link{optim}} for details and default.}
    }
  }

\details{
  \code{diffIRT} fits either the D-diffusion or the Q-diffusion IRT model to data by 
    minimizing -2 times the log marginal likelihood function using \code{\link{optim}}. 
    In the diffusion IRT model the traditional parameters from the diffusion model, 
    boundary separation and drift rate are decomposed into person and item parameters 
    (see van der Maas et al., 2011;  Tuerlinckx & De Boeck, 2005). 
    This results in: item boundary parameter \eqn{a[i]}, item drift parameter \eqn{v[i]}, 
    person boundary \eqn{gamma[p]}, and person drift \eqn{theta[p]}. The model 
    for the responses in the D-diffusion IRT model is then: \deqn{logit(x[p,i]) = gamma[p]/a[i]  x (theta[i]-v[i]),} 
    where \eqn{gamma[p] >= 0} and \eqn{ai[i] >= 0}. The Q-diffusion IRT model for the 
    responses is given by: \deqn{logit(x[p,i]) = gamma[p]/a[i]  x (theta[p]/v[i]).} 
    where \eqn{gamma[p] >= 0}, \eqn{a[i] >= 0}, \eqn{theta[p] >= 0}, and \eqn{theta[p] >= 0}. 
    As discussed in van der Maas et al. (2011), this setup makes the D-diffusion IRT model 
    suitable for personality data and the Q-diffusion IRT model for ability data. 
    The response times follow a distribution according to a Wiener process which includes 
    - in addition to the parameters above - a non-decision parameter for each 
    item, \eqn{Ter[i]}. For the random effects, \eqn{theta[p]} and \eqn{gamma[p]} scale parameters are estimated
    which are respectively, \eqn{omega[gamma]} and \eqn{omega[theta]}. The joint distribution of 
    responses and response times is evaluated using the approach 
    by Navarro & Fuss (2009). In this approach, the infinite sum in the density function
    of the diffusion model is being approximated with a maximum discrepancy of \code{epsilon}. This discrepancy can be
    set using the \code{eps} setting within the \code{control} argument. To facilitated numerical estimation, the natural logarithm 
    of the parameters that are strictly positive are estimated (i.e., \code{log(a[i])}, \code{log(Ter[i])}, 
    \code{log(omega[gamma])}, \code{log(omega[theta])}, and \code{log(v[i])} in the Q-diffusion model). However, in the output 
    the parameters are transformed back to their original scale. In addition, if \code{se=T}, 
    standard errors for the original parameters are calculated from the standard errors 
    of the transformed parameters using the delta method.
    
    Because gamma[p] and theta[p] in the Q-diffusion model and gamma[p] in the D-diffusion model can 
    only be positive, their population distribution is assumed to be normal on the log-scale. 
    As a consequence, gamma[p] and theta[p] follow a log-normal distribution. Thus, in the Q-diffusion model \code{gamma[p]} and \code{theta[p]} are
    log-normally distributed with scale parameters \code{omega[gamma]} and \code{omega[theta]} respectively. In case of the
    D-diffusion model, \code{theta[p]} is distributed log-normally with scale parameter \code{omega[gamma]}, and \code{theta[p]} 
    is distributed normally with scale parameter \code{omega[theta]} which equals the standard deviation.
    
    By default, starting values are calculated using the EZ-diffusion model (Wagemakers, van der Maas, & Grasman, 2007).
    
    For more details see Molenaar, Tuerlinkcx, & van der Maas (2013).
     
}
\value{
  An object of class \code{diffIRT} with values
   \item{ start.val}{ starting values}
   \item{ par}{ parameter estimates}
   \item{ par.log}{ log-transformed parameters (note that for the D-diffusion model v[i] is not transformed, see \bold{Details}.}
   \item{ std.err}{ if argument \code{se} equals \code{T}, it contains the standard errors 
      of the parameters}
   \item{ totLL }{ value of -2 times the log-marginal likelihood at convergence}
   \item{ npars}{ number of parameters in the model}
   \item{ AIC}{ AIC }
   \item{ BIC}{ BIC }
   \item{ sBIC}{ sample size adjusted BIC }
   \item{ DIC}{ DIC }
   \item{ subjLL}{ a vector of size \code{N} containing the individual contribution of each subject to the marginal likelihood of the model.}
   \item{ hist}{ a matrix containing the history of the optimisation procedure. First row contains the starting values. Next rows
        contain the successive optimizations using the different number of quadrature points as specified in \code{control}. Final
        4 columns contain \code{LL} (-2 times the log-marginal likelihood at convergence), \code{converg} (the convergence status returned by optim,
        \code{func} (the number of likelihood evaluations by \code{optim}), and \code{gradient} (the number of gradient evaluations by \code{optim}).}
   \item{ conv}{ integer convergence status returned by \link{optim}}
   \item{ nr_fail}{ the number of subjects for which the likelihood function was intractable (i.e., for these subjects, the log-likelihood approaches 
    minus infinity). This could be due to the starting values, or due to extremely fast response times.}
}
\references{ 

	Molenaar, D., Tuerlinkcx, F., & van der Maas, H.L.J. (2015). Fitting Diffusion Item Response Theory Models for Responses and Response Times Using the R Package diffIRT. 
	\emph{Journal of Statistical Software}, \bold{66(4)}, 1-34. URL \url{http://www.jstatsoft.org/v66/i04/}.

Navarro, D.J. & Fuss, I.G. (2009). Fast and accurate calculations for first-passagetimes in Wiener diffusion models.
\emph{Journal of Mathematical Psychology}, \bold{53}, 222-230.

Tuerlinckx, F., & De Boeck, P. (2005). Two interpretations of the discrimination parameter. \emph{Psychometrika}, \bold{70}, 629-650.

van der Maas, H.L.J., Molenaar, D., Maris, G., Kievit, R.A., & Borsboom, D. (2011). Cognitive	psychology meets psychometric theory: 
On the relation between process models for decision making and latent variable models for individual differences. 
\emph{Psychological Review}, \bold{118}, 339-356.

Wagenmakers, E. J., Van Der Maas, H. L., & Grasman, R. P. (2007). An EZ-diffusion model for response time and accuracy. 
\emph{Psychonomic Bulletin & Review}, \bold{14}, 3-22.

}
\author{ Dylan Molenaar \email{d.molenaar@uva.nl}}

\seealso{  
\code{\link{simdiff}} for simulating data according to the D-diffusion or Q-diffusion IRT model.
\code{\link{factest}} for estimation of factor scores (person drift rate and person boundary separation).
\code{\link{QQdiff}} and \link{RespFit} for model fit assessment.
\code{\link{summary.diffIRT}} for a overview of the model estimation results, including model fit statistics.
\code{\link{anova.diffIRT}} to conduct a likelihood ratio test between two nested diffIRT models.
\code{\link{coef.diffIRT}} to extract parameter estimates.
}
\examples{
\dontrun{
 # open extraversion data
data(extraversion)
x=extraversion[,1:10]
rt=extraversion[,11:20]
 
# fit an unconstrained D-diffusion model
res1=diffIRT(rt,x,model="D")          
         
# fit a model with equal item boundaries, a[i] using the manual setup
res2=diffIRT(rt,x,model="D",
  constrain=c(rep(1,10),2:11,12:21,22,23)) 

# fit a model with equal item boundaries, a[i] using the preprogrammed setup
res2=diffIRT(rt,x,model="D",
  constrain="ai.equal") 

# fit a model where all item drift parameters,vi, are fixed to 0. 
res3=diffIRT(rt,x,model="D",
  constrain=c(1:10,rep(0,10),11:20,21,22),
  start=c(rep(NA,10),rep(0,10),rep(NA,10),NA,NA))       

# fit a model without random person boundary parameters. 
res3=diffIRT(rt,x,model="D",
  constrain=c(1:30,0,31),
  start=c(rep(NA,30),0,NA))       

}
}  
\keyword{models}
\keyword{multivariate}
