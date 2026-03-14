\name{rotation}
\alias{rotation}
\docType{data}
\title{Responses and Response Times to 10 Mental Rotation Items}
\description{
These data comprise responses (first 10 columns labelled 'X[]') and 
response times (next 10 columns labelled 'T[]') of 121 subjects to 10 
binary mental rotation items. The data are part of the paper by van der Maas et al (2011).
These data are taken from a larger database published in Kievit (2010; see also Borst, Kievit, 
Thompson, & Kosslyn, 2011). Each item consists of a graphical display of two 3-dimensional objects. 
 The second object was either a rotated version of the first object, 
 or a rotated version of a different object. Subjects were asked 
 whether the second object was the same as the first object (yes/no). 
 The degree of rotation of the second object was either 50, 100, or 150 degrees. 
 Answers are coded to be correct (1) or false (0). Response times were recorded in seconds.
}
\format{
The specific rotation angles of the different items are:
\describe{
    \item{item 1}{'150'}
    \item{item 2}{'50'}
    \item{item 3}{'100'}
    \item{item 4}{'150'}
    \item{item 5}{'50'}
    \item{item 6}{'100'}
    \item{item 7}{'150'}
    \item{item 8}{'50'}
    \item{item 9}{'150'}
    \item{item 10}{'100'}
 }
}
\references{
Borst, G., Kievit, R. A., Thompson, W. L., & Kosslyn, S. M. (2011). Mental rotation is not easily cognitively penetrable. 
\emph{Journal of Cognitive Psychology}, \bold{23}, 60-75. 

Kievit, R. A. (2010). \emph{Representational inertia: The influence of associative knowledge on 3D mental transformations.} 
Unpublished manuscript, Department of Psychology, University of Amsterdam, The Netherlands.

	Molenaar, D., Tuerlinkcx, F., & van der Maas, H.L.J. (2015). Fitting Diffusion Item Response Theory Models for Responses and Response Times Using the R Package diffIRT. 
	\emph{Journal of Statistical Software}, \bold{66(4)}, 1-34. URL \url{http://www.jstatsoft.org/v66/i04/}.


van der Maas, H.L.J., Molenaar, D., Maris, G., Kievit, R.A., & Borsboom, D. (2011). Cognitive	psychology meets psychometric theory: 
On the relation between process models for decision making and latent variable models for individual differences. 
\emph{Psychological Review}, \bold{118}, 339-356.    
}
\examples{

data(rotation)
x=rotation[,1:10]             # responses, 0 for 'false', 1 for 'correct'
rt=rotation[,11:20]           # response times in seconds

}
\keyword{datasets}
