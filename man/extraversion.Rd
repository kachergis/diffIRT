\name{extraversion}
\alias{extraversion}
\docType{data}
\title{Responses and Response Times to 10 Extraversion Items}
\description{
These data comprise responses (first 10 columns labelled 'X[]') and response times (next 10 columns
labelled 'T[]') of 146 subjects to 10 extraversion items with a binary answer scale. Each item consists
of a particular habit, e.g., 'active' and 'noisy'. Subjects were asked to indicate whether (yes/no)
these habits are applicable to their personalities. Response times are in seconds.
}
\format{
The specific habits for the 10 items are (translated from Dutch):
\describe{
    \item{item 1}{'active'}
    \item{item 2}{'noisy'}
    \item{item 3}{'energetic'}
    \item{item 4}{'enthusiastic'}
    \item{item 5}{'impulsive'}
    \item{item 6}{'jovial'}
    \item{item 7}{'viable'}
    \item{item 8}{'eupeptic'}
    \item{item 9}{'communicative'}
    \item{item 10}{'spontaneous'}
 }
}
\references{
	Molenaar, D., Tuerlinkcx, F., & van der Maas, H.L.J. (2015). Fitting Diffusion Item Response Theory Models for Responses and Response Times Using the R Package diffIRT. 
	\emph{Journal of Statistical Software}, \bold{66(4)}, 1-34. URL \url{http://www.jstatsoft.org/v66/i04/}.

}
\examples{

data(extraversion)
x=extraversion[,1:10]             # responses, 0 for 'no', 1 for 'yes'
rt=extraversion[,11:20]           # response times in seconds

}
\keyword{datasets}
