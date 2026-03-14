\name{brightness}
\alias{brightness}
\docType{data}
\title{A Simulated Response Time Dataset according to an Experimental Design}
\description{
The data are simulated according to a design similar as that of a real 
brightness discrimination experiment by Ratcliff & Rouder (1998). In this experiment, 
the subject had to decide for a number of trials whether the brightness of a 
stimulus (a randomly generated array of pixels displayed on a computer screen) 
was either 'high' or 'low'. The true 
brightness of the stimuli were manipulated into a number of levels and 
administered with a speed instruction ("respond as fast as possible") and 
with an accuracy instruction ("respond as accurate as possible"). 
Present dataset was simulated according to a design with 6 different brightness levels and 
2 speed instructions resulting in 12 conditions. In the brightness data matrix, 
the first 12 columns are the responses and the next 12 columns are the response times. 
Each trial is assigned to a separate row with the response time of that trial in 
the corresponding column and \code{NA}'s on the remaining columns. Similarly 
for the responses. In addition, the data are arranged in such a way that the 
first 6 conditions are the speed instructed stimuli and the next 6 conditions 
are the corresponding accuracy instructed versions of these stimuli. See below for 
an example how to analyse these data using the \bold{diffIRT} package (taken from 
Molenaar, Tuerlinckx, & van der Maas, 2015).}

\references{
	Molenaar, D., Tuerlinkcx, F., & van der Maas, H.L.J. (2015). Fitting Diffusion Item Response Theory Models for Responses and Response Times Using the R Package diffIRT. 
	\emph{Journal of Statistical Software}, \bold{66(4)}, 1-34. URL \url{http://www.jstatsoft.org/v66/i04/}.


        Ratcliff, R., & Rouder, J. N. (1998). Modeling response times for two-choice 
	decisions. \emph{Psychological Science}, \bold{9(5)}, 347-356.
}
\examples{
data(brightness)
x=brightness[,1:12]
rt=brightness[,13:24]

\dontrun{
res = diffIRT(rt,x,model="D",constrain=c(rep(1,6),
rep(2,6),3:8,3:8,rep(9,12),0,10), start=c(rep(NA,36),0,NA))
}
}
\keyword{datasets}
