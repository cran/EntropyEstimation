\name{Tsallis.sd}
\alias{Tsallis.sd}
\title{Tsallis.sd}
\description{
Returns the estimated asymptotic standard deviation for the Z estimator of Tsallis Entropy. Note that this is also the asymptotic standard deviation of the plug-in estimator. See Zhang and Grabchak (2014a) for details.
}
\usage{
Tsallis.sd(x, r)
}
\arguments{
  \item{x}{
Vector of counts. Must be integer valued. Each entry represents the number of observations of a distinct letter.}
  \item{r}{ Order of Tsallis entropy. Must be a strictly positive real number. Not allowed to be 1, in that case use Entropy.sd instead.
  	}
}
\references{
Z. Zhang and M. Grabchak (2014a). Entropic representation and estimation of diversity indices. http://arxiv.org/abs/1403.3031.
}
\author{Lijuan Cao and Michael Grabchak}
\examples{
 x = c(1,3,7,4,8)
 Tsallis.sd(x,2)  
 
 data = rbinom(10,20,.5)
 counts = tabulate(as.factor(data))
 Tsallis.sd(counts,2)
}