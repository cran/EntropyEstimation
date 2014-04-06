\name{Entropy.sd}
\alias{Entropy.sd}
\title{Entropy.sd}
\description{
Returns the estimated asymptotic standard deviation for the Z estimator of Shannon's Entropy. Note that this is also the asymptotic standard deviation of the plug-in estimator. See Zhang and Grabchak (2014a) for details.
}
\usage{
Entropy.sd(x)
}
\arguments{
  \item{x}{Vector of counts. Must be integer valued. Each entry represents the number of observations of a distinct letter.}
}
\references{
Z. Zhang and M. Grabchak (2014a). Entropic representation and estimation of diversity indices. http://arxiv.org/abs/1403.3031.
}
\author{Lijuan Cao and Michael Grabchak}
\examples{
 x = c(1,3,7,4,8) # vector of counts
 Entropy.sd(x)  # Estimated standard deviation
 
 data = rbinom(10,20,.5)
 counts = tabulate(as.factor(data))
 Entropy.sd(counts)
}
