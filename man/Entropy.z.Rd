\name{Entropy.z}
\alias{Entropy.z}
\title{Entropy.z}
\description{
Returns the Z estimator of Shannon's Entropy. This estimator has exponentially decaying bias. See Zhang (2012), Zhang (2013), and Zhang and Grabchak (2014a) for details.
}
\usage{
Entropy.z(x)
}
\arguments{
  \item{x}{Vector of counts. Must be integer valued. Each entry represents the number of observations of a distinct letter.}
}
\references{
Z. Zhang (2012). Entropy estimation in Turing's' perspective. Neural Computation 24(5), 1368--1389.

Z. Zhang (2013). Asymptotic normality of an entropy estimator with asymptotically decaying bias. IEEE Transactions on Information Theory 59(1), 504--508.

Z. Zhang and M. Grabchak (2014a). Entropic representation and estimation of diversity indices. http://arxiv.org/abs/1403.3031.
}
\author{Lijuan Cao and Michael Grabchak}
\examples{
 x = c(1,3,7,4,8) 
 Entropy.z(x)  
 
 data = rbinom(10,20,.5)
 counts = tabulate(as.factor(data))
 Entropy.z(counts)
}
