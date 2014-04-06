\name{RenyiEq.z}
\alias{RenyiEq.z}
\title{RenyiEq.z}
\description{
Returns the Z estimator of Renyi's equivalent entropy. This estimator has exponentially decaying bias. When r=1 returns 1. See Zhang and Grabchak (2014a) for details.
}
\usage{
RenyiEq.z(x, r)
}
\arguments{
  \item{x}{
Vector of counts. Must be integer valued. Each entry represents the number of observations of a distinct letter.}
  \item{r}{
Order of Renyi's equivalent entropy. Must be a strictly positive real number.}
}
\references{
Z. Zhang and M. Grabchak (2014a). Entropic representation and estimation of diversity indices. http://arxiv.org/abs/1403.3031.}
\author{Lijuan Cao and Michael Grabchak}
\examples{
 x = c(1,3,7,4,8) 
 RenyiEq.z(x,2)  
 
 data = rbinom(10,20,.5)
 counts = tabulate(as.factor(data))
 RenyiEq.z(counts,2)
}
