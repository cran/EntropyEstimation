\name{Tsallis.z}
\alias{Tsallis.z}
\title{Tsallis.z}
\description{
Returns the Z estimator of Tsallis entropy. This is based on scaling and shifting the Z estimator of Renyi's equivalent entropy. When r=1 returns the Z estimator of Shannon's entropy. See Zhang and Grabchak (2014a) for details.}
\usage{
Tsallis.z(x, r)
}
\arguments{
  \item{x}{
Vector of counts. Must be integer valued. Each entry represents the number of observations of a distinct letter.}
  \item{r}{Order or Renyi's equivalent entropy this index is based on. Must be a strictly positive real number.}
}
\references{
Z. Zhang and M. Grabchak (2014a). Entropic representation and estimation of diversity indices. http://arxiv.org/abs/1403.3031.
}
\author{Lijuan Cao and Michael Grabchak}
\examples{
 x = c(1,3,7,4,8)
 Tsallis.z(x,2)  
 
 data = rbinom(10,20,.5)
 counts = tabulate(as.factor(data))
 Tsallis.z(counts,2)
}