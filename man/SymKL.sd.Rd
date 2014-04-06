\name{SymKL.sd}
\alias{SymKL.sd}
\title{SymKL.sd}
\description{
Returns the estimated asymptotic standard deviation for the Z estimator of Symmetrized Kullback-Leibler's divergence. Note that this is also the asymptotic standard deviation of the plug-in estimator. See Zhang and Grabchak (2014b) for details.}
\usage{
SymKL.sd(x, y)
}
\arguments{
  \item{x}{
Vector of counts from first distribution. Must be integer valued. Each entry represents the number of observations of a distinct letter.}
  \item{y}{
Vector of counts from second distribution. Must be integer valued. Each entry represents the number of observations of a distinct letter.}
}
\references{
Z. Zhang and M. Grabchak (2014b). Nonparametric Estimation of Kullback-Leibler Divergence.
}
\author{Lijuan Cao and Michael Grabchak}
\examples{
 x = c(1,3,7,4,8) # first vector of counts
 y = c(2,5,1,3,6) # second vector of counts
 SymKL.sd(x,y)  # Estimated standard deviation
}