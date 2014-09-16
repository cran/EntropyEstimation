\name{SymKL.z}
\alias{SymKL.z}
\title{SymKL.z}
\description{
Returns the Z estimator of Symetrized Kullback-Leibler Divergence, which has exponentialy decaying bias.  See Zhang and Grabchak (2014b) for details.}
\usage{
SymKL.z(x, y)
}
\arguments{
  \item{x}{
Vector of counts from first distribution. Must be integer valued. Each entry represents the number of observations of a distinct letter.}
  \item{y}{
Vector of counts from second distribution. Must be integer valued. Each entry represents the number of observations of a distinct letter.}
}
\references{
Z. Zhang and M. Grabchak (2014b). Nonparametric Estimation of Kullback-Leibler Divergence. Neural Computation, DOI 10.1162/NECO_a_00646.
}
\author{Lijuan Cao and Michael Grabchak}
\examples{
 x = c(1,3,7,4,8) 
 y = c(2,5,1,3,6) 
 SymKL.z(x,y) 
}