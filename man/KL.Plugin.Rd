\name{KL.Plugin}
\alias{KL.Plugin}
\title{KL.Plugin}
\description{
Returns the augmented plugin estimator of Kullback-Leibler Divergence.  See Zhang and Grabchak (2014b) for details.
}
\usage{
KL.Plugin(x, y)
}
\arguments{
  \item{x}{
Vector of counts from first distribution. Must be integer valued. Each entry represents the number of observations of a distinct letter.}
  \item{y}{
Vector of counts from second distribution. Must be integer valued. Each entry represents the number of observations of a distinct letter.}
}
\references{
Z. Zhang and M. Grabchak (2014b). Nonparametric Estimation of Kullback-Leibler Divergence. Neural Computation, 26(11): 2570-2593.
}
\author{Lijuan Cao and Michael Grabchak}
\examples{
 x = c(1,3,7,4,8) 
 y = c(2,5,1,3,6) 
 KL.Plugin(x,y)  
 KL.Plugin(y,x)  
}
