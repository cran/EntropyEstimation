\name{GenSimp.z}
\alias{GenSimp.z}
\title{GenSimp.z}
\description{
Returns the Z estimator of the generalized Simpson's index of order r, i.e. of the index sum_k p_k(1-p_k)^r. See Zhang and Zhou (2010) and Zhang and Grabchak (2014a) for details.
}
\usage{
GenSimp.z(x,r)
}
\arguments{
  \item{x}{
Vector of counts. Must be integer valued. Each entry represents the number of observations of a distinct letter.}
 \item{r}{Positive integer representing the order of the generalized Simpson's index. If a noninteger value is given then the integer part is taken. Must be strictly less than sum(x).}
}
\references{
Z. Zhang and M. Grabchak (2014a). Entropic representation and estimation of diversity indices. http://arxiv.org/abs/1403.3031.

Z. Zhang and J. Zhou (2010). Re-parameterization of multinomial distributions and diversity indices. Journal of Statistical Planning and Inference 140(7), 1731-1738.
}
\author{Lijuan Cao and Michael Grabchak}
\examples{
 x = c(1,3,7,4,8) 
 GenSimp.z(x,2)  
 
 data = rbinom(10,20,.5)
 counts = tabulate(as.factor(data))
 GenSimp.z(counts,2)
}
