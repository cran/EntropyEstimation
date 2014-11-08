\name{GenSimp.sd}
\alias{GenSimp.sd}
\title{GenSimp.sd}
\description{
Returns the estimated asymptotic standard deviation of the Z estimator of the generalized Simpson's index of order r, i.e. of the index sum_k p_k(1-p_k)^r. This estimate of the standard deviation is based on the formula in Zhang and Grabchak (2014a) and not the one in Zhang and Zhou (2010).
}
\usage{
GenSimp.sd(x, r)
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
 GenSimp.sd(x,2)  
 
 data = rbinom(10,20,.5)
 counts = tabulate(as.factor(data))
 GenSimp.sd(counts,2)
}