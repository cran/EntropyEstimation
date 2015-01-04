\name{MI.sd}
\alias{MI.sd}
\title{MI.sd}
\description{
Returns the estimated asymptotic standard deviation for the Z estimator of mutual information. Note that this is also the asymptotic standard deviation of the plug-in estimator. See Zhang and Zheng (2014) for details.}
\usage{
MI.sd(y)
}
\arguments{
  \item{y}{
Matrix of counts. Must be integer valued. Each entry represents the number of observations of a distinct combination of letters from the two alphabets.}
}
\references{
Z. Zhang and L. Zheng (2014). A Mutual Information Estimator with Exponentially Decaying Bias. 
}
\author{Lijuan Cao and Michael Grabchak}
\examples{
 x = matrix(c(0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
       0, 0, 0, 1, 0, 0, 1, 1, 0, 1,
       0, 0, 0, 2, 1, 0, 1, 0, 0, 1,
       0, 0, 0, 1, 1, 2, 0, 0, 0, 0,
       0, 0, 0, 3, 6, 2, 2, 0, 0, 0,
       2, 0, 2, 5, 6, 5, 1, 0, 0, 0,
       0, 0, 4, 6, 11, 5, 1, 1, 0, 1,
       0, 0, 5, 10, 21, 7, 5, 1, 0, 1,
       0, 0, 7, 11, 9, 6, 3, 0, 0, 1,
       0, 0, 4, 10, 6, 5, 1, 0, 0, 0),10,10,byrow=TRUE)
MI.sd(x)  

x = rbinom(100,20,.5)
y = rbinom(100,20,.5)
MI.sd(table(x,y))
}
