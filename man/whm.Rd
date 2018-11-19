\name{whm}
\alias{whm}

\title{
Weighted M-estimation
}
\description{
This function performs a weighted M-estimation described by Carroll and Ruppert (1982) with the Huber loss function. First, a M-estimation is performed on the data assuming that the variance is constant. The residuals of this model are used to robustly estimate the variance parameter. Then, a weighted M-estimation with variance as weight is used to update the regression parameters. These steps are iterated until desired convergence. 
 }
\usage{
whm(yy, xx, var.function = "power", tuning.para = 1.345, ite = 5)
}

\arguments{
  \item{yy}{Vector representing the response variable
}
  \item{xx}{Design matrix of the covariates including the intercept in the first column
}
\item{var.function}{Assumed function for the variance. "\code{power}" function corresponds to  \eqn{\sqrt(Var) = \sigma = \phi |x^T \beta|^{\gamma}} and "\code{exponential}" to \eqn{\sqrt(Var) = \sigma = \phi e^{\gamma |x^T \beta|}}.

}
\item{tuning.para}{Value of the tuning parameter associated with the loss function.  

}
\item{ite}{Number of iterations for the estimation procedure. 

}
}
%\details{
%efgdgd
%}
\value{
The function returns a list including

\item{esti}{ Value of the robust estimate}
\item{Std.Error}{ Standard error of the robust estimate}
\item{tunning}{ Optimum tunning parameter}
\item{R2}{ R-squared value}
}
\references{

Carroll, R. J., & Ruppert, D. (1982). Robust estimation in heteroscedastic linear models. The annals of statistics, 429-441.


}
\author{
Aurelien Callens, You-Gan Wang, Benoit Liquet.
}


\seealso{
\code{rlm} function from  package \code{MASS}
}
\examples{

library(MASS)
data(stackloss)

LS <- lm(stack.loss ~ stack.x)
RB <- rlm(stack.loss ~ stack.x, psi = psi.huber, k = 1.345)

yy <- stack.loss 
xx <- model.matrix(stack.loss ~ stack.x)

#With power function as variance function
WHM_p <- whm(yy, xx, var.function = "power", tuning.para = 1.345, ite = 5)

#With exponential function as variance function
WHM_e <- whm(yy, xx, var.function = "exponential", tuning.para = 1.345, ite = 5)


}
\keyword{regression}
