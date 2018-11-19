\name{rlmDD}
\alias{rlmDD}

\title{
Data driven robust methods
}
\description{
Robust estimation often relies on a dispersion function that is more slowly varying at large values than the squared function. However, the choice of tuning constant in dispersion
function may impact the estimation efficiency to a great extent. For a given family of dispersion functions, we suggest obtaining the `best' tuning constant from the data so that the asymptotic efficiency is maximized.

This library provides a robust linear regression with a tuning parameter being automatically chosen  to provide the necessary resistance against outliers. The robust (loss) functions include the Huber, Tukey bisquare and the exponential loss.
 }
\usage{
rlmDD(yy, xx, beta0, betaR, method, plot)
}

\arguments{
  \item{yy}{Vector representing the response variable
}
  \item{xx}{Design matrix of the covariates inclunding the intercept in the first column
}
  \item{beta0}{Initial parameter estimate using \code{lm}

}
  \item{betaR}{Robust estimate of beta with a fixed tuning constant using \code{rlm}

}
  \item{method}{Huber, Bisquare or Exponential

}
  \item{plot}{"Y" gives a plot: the efficiency factor versus a range of tunning parameter values.

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
Wang, Y-G., Lin, X., Zhu, M., & Bai, Z. (2007). Robust estimation using the Huber function with a data-dependent tuning constant. Journal of Computational and Graphical Statistics, 16(2), 468-481.

Wang, X., Jiang, Y., Huang, M., & Zhang, H. (2013).  Robust variable selection with exponential squared loss. Journal of the American Statistical Association, 108, 632-643.

Wang, N., Wang, Y-G., Hu, S., Hu, Z. H., Xu, J., Tang, H., & Jin, G. (2018). Robust Regression with Data-Dependent Regularization Parameters and Autoregressive Temporal Correlations. Environmental Modeling & Assessment, in press.


}
\author{
You-Gan Wang, Na Wang
}


\seealso{
\code{rlm} function from  package \code{MASS}
}
\examples{

library(MASS)
data(stackloss)

LS <- lm(stack.loss ~ stack.x)
RB <- rlm(stack.loss ~ stack.x, psi = psi.huber, k = 1.345)
DD1 <- rlmDD(stack.loss, stack.x, LS$coef, RB$coef, method = "Huber", plot = "Y")

LS <- lm(stack.loss ~ stack.x)
RB <- rlm(stack.loss ~ stack.x, psi = psi.bisquare, c = 4.685)
DD2 <- rlmDD(stack.loss, stack.x, LS$coef, RB$coef, method = "Bisquare", plot = "Y")

LS <- lm(stack.loss ~ stack.x)
RB <- rlm(stack.loss ~ stack.x, psi = psi.huber, k = 1.345)
DD3 <- rlmDD(stack.loss, stack.x, LS$coef, RB$coef, method = "Exponential", plot = "Y")


## Plasma dataset


data(plasma)

y <- plasma$y
x <- cbind(plasma$calories, plasma$dietary)

LS <- lm(y ~ x)
RB <- rlm(y ~ x, psi = psi.huber, k = 1.345)
DD.h <- rlmDD(y, x, LS$coef, RB$coef, method = "Huber", plot = "Y")

LS <- lm(y ~ x)
RB <- rlm(y ~ x, psi = psi.bisquare, c = 4.685)
DD.b <- rlmDD(y, x, LS$coef, RB$coef, method = "Bisquare", plot = "Y")

LS <- lm(y ~ x)
RB <- rlm(y ~ x, psi = psi.huber, k = 1.345)
DD.e <- rlmDD(y, x, LS$coef, RB$coef, method = "Exponential", plot = "Y")


}
\keyword{regression}
