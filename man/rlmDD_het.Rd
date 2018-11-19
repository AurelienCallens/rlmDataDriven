\name{rlmDD_het}
\alias{rlmDD_het}

\title{
  Robust estimation for autoregressive models with heterogeneity
}
\description{
Performs robust regression for autoregressive models with heterogeneity. 

First, a M-estimation is performed on the data assuming that the variance is constant. The residuals of this model are used to robustly estimate the variance parameter. Then, a weighted M-estimation with variance as weight is used to update the regression parameters. These steps are repeated for different values of tuning parameter. The best tuning parameter is the one which minimizes the variance of the estimator. 

Finally, lagged term are built and added to the regression model therefore accounting for temporal correlations. The loss function used is Huber's function.
}
\usage{
rlmDD_het(yy, xx, var.function = c("power", "exponential"), phi.par = TRUE, 
tuning.para = NULL, step = 0.1, n.lag = NULL, print.summary = TRUE)
}

\arguments{
\item{yy}{Vector representing the response variable
}
\item{xx}{Design matrix of the covariates (including the intercept)
}
\item{var.function}{Assumed function for the variance. "\code{power}" function corresponds to  \eqn{\sqrt{Var}=\sigma=\phi|x^T \beta|^{\gamma}}{\sqrt(var) = \sigma = \phi |x^T \beta|^ \gamma} and "\code{exponential}" to \eqn{\sqrt{Var} = \sigma = \phi e^{\gamma |x^T \beta|}}{\sqrt(Var) = \sigma = \phi e^(\gamma |x^T \beta|)}.

}
\item{phi.par}{If \code{TRUE}, the function estimate the phi parameter. If \code{FALSE}, phi is assumed equal to 1.

}
\item{tuning.para}{If \code{NULL}, the function will run the estimation procedure for a range of value between 0 and 3 and will select the tuning parameter that minimizes the variance of the estimates. The user can also indicate a value of tuning parameter: in this case the estimation procedure will be evaluated once with the selected value of the tuning parameter.  

}
\item{step}{Only works when \code{tuning.para = NULL}, indicates the increment of the tuning parameter sequence (between 0 and 3) tested by the function. It will determine the precision of the tuning parameter. Caution : a smaller value indicates a larger number of value tested, resulting in a longer computing time.

}
\item{n.lag}{If \code{NULL}, a pAcf plot of the resiuduals will appear and you will have to indicate the number of lags the method has to include. The user can also give an integer corresponding to the number of lags desired.

}
\item{print.summary}{If \code{TRUE}, prints a summary of the estimates.
}


}
%\details{
%efgdgd
%}
\value{
The function returns a list including

\item{coefficients}{ Value of the robust estimates}

\item{residuals}{Residuals of the model.}

\item{p_residuals}{Pearson residuals of the model.}

\item{r_residuals}{Robust pearson residuals of the model : \eqn{\psi(p_residuals,c)}} with 
\eqn{\psi} the derivative of the loss function and c the chosen tuning parameter.

\item{fitted values}{Fitted values obtained with the robust method}

\item{vcov}{Variance-covariance matrix of the estimates}

\item{summary}{Summary of the model including: values, standard errors and z-values of the estimates}

\item{model}{Design matrix of the model}

\item{tuningpara}{When \code{tuning.para = NULL}, list containing the optimal tuning parameter, all the values of tuning parameter tested and their associated variance obtained.}

\item{varpara}{Estimates of the variance parameters}
}
\references{
Callens, A., Wang, Y-G., Fu, L. & Liquet, B. (2018). Robust estimation for autoregressive models with heterogeneity. Submitted.


}
\author{
Aurelien Callens, You-Gan Wang, Benoit Liquet
}


\seealso{
\code{rlm} function from  package \code{MASS}
}
\examples{

library(tseries)
data(ice.river)
xx <- model.matrix(flow.vat ~ prec + temp, data = ice.river)
yy <- flow.jok

least_square <- lm(flow.vat ~ prec + temp, data = ice.river)
pacf(least_square$residuals)
qqnorm(least_square$residuals)
qqline(least_square$residuals, col = "red", lwd = 2)

#With choice of optimal tuning parameter and 2 lags. 
#Note that if lag = NULL, a Pacf plot will appear to help you choose 
#the number of lags  and you will need to input this number in the console. 

model_1 <- rlmDD_het(yy, xx, var.function = "exponential", 
                                    tuning.para = NULL, n.lag = 2)

pacf(model_1$p_residuals)     
qqnorm(model_1$r_residuals)
qqline(model_1$r_residuals, col = "red", lwd = 2)

#For fixed number of lags and tuning parameter
model_2 <- rlmDD_het(yy, xx, var.function = "exponential", 
                                  tuning.para = 1.345, n.lag = 2)



}
\keyword{regression}
