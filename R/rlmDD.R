
rlmDD <-function(yy,xx,beta0,betaR,method =c("Huber","Bisquare","Exponential"),plot)
{
  X <- matrix(c(rep(1, length(yy)), xx), ncol = dim(xx)[2]+1)
  u <- yy-X%*%beta0
  sigma0 <- median(abs(u[!is.na(u)]))/0.6745      # drop NA
  uu <- yy-X%*%betaR
  ri <- (uu[!is.na(uu)])/sigma0            # drop NA
  newc <- eff(ri, method, plot)


  ## MM estimation of the location for R2 calculation
  rlm0 <- rlm(yy~1, method = "MM")
  mu <- rlm0$coef

  method <- match.arg(method)

  if (method =="Huber"){
    rlm1 <- rlm(yy~xx, psi = psi.huber, k = newc)
    sig <- rlm1$s
    A <- rho.h((yy-mu)/sig, newc)
    B <- rho.h(rlm1$res/sig, newc)
    r2 <- 1-sum(B)/sum(A)
    list(esti = rlm1, Std.Error = summary(rlm1)$coef[,2],
         tunning = newc, R2 = r2)
  }
  else if (method=="Bisquare"){
    rlm1 <- rlm(yy~xx, psi = psi.bisquare, c=newc)
    sig <- rlm1$s
    A <- rho.b((yy-mu)/sig,newc)
    B <- rho.b(rlm1$res/sig,newc)
    r2 <- 1-sum(B)/sum(A)
    list(esti = rlm1, Std.Error = summary(rlm1)$coef[,2],
         tunning = newc, R2 = r2)
  }
  else if (method=="Exponential")
  {
    m <- rlm(yy~xx, method="MM", maxit = 1000)
    re_mm <- data.frame(m[1])$coefficients
    X <- cbind(1,xx)
    beta <- re_mm
    ESL_O(X, xx, yy, re_mm, newc, maxit=500, toler=1e-6)
  }
}
