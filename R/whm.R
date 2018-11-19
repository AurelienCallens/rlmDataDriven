whm <- function(yy, xx, var.function = "power", tuning.para = 1.345, ite = 5){
  
  c <- tuning.para
  
  Y <- as.matrix(yy)
  n <- length(Y) 
  z <- c()
  rlm.model <- rlm(yy ~ xx -1, k = c)
  X <- as.matrix(rlm.model$x) 
  
  switch(var.function,
         "power"={
           est <- function(x){sum(chi((Y-mu)/(phi*abs(mu)^x))*log(abs(mu)))}
           var.func <- function(mu, esti, phi){phi*abs(mu)^esti}
           est_phi <- function(x){sum(chi((Y-mu)/(x*abs(mu)^esti)))}
         },
         "exponential"={
           est <- function(x){sum(chi((Y-mu)/(phi*exp(x*abs(mu))))*abs(mu))}
           var.func <- function(mu, esti, phi){phi*exp(esti*abs(mu))}
           est_phi <- function(x){sum(chi((Y-mu)/(x*exp(esti*abs(mu)))))}
         },
         stop("Wrong function name")
  )
  
  
  
  mu <- rlm.model$fitted.values
  phi <- 1
  
  for(i in 1: ite){
    esti<- uniroot(est,c(-2,2))$root
    
    phi <- uniroot(est_phi,c(0,2))$root # median(abs(Y-mu)/(var.func(mu, esti, phi = 1)))/0.6745
    
    pearson_res <- (Y-mu) / (var.func(mu, esti, phi))
    
    w_i <- as.vector((psi.Huber(pearson_res, c)/(pearson_res)))
    
    W <- pmin(as.vector(w_i / (var.func(mu, esti, phi)))^2, 1)
    
    B_est <- solve(crossprod(X,W*X), crossprod(X,W*Y))
    
    mu <-  X %*% B_est
  }
  pearson_res <- (Y-mu) / (var.func(mu, esti, phi))
  
  
  p<- length(B_est)
  psip <- dpsi.Huber(pearson_res,c)
  mn <-mean(psip)
  K <- 1 + p/n * var(psip)/(mn^2)
  CC <- solve(t(X)%*%  diag(as.vector(1 / (var.func(mu, esti, phi))))^2 %*% X) 
  vcov.mod<-(1/(n-p))*sum(psi.Huber(pearson_res,c)^2)/mn^2*K^2*CC
  
  sum_table <- data.frame(Estimate = B_est, Std.Error=sqrt(diag(vcov.mod)))
  sum_table$t_value <- sum_table$Estimate/sum_table$Std.Error
  sum_table<- round(sum_table[,-5], digits = 5)
  
  z$coefficients <- B_est
  z$residuals <-  (Y-mu) 
  z$fitted.values <- mu
  z$vcov <- vcov.mod
  z$summary <-sum_table
  z$model <- cbind(yy,X)
  z$tuningpara <- list(c = c)
  z$varpara <-data.frame(gamma=esti, phi= phi)
  cat("Summary : \n ")
  print(z$summary)
  cat("\n Tuning parameter:",c)
  z
}
