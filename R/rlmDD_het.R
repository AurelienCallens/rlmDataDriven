rlmDD_het <- function (yy, xx, var.function = c("power", "exponential"), 
                           phi.par = TRUE, tuning.para = NULL, step = 0.1, n.lag = NULL, 
                           print.summary = TRUE){
  Y <- as.matrix(yy)
  n <- length(Y)
  z <- c()
  rlm.model <- rlm(yy ~ xx - 1)
  X <- as.matrix(rlm.model$x)
  switch(var.function, power = {
    est <- function(x) {
      sum(rlmDataDriven::chi((Y - mu)/(phi * abs(mu)^x)) * log(abs(mu)))
    }
    est_lag <- function(x) {
      sum(rlmDataDriven::chi((Y[-c(1:nlag)] - mu)/(phi * abs(mu)^x)) * 
            log(abs(mu)))
    }
    var.func <- function(mu, esti, phi) {
      phi * abs(mu)^esti
    }
  }, exponential = {
    est <- function(x) {
      sum(rlmDataDriven::chi((Y - mu)/(phi * exp(x * abs(mu)))) * abs(mu))
    }
    est_lag <- function(x) {
      sum(rlmDataDriven::chi((Y[-c(1:nlag)] - mu)/(phi * exp(x * abs(mu)))) * 
            abs(mu))
    }
    var.func <- function(mu, esti, phi) {
      phi * exp(esti * abs(mu))
    }
  }, stop("Wrong function name"))
  if (is.null(tuning.para)) {
    c_i <- seq(from = 0.1, to = 3, by = step)
    indic <- c()
    if (print.summary) 
      pb <- txtProgressBar(min = 0, max = length(c_i), 
                           style = 3)
    for (i in 1:length(c_i)) {
      c <- c_i[i]
      mu <- rlm.model$fitted.values
      phi <- 1
      esti <- uniroot(est, c(-4, 4))$root
      if (phi.par == T) {
        phi <- median(abs(Y - mu)/(var.func(mu, esti, 
                                            phi = 1)))/0.6745
      }
      else {
        phi <- 1
      }
      pearson_res <- (Y - mu)/(var.func(mu, esti, phi))
      w_i <- as.vector((psi.Huber(pearson_res, c)/(pearson_res)))
      W <- pmin(as.vector(w_i/(var.func(mu, esti, phi)))^2, 
                1)
      B_est <- solve(crossprod(X, W * X), crossprod(X, 
                                                    W * Y))
      mu <- X %*% B_est
      pearson_res <- (Y - mu)/(var.func(mu, esti, phi))
      psip <- rlmDataDriven::dpsi.Huber(pearson_res, c)
      mn <- mean(psip)
      p <- as.numeric(length(B_est))
      K <- 1 + p/n * var(psip)/(mn^2)
      Kc <- (1/(n - p)) * sum(psi.Huber(pearson_res, c)^2)/mn^2 * 
        K^2
      indic[i] <- sum(diag(Kc * solve(t(X) %*% diag(as.vector(1/(var.func(mu, 
                                                                          esti, phi))))^2 %*% X)))
      if (print.summary) {
        setTxtProgressBar(pb, i)
      }
    }
    c <- c_i[which.min(indic)]
  }else{
    c <- tuning.para
  }
  mu <- rlm.model$fitted.values
  phi <- 1
  esti <- uniroot(est, c(-4, 4))$root
  if (phi.par == T) {
    phi <- median(abs(Y - mu)/(var.func(mu, esti, phi = 1)))/0.6745
  }else {
    phi <- 1
  }
  pearson_res <- (Y - mu)/(var.func(mu, esti, phi))
  w_i <- as.vector((psi.Huber(pearson_res, c)/(pearson_res)))
  W <- pmin(as.vector(w_i/(var.func(mu, esti, phi)))^2, 1)
  B_est <- solve(crossprod(X, W * X), crossprod(X, W * Y))
  mu <- X %*% B_est
  pearson_res <- (Y - mu)/(var.func(mu, esti, phi))
  if (is.null(n.lag)) {
    pacf(psi.Huber(pearson_res, c), main = "", ylim = c(0, 
                                                        1))
    cat("\n")
    cat("Select the number of lags: ")
    nlag <- as.integer(readLines(n = 1))
  }
  else {
    nlag <- n.lag
  }
  if (nlag == 0) {
    stop("Number of lag must be at least 1")
  }
  
  list.l <- create_lag(pearson_res, n.lag = nlag)
  
  for (a in 1:nlag) assign(paste0("lag", a), list.l[[a]])
  dat.l <- as.data.frame(list.l) * (var.func(mu, esti, phi))
  dat.ln <- matrix(c(rep(1, nlag)), ncol = 1)
  row.names(dat.ln) <- names(dat.l)
  B_est <- as.matrix(rbind(B_est, dat.ln))
  X <- as.matrix(cbind(X, dat.l))
  B_est <- solve(crossprod(X[-c(1:nlag), ], W[-c(1:nlag)] * 
                             X[-c(1:nlag), ]), crossprod(X[-c(1:nlag), ], W[-c(1:nlag)] * 
                                                           Y[-c(1:nlag)]))
  mu <- X[-c(1:nlag), ] %*% B_est

  pearson_res <- (Y[-c(1:nlag)] - mu)/(var.func(mu, esti, 
                                                phi))

  p <- length(B_est)
  psip <- dpsi.Huber(pearson_res, c)
  mn <- mean(psip)
  K <- 1 + p/n * var(psip)/(mn^2)
  CC <- solve(t(X[-c(1:nlag), ]) %*% diag(as.vector(1/(var.func(mu, 
                                                                esti, phi))))^2 %*% X[-c(1:nlag), ])
  vcov.mod <- (1/(n - p)) * sum(psi.Huber(pearson_res, c)^2)/mn^2 * 
    K^2 * CC
  sum_table <- data.frame(Estimate = B_est, Std.Error = sqrt(diag(vcov.mod)))
  sum_table$t_value <- sum_table$Estimate/sum_table$Std.Error
  sum_table <- round(sum_table[, -5], digits = 5)
  z$coefficients <- B_est
  z$residuals <- (Y[-c(1:nlag)] - mu)
  z$fitted.values <- mu
  z$vcov <- vcov.mod
  z$summary <- sum_table
  z$model <- cbind(yy, X)
  z$p_residuals <- pearson_res
  z$r_residuals <- psi.Huber((Y[-c(1:nlag)] - mu), c = c)
  z$tuningpara <- if (is.null(tuning.para)) {
    list(optimal = c, tested_values = c_i, obtained_var = indic)
  }
  else {
    list(c = c)
  }
  z$varpara <- data.frame(gamma = esti, phi = phi)
  if (print.summary) {
    cat("Summary : \n ")
    print(z$summary)
    if (is.null(tuning.para)) {
      cat("\n Optimal tuning parameter:", c)
    }
    else {
      cat("\n Tuning parameter:", c)
    }
  }
  return(z)
}
