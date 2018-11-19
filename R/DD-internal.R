##huber's function
psi.Huber<-function(r,c)
{
  newr<-pmin(c,pmax(r,-c))
  newr
}

##derivative of huber's function
dpsi.Huber<-function(r,c)
{
  newr<-as.numeric(abs(r)<=c)
  newr
}


##bisquare's function
psi.Tukey<-function(r,c)
{
  r<-pmin(c,pmax(r,-c))
  newr<-r*(1-(r/c)^2)^2
  newr
}

##derivative of bisquare's function
dpsi.Tukey<-function(r,c)
{
  r<-pmin(c,pmax(r,-c))
  newr<-(1-(r/c)^2)*(1-5*(r/c)^2)
  newr
}

##modified huber's function
psi.Exp<-function(r,c)
{
  newr<-r*exp(-(r/c)^2)
  newr
}

##derivative of modified huber's function
dpsi.Exp<-function(r,c)
{
  newr<-(1-(2*r^2)/c^2)*exp(-(r/c)^2)
  newr
}


##find the most efficient tuning constant
eff<-function(r,method,plot)
{
  nn<-  length(r)
  if (method =="Huber")
  {
    tau<-(1:30)/10
    new<-unlist(lapply(tau,FUN=function(x){(sum(dpsi.Huber(r,x)))^2/(nn*sum(psi.Huber(r,x)^2))}))
  }
  else if (method=="Bisquare")
  {
    tau<-(15.48:59.48)/10
    new<-unlist(lapply(tau,FUN=function(x){(sum(dpsi.Tukey(r,x)))^2/(nn*sum(psi.Tukey(r,x)^2))}))
  }
  else if (method=="Exponential")
  {
    tau<-(5:100)/10
    new<-unlist(lapply(tau,FUN=function(x){(sum(dpsi.Exp(r,x)))^2/(nn*sum(psi.Exp(r,x)^2))}))
  }

  if (plot =="Y")
  {
    plot(new, type="o", xaxt = "n", xlab="Tunning parameter",
        ylab="Efficiency",ylim=c(0,max(new)*1.3))
    axis(1, at=seq_along(tau), labels=tau)
    if (method =="Exponential") {
      legend(length(tau)*0.8,max(new)*1.3,c("Exp"),bty = "n")
    } else {
      legend(length(tau)*0.8,max(new)*1.3,c(method),bty = "n")
    }
  }

  xx<- order(new)
  return(tau[xx[length(tau)]])

}


rho.h<-function(u,c){
  phi<-0;
  x1<-(abs(u)<=c)
  phi<-x1* (u^2/2) + (1-x1)*(abs(u)*c-c^2/2)
  phi
}

rho.b<-function(u,c){
  phi<-0;
  x1<-(abs(u)<=c)
  phi<-x1*(1-(1-(u/c)^2)^3) + (1-x1)
  phi
}

rho.e<-function(u,c){
  1-exp(-(u/c)^2)
}


ESL_O<-function(x,xx,y,beta,newc,maxit=500, toler=1e-6)
  {
  it=0;delta=1;
  while (delta>toler && it<maxit){
    it<-it+1;
    beta0<-beta;
    w<-0;
    e<-y-x%*%beta0;
    S_n<-max(mad(e),1e-6)
    u<-e/S_n
    for(i in seq_along(u)){
      w[i]<-exp(-(u[i]/newc)^2);
    }

    W<-diag(w);
    #beta<-solve(t(x)%*%W%*%x)%*%t(x)%*%W%*%y;
    rlm1<-rlm(y~xx,psi=psi.huber,weights=w,k=newc)
    beta<-rlm1$coef
    r2<-summary(lm(y~xx,weights=w))$r.squared
    delta<-sqrt(sum((beta-beta0)^2));
  }
  list(esti=rlm1, Std.Error=summary(rlm1)$coef[,2], weights=w, tunning=newc, R2=r2)
}


chi<- function(obj){
  pmin(obj^2/1.041^2,1)-0.5
}

create_lag <- function(vec, n.lag){
  list.lag<-list()
  for(i in 1: n.lag){
    lagged.vec <- c(rep(NA,i),vec)[1:length(vec)]
    list.lag[[paste0("lag", i)]] <- lagged.vec
  }
  return(list.lag)
  
}

