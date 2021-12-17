## -----------------------------------------------------------------------------
library(StatComp21038)
data("Data")
attach(Data)

## ----eval=FALSE---------------------------------------------------------------
#  double sumC (NumericVector x, int p){
#    double total = 0;
#    for(int i=0; i < p; i++){
#      total += x[i];
#    }
#    return total;
#  }

## ----eval=FALSE---------------------------------------------------------------
#  double loglikeli(double lambda, NumericVector tauinv2, int p) {
#      return(p*log(pow(lambda,2))-pow(lambda,2)/2*sumC(1/tauinv2,p));
#  }

## ----eval=FALSE---------------------------------------------------------------
#  NumericVector EmpiricalBayes(double L0, double tol, NumericVector tauinv2, int p, double lambda){
#      NumericVector value(2);
#      double lambda1;
#      double lambda2;
#      double L1;
#      double L2;
#      lambda1 = sqrt((p*2)/sumC(1/tauinv2,p));
#      L1 = Q(lambda,tauinv2,p);
#      if(abs(L0-L1)>tol){
#          lambda2 = lambda1;
#          L2 = L1;
#      }
#      else{
#          lambda2 = lambda;
#          L2 = L0;
#      }
#      value[0]=lambda2;
#      value[1]=L2;
#      return(value);
#  }

## ----eval=FALSE---------------------------------------------------------------
#  double HyperpriorBayes(NumericVector tauinv2, int r, int d, int p){
#      double lambda2;
#      int sh;
#      double sc;
#      sh = p+r;
#      sc = sumC(1/tauinv2,p)/2 + d;
#      lambda2 =sqrt(rgamma(1, sh, sc)[0]);
#      return(lambda2);
#  }

## ----eval = TRUE--------------------------------------------------------------
BayesianLasso<-function(x,y,center = T, scale = T,n.max=10000,E=TRUE,r=1,d=1){
  #beginning
  require("MASS")    # for the inverse of matrix
  require("mvtnorm") #multivarate normal
  require("pscl")    #inverse gamma
  require("statmod") #inverse gaussian 
  require("Rcpp")    #for we can use cpp function in the function
  #sourceCpp(paste0("//Users/apple/Desktop/StatComp/src/","StatComp.cpp"))
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  #centring and scaling
  ybar<-mean(y)
  yc<-y-ybar
  if(center==T){
    if(scale==T){
      xc<-x
    }
    else{
      xc<-scale(x,scale=T)
    }
  }
  else{
    if(scale==T){
      xc<-scale(x,center = T,scale = F)
    }
    else{
      xc<-scale(x,center = T,scale = F) 
    }
  }
  xbar<-apply(x,2,mean)
  xtx<-t(xc)%*%xc
  xty<-t(xc)%*%yc
  # initial
  beta <- drop(ginv(xtx)%*%xty)
  resid <- yc - xc %*% beta
  sigma2 <- (t(resid) %*% resid)/(n-p)
  tauinv2 <- 1 / (beta * beta)
  lambda <- p * sqrt(sigma2) / sum(abs(beta))
  #the loglike function prepare for the mento carlo EM
  L0 <- p*log(lambda^2)-lambda^2/2*sum(1/tauinv2) 
  Beta <- matrix(0, n.max, p)
  Sigma2 <- numeric(n.max)
  Tau2 <- matrix(0, n.max, p)
  Lambda <- numeric(n.max)
  # gibbs sampling
  for(i in 1:n.max){
    Dinv <- diag(tauinv2)
    Ainv <- ginv(xtx + Dinv)
    beta.m <- Ainv %*% xty
    Sigma <- as.numeric(sigma2) * Ainv
    # update beta
    if(det(Sigma)<0.001){
      beta <- rmvnorm(1, beta.m, Sigma, "svd") # for we can consider the matrix is singular
    } else{
      beta <- rmvnorm(1, beta.m, Sigma, "chol")
    }
    Beta[i,]<-beta
    # update sigma2
    sig.a <- (n-1)/2+p/2
    resid <- yc - xc %*% t(beta)
    sig.b <- t(resid) %*% resid/2 + beta%*% Dinv%*% t(beta)/2
    sigma2 <- rigamma(1, alpha=sig.a, beta=sig.b) 
    Sigma2[i] <- sigma2
    # update tau2
    beta<-as.vector(beta)
    mu.t<-numeric(p)
    for(k in 1:p){
      mu.t[k]<-sqrt(lambda^2 * sigma2 / beta[k]^2)
    }
    lambda.t <- lambda^2
    tauinv2 <- rinvgauss(p, mean=mu.t, shape=lambda.t)
    Tau2[i, ] <- 1/tauinv2
    if(E==T){
      tol<-10^(-3)
      lambda<-EmpiricalBayes(L0,tol,tauinv2,p,lambda)[1]
      L0<-EmpiricalBayes(L0,tol,tauinv2,p,lambda)[2]
    }
    else{
      lambda<-HyperpriorBayes(tauinv2,r,d,p)
    }
    Lambda[i]<-lambda
  }
  
  
  list(beta=round(apply(Beta[seq(round(n.max/2), n.max),],2, median),3),
       tau2=round(apply(Tau2[seq(round(n.max/2), n.max),],2,  median),3),
       sigma2=round(median(Sigma2[seq(round(n.max/2), n.max)]),3),
       lambda=round(median(Lambda[seq(round(n.max/2), n.max)]),3))
}

## -----------------------------------------------------------------------------
result<-BayesianLasso(Data$diabetes.x,Data$diabetes.y,T,T)
proc.time()

## -----------------------------------------------------------------------------
resid<-Data$diabetes.y-mean(Data$diabetes.y)-Data$diabetes.x%*%result$beta
n<-nrow(Data$diabetes.x)
t(resid)%*%resid/n

## -----------------------------------------------------------------------------
library(monomvn)
n.max=1e4
result1<-blasso(Data$diabetes.x,Data$diabetes.y,T=10000)
proc.time()
result2<-round(apply(result1$beta[seq(round(n.max/2), n.max),],2, median),3)
result2

## -----------------------------------------------------------------------------
resid<-Data$diabetes.y-mean(Data$diabetes.y)-Data$diabetes.x%*%result2
n<-nrow(Data$diabetes.x)
t(resid)%*%resid/n

## -----------------------------------------------------------------------------
Bayespredmse<-function(x,y,x1=as.matrix(0,1,1), center = T, scale = T,
                                        n.max=10000,E=TRUE,r=1,d=1,pred = T){
  if(pred==T){
    beta<-BayesianLasso(x,y,center,scale,n.max,E,r,d)$beta
    if(center==T){
      if(scale==T){
        xc<-x1
      }
      else{
        xc<-scale(x1,scale=T)
      }
    }
    else{
      if(scale==T){
        xc<-scale(x1,center = T,scale = F)
      }
      else{
        xc<-scale(x1,center = T,scale = F) 
      }
    }
    ypred<-xc%*%beta
    return(ypred)
  }
  else{
    n<-nrow(x)
    yc<-y-mean(y)
    result<-BayesianLasso(x,y,center,scale,n.max,E,r,d)
    if(center==T){
      if(scale==T){
        xc<-x
      }
      else{
        xc<-scale(x,scale=T)
      }
    }
    else{
      if(scale==T){
        xc<-scale(x,center = T,scale = F)
      }
      else{
        xc<-scale(x,center = T,scale = F) 
      }
    }
    resid<-yc-xc%*%result$beta
    mse1<-t(resid)%*%resid/n
    return(mse1)
  }
}


## -----------------------------------------------------------------------------
x<-Data$diabetes.x
y<-Data$diabetes.y
Bayespredmse(x,y,pred = F)

