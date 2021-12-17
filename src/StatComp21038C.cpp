#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param thin the number of between-sample random numbers
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' rnC <- gibbsC(100,10)
//' par(mfrow=c(2,1));
//' plot(rnC[,1],type='l')
//' plot(rnC[,2],type='l')
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N, int thin) {
  NumericMatrix mat(N, 2);
  double x = 1, y = 0.1;
  int a=1, b=5, n=10;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rbinom(1, n, y)[0];
      y = rbeta(1, x+a, n-x+b)[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}

//' @title A function that can sum the value of a vector
//' @description A function that can sum the value of a vector
//' @param x the vector
//' @param p the length of x or the length we want to caculate
//' @return a double value
//' @examples
//' \dontrun{
//' r<-as.numeric(-1,-2,3)
//' p<-length(r)
//' sumC(r,p)
//' }
//' @export
// [[Rcpp::export]]
double sumC (NumericVector x, int p){
  double total = 0;
  for(int i=0; i < p; i++){
    total += x[i];
  }
  return total;
}

//' @title The loglikelihood function of lambda
//' @description The loglikelihood function of lambda
//' @param lambda a double value
//' @param tauinv2 the vector
//' @param p the length of tauinv2
//' @return a double value
//' @examples
//' \dontrun{
//' lambda<-1;tauinv2<-c(1,2,3);p<-length(tauinv2)
//' loglikeli(lambda,tauinv2,p)
//' }
//' @export
// [[Rcpp::export]]
double loglikeli(double lambda, NumericVector tauinv2, int p) {
    return(p*log(pow(lambda,2))-pow(lambda,2)/2*sumC(1/tauinv2,p));
}

//' @title Empirical Bayes method for lambda
//' @ description it is an empirical bayes method to get the lambda, and we use mento carlo EM
//' @param L0 a double value
//' @param tol a double value
//' @param tauinv2 a vector
//' @param p an int value
//' @param lambda a double value
//' @return a numericvector
//' @examples
//' \dontrun{
//' L0<-0.1;tol<-10^(-3);tauinv2<-c(1,2,3);p<-3;lambda<-1
//' EmpiricalBayes(L0,tol,tauinv2,p,lambda)
//' }
//' @export
// [[Rcpp::export]]
NumericVector EmpiricalBayes(double L0, double tol, NumericVector tauinv2, int p, double lambda){
    NumericVector value(2);
    double lambda1;
    double lambda2;
    double L1;
    double L2;
    lambda1 = sqrt((p*2)/sumC(1/tauinv2,p));
    L1 = loglikeli(lambda,tauinv2,p);
    if(abs(L0-L1)>tol){
        lambda2 = lambda1;
        L2 = L1;
    }
    else{
        lambda2 = lambda;
        L2 = L0;
    }
    value[0]=lambda2;
    value[1]=L2;
    return(value);
}


//' @title Hyperprior method for lambda
//' @description Hyperprior method for lambda
//' @param tauinv2 a vector
//' @param r an int value
//' @param d an int value
//' @param p an int vale
//' @return a double value
//' @examples
//' \dontrun{
//' tauinv2<-c(1,2,3);r<-1;d<-1;p<-3
//' HyperpriorBayes(tauinv2,r,d,p)
//' }
//' @export
// [[Rcpp::export]]

double HyperpriorBayes(NumericVector tauinv2, int r, int d, int p){
    double lambda2;
    int sh;
    double sc;
    sh = p+r;
    sc = sumC(1/tauinv2,p)/2 + d;
    lambda2 =sqrt(rgamma(1, sh, sc)[0]);
    return(lambda2);
}
























