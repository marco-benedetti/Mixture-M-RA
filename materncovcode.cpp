
#include <Rcpp.h>
#include <boost/math/special_functions/bessel.hpp>  // included in BH  

// [[Rcpp::depends(BH)]]    

using namespace Rcpp;

// [[Rcpp::export]]   
double besselKcpp(double x, double nu){
  return boost::math::cyl_bessel_k(nu,x);
}

//[[Rcpp::export]]
NumericMatrix maternCovcpp(NumericMatrix distmat,double sigma2,double phi,double nu){
  int m = distmat.nrow();
  int n = distmat.ncol();
  NumericMatrix Sigma(m,n);
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      if(distmat(i,j)==0){
        Sigma(i,j) = sigma2;
      }
      else{
        Sigma(i,j) = sigma2*((pow((distmat(i,j)/phi),nu)))*besselKcpp(distmat(i,j)/phi,nu)/(pow(2,(nu-1))*tgamma(nu));
      }
    }
  }
  return Sigma;
}


