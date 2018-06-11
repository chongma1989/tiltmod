#ifndef __UBMM__
#define __UBMM__

// We can now use the BH package
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/special_functions/trigamma.hpp>

using namespace Rcpp;
using namespace boost::math;

//[[Rcpp::plugins("cpp11")]]
////[[Rcpp::interfaces(r,cpp)]]

//' Fit a two-point mixture of Beta distributions
//' @param x A vector of numeric values
//' @param w A vector of two numeric values, representing the weights  
//'          of two Beta distributions. Default values are 0.5, respectively.
//' @param a0 Initial values of the alpha and beta for Beta distribution f0. 
//'        Default values are 1 and 1, respectively.
//' @param a1 Initial values of the alpha and beta for Beta distribution f1. 
//'        Default values are 0.5 and 0.5, respectively. 
//' @param precision The tolerance for convergence. Default value is 
//'        1e-6.
//' @param MaxIter The maximum iteration for the EM algorithm. Default value
//'        is 10000L.  
//' @return A list of four components, including the converged weight,
//'         parameters for Beta distribution f0, parameters for Beta distribution
//'         f1, and the convergence iteration, respectively. 
//'                  
//' @examples 
//' x0=rbeta(900,0.8,0.8)
//' x1=rbeta(100,0.2,0.2)
//' \dontrun{
//'    MBM(c(x0,x1),w=c(0.8,0.2),a0=c(1,1),a1=c(0.5,0.5))
//' }
//' 
//' @export
//[[Rcpp::export]]
Rcpp::List MBM(const NumericVector & x, 
               const NumericVector & w = NumericVector::create(), 
               const NumericVector & a0 = NumericVector::create(), 
               const NumericVector & a1 = NumericVector::create(),
               const double & precision = 1e-6,
               const int & MaxIter=10000) {
  //Initialization
  NumericVector old_w = w.size()==0 ? NumericVector::create(0.5,0.5) : w;
  NumericVector old_a0 = a0.size()==0 ? NumericVector::create(1,1) : a0;
  NumericVector old_a1 = a1.size()==0 ? NumericVector::create(0.5,0.5) : a1;
  NumericVector new_w(2),new_a0(2),new_a1(2);
  NumericVector logx=Rcpp::log(x),logxc=Rcpp::log(1-x);
 
  int itr=0;
  while(itr<MaxIter){
    NumericVector old_dbeta0=Rcpp::dbeta(x,old_a0[0],old_a0[1]);
    NumericVector old_dbeta1=Rcpp::dbeta(x,old_a1[0],old_a1[1]);
    NumericVector weight=old_w[1]*old_dbeta1/(old_w[0]*old_dbeta0+old_w[1]*old_dbeta1);
   
    //update weights
    new_w[0]=1-Rcpp::mean(weight);
    new_w[1]=1-new_w[0];
    
    //update parameters in Beta distribution
    new_a0[0]=old_a0[0]-((digamma(old_a0[0]+old_a0[1])-digamma(old_a0[0]))*new_w[0]+
      Rcpp::mean((1-weight)*logx))/((trigamma(old_a0[0]+old_a0[1])-trigamma(old_a0[0]))*new_w[0]);
    
    new_a0[1]=old_a0[1]-((digamma(new_a0[0]+old_a0[1])-digamma(old_a0[1]))*new_w[0]+
      Rcpp::mean((1-weight)*logxc))/((trigamma(new_a0[0]+old_a0[1])-trigamma(old_a0[1]))*new_w[0]);
    
    new_a1[0]=old_a1[0]-((digamma(old_a1[0]+old_a1[1])-digamma(old_a1[0]))*new_w[1]+
      Rcpp::mean(weight*logx))/((trigamma(old_a1[0]+old_a1[1])-trigamma(old_a1[0]))*new_w[1]);
    
    new_a1[1]=old_a1[1]-((digamma(new_a1[0]+old_a1[1])-digamma(old_a1[1]))*new_w[1]+
      Rcpp::mean(weight*logxc))/((trigamma(new_a1[0]+old_a1[1])-trigamma(old_a1[1]))*new_w[1]);
    
    NumericVector new_dbeta0=Rcpp::dbeta(x,new_a0[0],new_a0[1]);
    NumericVector new_dbeta1=Rcpp::dbeta(x,new_a1[0],new_a1[1]);
    double loldsum=Rcpp::sum(Rcpp::log(old_w[0]*old_dbeta0+old_w[1]*old_dbeta1));
    double lnewsum=Rcpp::sum(Rcpp::log(new_w[0]*new_dbeta0+old_w[1]*new_dbeta1));
    
    if(std::fabs(loldsum-lnewsum)<precision){
      break;
    }else{
      old_w=new_w;
      old_a0=new_a0;
      old_a1=new_a1;
      itr=itr+1;
    }
  }
  return Rcpp::List::create(Rcpp::Named("Weight")=new_w,
                            Rcpp::Named("BetaPar0")=new_a0,
                            Rcpp::Named("BetaPar1")=new_a1,
                            Rcpp::Named("Iterations")=itr);
}


//' Fit a mixture of uniform and Beta distribution
//' @param x A vector of numeric values
//' @param w A vector of two numeric values, representing the weights  
//'          of the uniform and Beta distributions. Default 
//'          values are 0.5, respectively.
//' @param a Initial values of the alpha and beta for the Beta
//'          distribution. Defaults are obtained from MOM estimators.
//' @param precision The tolerance for convergence. Default value is 
//'        1e-8.
//' @param MaxIter The maximum iteration for the EM algorithm. Default value
//'        is 10000L.  
//' @return A list of three components, including the converged weight,
//'         parameters for Beta distribution, and the convergence iteration, 
//'         respectively. 
//'                  
//' @examples 
//' x0=runif(900) 
//' x1=rbeta(100,0.5,0.5)
//' UBMM(c(x0,x1),w=c(0.8,0.2),a=c(0.7,0.8))
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List UBMM(const NumericVector & x, 
                const NumericVector & w = NumericVector::create(), 
                const NumericVector & a = NumericVector::create(), 
                const double & precision = 1e-8,
                const int & MaxIter=10000) {
  //Initialization
  NumericVector init_w(2), init_a(2);
  init_a[0]=Rcpp::mean(x)*(Rcpp::mean(x)*(1-Rcpp::mean(x))/Rcpp::var(x)-1);
  init_a[1]=init_a[0]*(1/Rcpp::mean(x)-1);
  init_w[0]=Rcpp::mean(1/(1+Rcpp::dbeta(x,init_a[0],init_a[1])));
  init_w[1]=1-init_w[0];
  
  NumericVector old_w = w.size()==0 ? init_w : w;
  NumericVector old_a = a.size()==0 ? init_a : a;
  NumericVector new_w(2),new_a(2);
  NumericVector logx=Rcpp::log(x),logxc=Rcpp::log(1-x);
  
  int itr=0;
  while(itr<MaxIter){
    NumericVector old_dbeta=Rcpp::dbeta(x,old_a[0],old_a[1]);
    NumericVector weight=old_w[1]*old_dbeta/(old_w[0]+old_w[1]*old_dbeta);
    
    //update weights
    new_w[0]=1-Rcpp::mean(weight);
    new_w[1]=1-new_w[0];
   
    //update parameters in Beta distribution
    new_a[0]=old_a[0]-((digamma(old_a[0]+old_a[1])-digamma(old_a[0]))*new_w[1]+
      Rcpp::mean(weight*logx))/((trigamma(old_a[0]+old_a[1])-trigamma(old_a[0]))*new_w[1]);
    
    new_a[1]=old_a[1]-((digamma(new_a[0]+old_a[1])-digamma(old_a[1]))*new_w[1]+
      Rcpp::mean(weight*logxc))/((trigamma(new_a[0]+old_a[1])-trigamma(old_a[1]))*new_w[1]);
    
    NumericVector new_dbeta=Rcpp::dbeta(x,new_a[0],new_a[1]);
    double loldsum=Rcpp::sum(Rcpp::log(old_w[0]+old_w[1]*old_dbeta));
    double lnewsum=Rcpp::sum(Rcpp::log(new_w[0]+old_w[1]*new_dbeta));
    
    if(std::fabs(loldsum-lnewsum)<precision){
      break;
    }else{
      old_a=new_a;
      old_w=new_w;
      itr=itr+1;
    }
  }
  return Rcpp::List::create(Rcpp::Named("Weight")=new_w,
                            Rcpp::Named("BetaPar")=new_a,
                            Rcpp::Named("Iterations")=itr);
}

#endif 
