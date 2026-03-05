#include <Rcpp.h>
using namespace Rcpp;

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]


#include "bagelR.h"


// [[Rcpp::export]]
int f(const int& x)
{
  return 2*x;
}


bagelR::bagelR(const double& _data)
{
  data = _data;
}

double bagelR::doit(const double& x)
{
  return data + x;
}


std::list<int> bagelR::taus()
{
  // just send a fixed number to test latency
  std::list<int> tau_values;
  for(int tau = 0; tau < 1000; tau++)
    {
      tau_values.push_front(tau);
    }
  return tau_values; 
}


int bagelR::feature_vectors(const std::vector<matrix>& fvs)
{
  int n = fvs.size();
  std::map<int,matrix> M;
  int tau = 0;
  for(auto& m : fvs)
    {
      M[tau] = m;
      tau = tau + 1;
    }
  
  for(int tau = 0; tau < n; tau++)
    {
      matrix k = M[tau]; 
    }
  
  return 0;
}

/*
  // just send a fixed number to test latency
int bagelR::feature_vectors(Rcpp::List& fvs)
{
  matrix m = Rcpp::as<matrix>(fvs[0]);
  
  int n = fvs.size();
  std::map<int,matrix> M;
  int tau = 0;
  for(auto& m : fvs)
    {
      M[tau] = m;
    }
  
    
  return 0;
}
*/




RCPP_MODULE(bagelR) 
{
  class_<bagelR >("bagelR")
  .constructor<double>()
  .method("doit", &bagelR::doit)
  .method("taus", &bagelR::taus)
  .method("feature_vectors", &bagelR::feature_vectors)
;
}



