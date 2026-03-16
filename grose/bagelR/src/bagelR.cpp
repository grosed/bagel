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

// temporary testing using c++ version with example_2
#include "example_2.h"


bagelR::bagelR(const probability_type& p0,
	       const probability_type& p,
	       const real_type& s)
{
  sp_bagel = std::make_shared<bagel_type>(prior_example_2,feature_vector_example_2,p0,p,s);
  std::cout << "p0 is : " << p0 << std::endl;
  std::cout << "p is : " << p << std::endl;
  std::cout << "s is : " << s << std::endl;
}


real_type bagelR::update(const real_type& x)
{
  sp_bagel -> update(x);
  return sp_bagel -> weight_0_t();
}




/*
bagelR::bagelR(const probability_type& _data)
{
  data = _data;
}
*/

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
    // .constructor<const probability_type&,const probability_type&,const real_type&>()
  .constructor<probability_type,probability_type,real_type>()
  .method("doit", &bagelR::doit)
  .method("taus", &bagelR::taus)
  .method("feature_vectors", &bagelR::feature_vectors)
  .method("update", &bagelR::update)
;
}



