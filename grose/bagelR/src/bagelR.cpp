#include <Rcpp.h>
using namespace Rcpp;

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]


#include "bagelR.h"

#include <algorithm>




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

  // std::function<const double& (const int&)> G = std::bind(&B::F, b, std::placeholders::_1);

  // std::function<const matrix& (const int&,const int&)> G = std::bind(&bagelR::H, *this, std::placeholders::_1,std::placeholders::_2);

  feature_vector_function_type G =   std::bind(&bagelR::H, this, std::placeholders::_1,std::placeholders::_2);

  
  // sp_bagel = std::make_shared<bagel_type>(prior_example_2,feature_vector_example_2,p0,p,s);
  sp_bagel = std::make_shared<bagel_type>(prior_example_2,G,p0,p,s);
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
  std::list<int> ltaus;
  std::transform(sp_bagel->particles.begin(),
		 sp_bagel->particles.end(),
		 std::inserter(ltaus,ltaus.end()),
		 [](auto& particle){return particle.tau;});  
  return ltaus;  
}

const matrix& bagelR::H(const int& t,const int& tau)
{
  return M[tau];
}


int bagelR::feature_vectors(const std::vector<int>& taus, const std::list<matrix>& fvs)
{ 
  M.clear();
  std::transform(taus.begin(),
		 taus.end(),
		 fvs.begin(),
		 std::inserter(M,M.end()),
		 [](auto& tau,auto& m){return std::make_pair(tau,m);});
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



