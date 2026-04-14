#include <Rcpp.h>
using namespace Rcpp;

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]


#include "bagelR.h"
#include <algorithm>

// temporary testing using c++ version with example_2
#include "example_2.h"


bagelR::bagelR(const probability_type& p0,
	       const probability_type& p,
	       const real_type& s,
	       const int& n)
{
  feature_vector_function_type G_feature_vector =   std::bind(&bagelR::feature_vector_from_R, this, std::placeholders::_1,std::placeholders::_2);
  prior_function_type G_prior =   std::bind(&bagelR::prior_from_R, this, std::placeholders::_1);
  sp_bagel = std::make_shared<bagel_type>(G_prior,G_feature_vector,p0,p,s,n);
}


real_type bagelR::update(const real_type& x)
{
  sp_bagel -> update(x);
  return sp_bagel -> weight_0_t();
}


std::list<double> bagelR::get_weights()
{

  std::list<double> lweights;
  std::transform(sp_bagel->particles.begin(),
		 sp_bagel->particles.end(),
		 std::inserter(lweights,lweights.end()),
		 [](auto& particle){return particle.weight;});  
  return lweights;  
}
  
double bagelR::get_time()
{
  return sp_bagel -> t;
}


std::list<int> bagelR::get_taus()
{
  std::list<int> ltaus;
  std::transform(sp_bagel->particles.begin(),
		 sp_bagel->particles.end(),
		 std::inserter(ltaus,ltaus.end()),
		 [](auto& particle){return particle.tau;});  
  return ltaus;  
}

const matrix& bagelR::feature_vector_from_R(const int& t,const int& tau)
{
  return M[tau];
}

const prior_type bagelR::prior_from_R(const int& t)
{
  return P[t];
}


void bagelR::set_feature_vectors(const std::vector<int>& taus_from_R, const std::list<matrix>& feature_vectors_from_R)
{ 
  M.clear();
  std::transform(taus_from_R.begin(),
		 taus_from_R.end(),
		 feature_vectors_from_R.begin(),
		 std::inserter(M,M.end()),
		 [](auto& tau,auto& m){return std::make_pair(tau,m);});
}

void bagelR::set_priors(const std::vector<int>& ts_from_R,
			const std::list<matrix>& prior_mus_from_R,
			const std::list<matrix>& prior_sigmas_from_R)
{
  P.clear();
  
  // zip into prior_types
  std::list<prior_type> priors;
  std::transform(prior_mus_from_R.begin(),
		 prior_mus_from_R.end(),
		 prior_sigmas_from_R.begin(),
		 std::inserter(priors,priors.end()),
		 [](auto& mu,auto& sigma){prior_type prior;
		                          prior.mu = mu;
					  prior.sigma = sigma;
					  return prior;});
  std::transform(ts_from_R.begin(),
		 ts_from_R.end(),
		 priors.begin(),
		 std::inserter(P,P.end()),
		 [](auto& t,auto& prior){return std::make_pair(t,prior);});
}



std::list<std::list<double> > bagelR::get_ratios()
{
  return sp_bagel -> ratios();
}


RCPP_MODULE(bagelR) 
{
  class_<bagelR >("bagelR")
    .constructor<probability_type,probability_type,real_type,int>()
  .method("get_time", &bagelR::get_time)
  .method("get_taus", &bagelR::get_taus)
  .method("get_weights", &bagelR::get_weights)
  .method("set_feature_vectors", &bagelR::set_feature_vectors)
  .method("set_priors", &bagelR::set_priors)
  .method("update", &bagelR::update)
  .method("get_ratios", &bagelR::get_ratios)
;
}



