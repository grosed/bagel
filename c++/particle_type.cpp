
#include "particle_type.h"



particle_type::particle_type() {};

particle_type& particle_type::operator=(const particle_type& other)
{
  prior_function = other.prior_function;
  feature_vector_function = other.feature_vector_function;
  tau = other.tau;
  p0 = other.p0;
  p = other.p;
  s = other.s;
  post = other.post;
  weight = other.weight;
  ratios = other.ratios;
  return *this;
}

particle_type::particle_type(const prior_function_type& _prior_function,
			     const feature_vector_function_type& _feature_vector_function,
			     const tau_type& _tau,
			     const real_type& _p0,
			     const real_type& _p,
			     const real_type& _s) 
{
  prior_function = _prior_function;
  feature_vector_function = _feature_vector_function;
  tau = _tau;
  p0 = _p0;
  p = _p;
  s = _s;
  weight = 1.0;
  ratios.push_front(1.0);
}


