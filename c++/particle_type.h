
#ifndef ___PARTICLE_TYPE_H___
#define ___PARTICLE_TYPE_H___


#include "prior_type.h"
#include "post_type.h"
#include "feature_vector_type.h"
#include "real_type.h"
#include "tau_type.h"
#include "ratio_type.h"
#include "noise_type.h"
#include "plurality_type.h"
#include <list>




template<typename noise, plurality_type plurality>
struct particle_type
{
  noise noise_structure;
  prior_function_type prior_function;
  feature_vector_function_type feature_vector_function;
  tau_type tau;
  real_type p0;
  real_type p;
  real_type s;
  real_type lst_nu;
  real_type lst_iota;
  post_type post;
  real_type weight;
  std::list<ratio_type> ratios;
  bool known_variance;

  particle_type() {};


  particle_type& operator=(const particle_type& other)
    {
      prior_function = other.prior_function;
      feature_vector_function = other.feature_vector_function;
      noise_structure = other.noise_structure;
      tau = other.tau;
      p0 = other.p0;
      p = other.p;
      s = other.s;
      lst_nu = other.lst_nu;
      lst_iota = other.lst_iota;
      known_variance = other.known_variance;
      post = other.post;
      weight = other.weight;
      ratios = other.ratios;
      return *this;
    }

    particle_type(const prior_function_type& _prior_function,
		  const feature_vector_function_type& _feature_vector_function,
		  const noise& _noise_structure,
		  const tau_type& _tau,
		  const real_type& _p0,
		  const real_type& _p,
		  const real_type& _s,
		  const real_type& _lst_nu,
		  const real_type& _lst_iota,
		  const bool& _known_variance) 
  {
    prior_function = _prior_function;
    feature_vector_function = _feature_vector_function;
    noise_structure = _noise_structure;
    tau = _tau;
    p0 = _p0;
    p = _p;
    s = _s;
    lst_nu = _lst_nu;
    lst_iota = _lst_iota;
    known_variance = _known_variance;
    weight = 1.0;
    ratios.push_front(1.0);
  }

  
};












#endif

