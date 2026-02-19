
#ifndef ___PARTICLE_H___
#define ___PARTICLE_H___


#include "prior.h"
#include "post.h"
#include "feature_vector.h"
#include "real.h"
#include "tau.h"


struct particle_type
{
  prior_function_type prior_function;
  feature_vector_function_type feature_vector_function;
  tau_type tau;
  real_type p0;
  real_type p;
  real_type s;
  post_type post;
  real_type weight;
  particle_type();
  particle_type(const prior_function_type&,
		const feature_vector_function_type&,
		const tau_type&,
		const real_type&,
		const real_type&,
		const real_type&);
  particle_type& operator=(const particle_type&);
};





#endif

