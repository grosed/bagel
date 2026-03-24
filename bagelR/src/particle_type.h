
#ifndef ___PARTICLE_TYPE_H___
#define ___PARTICLE_TYPE_H___


#include "prior_type.h"
#include "post_type.h"
#include "feature_vector_type.h"
#include "real_type.h"
#include "tau_type.h"


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

