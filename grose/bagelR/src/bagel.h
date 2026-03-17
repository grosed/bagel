#ifndef ___BAGEL_H___
#define ___BAGEL_H___

#include "probability.h"
#include "time.h"
#include "real.h"
#include "tau.h"
#include "time.h"
#include "particle.h"
#include "theorem_2.h"
#include "theorem_3.h"
#include "theorem_1.h"
#include "theorem_4.h"
#include "dnorm.h"

#include <list>


struct bagel_type
{

  bagel_type(const prior_function_type&,
	     const feature_vector_function_type&,
	     const probability_type&,
	     const probability_type&,
	     const real_type&);

  bagel_type& update(const real_type&);
  real_type weight_0_t() const;

  
  // private:
  
  time_type t;
  probability_type p = 1.0;
  probability_type p0 = 0.9;
  real_type s = 1.0;
  particle_type initial_particle;
  std::list<particle_type> particles;

};

real_type weight_0_t(const bagel_type&);
bagel_type& update(bagel_type&, const real_type&);







#endif
