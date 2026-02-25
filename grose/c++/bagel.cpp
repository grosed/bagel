
#include "bagel.h"


bagel_type::bagel_type(const prior_function_type& prior_function,
		       const feature_vector_function_type& feature_vector_function,
		       const probability_type& p0,
		       const probability_type& p,
		       const real_type& s)
{
  t = 1;
  initial_particle = particle_type(prior_function,feature_vector_function,0,p0,p,s);
}


real_type bagel_type::weight_0_t() const
{
  return particles.front().weight;
}

real_type weight_0_t(const bagel_type& bagel)
{
  return bagel.weight_0_t();
}


bagel_type& update(bagel_type& bagel, const real_type& y)
{
  return bagel.update(y);
}

bagel_type& bagel_type::update(const real_type& y)
{
  if(particles.size() == 1000)
    {
      auto it = particles.begin();
      advance(it, 1);
      particles.erase(it);
    }
  
  if(t == 1)
    {
      initial_particle = theorem_2(initial_particle,t);
      initial_particle = theorem_1(initial_particle,initial_particle,t);
      initial_particle = theorem_4(initial_particle,t,y);
      initial_particle = theorem_3(initial_particle,t,y);
      initial_particle.weight = p0; 
      particles.push_back(initial_particle);
    }
  else
    {
      // new particle
      initial_particle = particles.front();
      particles.push_back(theorem_2(initial_particle,t));
      real_type sum_of_weights = 0.0;      
      for(auto& p : particles)
	{
	  p = theorem_1(initial_particle,p,t);
	}
      
      for(auto& p : particles)
	{
	  p = theorem_4(p,t,y);
	  p = theorem_3(p,t,y);
	  sum_of_weights += p.weight;
	}
      for(auto& p : particles)
	{
	  p.weight = p.weight/sum_of_weights;
	}	    
    }
  t = t + 1;
  return *this;
}

