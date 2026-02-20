#include "theorem_1.h"
#include "real.h"

particle_type& theorem_1(particle_type& particle_1,particle_type& particle_t, const time_type& t)
{
  real_type weight;
  if(particle_t.tau == particle_1.tau) // dealing with initil particle
    {
      if(t == 1)
	{
	  weight = particle_1.p0;
	}
      else
	{
	  weight = particle_t.weight;
	}
      particle_t.weight = weight;
      return particle_t;
    }
  if(particle_t.tau = t - 1)
    {
      weight = particle_1.weight*(1.0 - particle_t.p0)/(particle_t.p0 * real_type(t -1));
    }
  else
    {
      if(t == 2)
	{
	  weight = 1.0 - particle_t.p0;
	}
      else
	{
	  weight = particle_t.weight*real_type(t-2)/real_type(t-1);
	}
    }
  particle_t.weight = weight;
  return particle_t;
}
