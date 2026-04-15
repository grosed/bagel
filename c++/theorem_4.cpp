#include "theorem_4.h"
#include "dnorm.h"
#include "dlst.h"

#include <cmath>

particle_type& theorem_4(particle_type& particle_t, const time_type& t, const real_type& y)
{

  // use temporary objects for now - optimise later
  matrix sigma_post = particle_t.post.sigma;
  matrix mu_post = particle_t.post.mu;
  matrix H_t_tau  = particle_t.feature_vector_function(t,particle_t.tau);
  real_type sigma = particle_t.s; // maybe change s for sigma

  // moving from matrices to scalars like this seems like code stink
  matrix temp = H_t_tau.transpose() * sigma_post * H_t_tau;
  real_type mu_pred = (H_t_tau.transpose()*mu_post)(0,0);
  
  if(particle_t.known_variance == true)
    {
      real_type var_pred = sigma*sigma*(1.0 + temp(0,0));
      particle_t.weight = particle_t.weight * dnorm(y,mu_pred,std::sqrt(var_pred));
    }

  else
    {
      real_type lst_nu = particle_t.lst_nu;
      real_type lst_iota = particle_t.lst_iota;
      real_type var_pred = (lst_iota/lst_nu)*(1.0 + temp(0,0));
      particle_t.weight = particle_t.weight * dlst(y,2*lst_nu,mu_pred,var_pred);
    }
  
  return(particle_t);
}



