

#ifndef ___THEOREM_4_H___
#define ___THEOREM_4_H___

#include "particle_type.h"
#include "time_type.h"
#include "real_type.h"
#include "dnorm.h"
#include "dlst.h"
#include "plurality_type.h"


template<typename noise, plurality_type plurality>
real_type predict_mu(const particle_type<noise,plurality>& particle_t,const time_type& t)
{

  // use temporary objects for now - optimise later
  matrix sigma_post = particle_t.post.sigma;
  matrix mu_post = particle_t.post.mu;
  matrix H_t_tau  = particle_t.feature_vector_function(t,particle_t.tau);

  // moving from matrices to scalars like this seems like code stink
  matrix temp = H_t_tau.transpose() * sigma_post * H_t_tau;
  real_type mu_pred = (H_t_tau.transpose()*mu_post)(0,0);
  return mu_pred;
}


template <typename noise_type, plurality_type plurality>
requires requires { requires std::same_as<noise_type,unknown_variance>; }
real_type predict_noise(const particle_type<noise_type,plurality>& particle_t,const time_type& t, const real_type& y)
{
  // use temporary objects for now - optimise later
  matrix sigma_post = particle_t.post.sigma;
  matrix mu_post = particle_t.post.mu;
  matrix H_t_tau  = particle_t.feature_vector_function(t,particle_t.tau);
  // moving from matrices to scalars like this seems like code stink
  matrix temp = H_t_tau.transpose() * sigma_post * H_t_tau;
  real_type mu_pred = predict_mu(particle_t,t);
  real_type lst_nu = particle_t.noise_structure.nu;
  real_type lst_iota = particle_t.noise_structure.iota;
  
  real_type var_pred = (lst_iota/lst_nu)*(1.0 + temp(0,0));
  return dlst(y,2*lst_nu,mu_pred,var_pred);
}

template <typename  noise_type, plurality_type plurality>
requires requires { requires std::same_as<noise_type,known_variance>; }
real_type predict_noise(const particle_type<noise_type,plurality>& particle_t,const time_type& t, const real_type& y)
{
  // use temporary objects for now - optimise later
  matrix sigma_post = particle_t.post.sigma;
  matrix mu_post = particle_t.post.mu;
  matrix H_t_tau  = particle_t.feature_vector_function(t,particle_t.tau);
  real_type sigma = particle_t.noise_structure.sigma;   
  // moving from matrices to scalars like this seems like code stink
  matrix temp = H_t_tau.transpose() * sigma_post * H_t_tau;
  real_type mu_pred = predict_mu(particle_t,t);
  real_type var_pred = sigma*sigma*(1.0 + temp(0,0));
  return dnorm(y,mu_pred,std::sqrt(var_pred));
}


template<typename noise, plurality_type plurality>
particle_type<noise,plurality>& theorem_4(particle_type<noise,plurality>& particle_t, const time_type& t, const real_type& y)
{
  particle_t.weight = particle_t.weight * predict_noise(particle_t,t,y); 
  return(particle_t);
}


#endif
