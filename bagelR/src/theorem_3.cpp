
#include "theorem_3.h"
#include "matrix_type.h"
#include "particle_type.h"
#include "tau_type.h"


// modifies exiting particles - mutating
particle_type& theorem_3(particle_type& particle_t, const time_type& t, const real_type& val)
{
  matrix y(1,1);
  y(0,0) = val;
  prior_type prior_t = particle_t.prior_function(t);
  matrix mu = particle_t.post.mu;
  matrix sigma = particle_t.post.sigma;
  tau_type tau = particle_t.tau;
  matrix H_t = particle_t.feature_vector_function(t,tau);
  matrix I = matrix::Identity(1,1);
  matrix e = y - H_t.transpose() * mu;  
  matrix Q = H_t.transpose() * sigma * H_t + I;
  real_type rQ = Q(0,0);
  matrix A = sigma * H_t / rQ;
  sigma = sigma - A * A.transpose() * rQ;
  mu = mu + A * e;
  particle_t.post.mu = mu;
  particle_t.post.sigma = sigma;
  return particle_t;
}



