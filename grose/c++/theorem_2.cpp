

#include "theorem_2.h"


particle_type theorem_2(const particle_type& particle_0, const time_type& t)
{

  particle_type particle_t = particle_0;

  if(t == 1)
    {
      particle_t.post = particle_post.prior_function(t);
      particle_t.tau= t - 1;
      return particle_t;
    }


  prior_type prior_1 = particle_0.prior_function(1);
  prior_type prior_t = particle_0.prior_function(t);
  
  int d1 = prior_1.mu.rows();
  int d2 = prior_1.sigma.rows();
  int n = d1 + d2;


  sigma_beta_beta_t <- prior_t.sigma.block(0,0,d1,d1);     
  sigma_beta_gamma_t <- prior_t.sigma.block(0,d1,d1,d2);
  sigma_gamma_beta_t <- prior_t.sigma.block(d1,0,d2,d1);      [(d1+1):n,1:d1]
   		sigma_gamma_gamma_t <- prior_t.sigma[(d1+1):n,(d1+1):n]


		  
   		mu_beta_t <- as.matrix(prior.mu.t[1:d1,1])
   		mu_gamma_t <- as.matrix(prior.mu.t[(d1+1):n,1])


  

   		sigma.beta.beta.t <- prior.sigma.t[1:d1,1:d1]
   		sigma.beta.gamma.t <- prior.sigma.t[1:d1,(d1+1):n]
   		sigma.gamma.beta.t <- prior.sigma.t[(d1+1):n,1:d1]
   		sigma.gamma.gamma.t <- prior.sigma.t[(d1+1):n,(d1+1):n]
   		mu.beta.t <- as.matrix(prior.mu.t[1:d1,1])
   		mu.gamma.t <- as.matrix(prior.mu.t[(d1+1):n,1])

  
  
  
  	        prior.1 <- object@prior(1L)
   		prior.mu.1 <- prior.1$mu
   		prior.sigma.1 <- prior.1$sigma
  


  
}

