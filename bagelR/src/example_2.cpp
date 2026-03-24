
#include "example_2.h"


prior_type prior_example_2(const int& t)
{
  prior_type prior;
  if(t == 1)
    {	
      prior.mu = matrix(2,1);
      prior.mu << 0.0, 0.0;
      prior.sigma = Eigen::MatrixXd::Identity(2,2);
    }
  else
    {
      prior.mu = matrix(4,1);
      prior.mu << 0.0, 0.0, 0.0, 0.0;
      double rt = double(t);
      matrix sigma_gamma_gamma(2,2);
      sigma_gamma_gamma << std::pow(rt-1,2), -(rt-1), -(rt-1), 1.0;
      prior.sigma = Eigen::MatrixXd::Identity(4,4);
      prior.sigma.block(2,2,2,2) = sigma_gamma_gamma;
    }
  return prior;
}

feature_vector_type feature_vector_example_2(const int& t, const int& tau)
{
  feature_vector_type h(4,1);
  double rt = double(t);
  h << 1.0, t, 1.0, t;
  if(t <= tau)
    {
      h.block(2,0,2,1) << 0.0, 0.0; 
    }
  if(tau == 0)
    {
      return h.block(0,0,2,1);
    }
  return h;
}










