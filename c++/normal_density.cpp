
#include "normal_density.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>

real_type normal_density(const real_type& x,const real_type& mu,const real_type& sd)
{
  boost::math::normal Z; // standard normal distribution
  return boost::math::pdf(Z,(x - mu)/sd)/sd;
}

real_type log_normal_density(const real_type& x,const real_type& mu,const real_type& sd)
{
  real_type z = (x - mu)/sd;
  return -0.5*z*z - 0.5*std::log(2.0*boost::math::constants::pi<real_type>()) - std::log(sd);
}

