
#include "dnorm.h"
#include <cmath>

real_type dnorm(const matrix& x, const matrix& mu, const matrix& sigma)
{
    int n = x.size();
    matrix diff = x - mu;
    real_type exponent = -0.5 * (diff.transpose() * sigma.inverse() * diff)(0, 0);
    real_type coeff = 1.0 / (std::pow(2 * M_PI, n / 2.0) * std::sqrt(sigma.determinant()));
    return coeff * std::exp(exponent);
}


