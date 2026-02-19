
#ifndef ___PRIOR_H___
#define ___PRIOR_H___

#include <functional>
#include "matrix.h"
 

struct prior_type
{
  matrix mu;
  matrix sigma;
  prior_type& operator=(const prior_type&);
};

typedef std::function<prior_type (const int&)> prior_function_type ;


#endif
