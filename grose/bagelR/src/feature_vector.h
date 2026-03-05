

#ifndef ___FEATURE_VECTOR_H___
#define ___FEATURE_VECTOR_H___

#include <functional>
#include "matrix.h"

typedef matrix feature_vector_type;
typedef std::function<feature_vector_type (const int&,const int&)> feature_vector_function_type ;

#endif
