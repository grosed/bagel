#ifndef ___EXAMPLES_H___
#define ___EXAMPLES_H___

#include "model_type.h"

prior_type prior_example_2(const int&);
feature_vector_type feature_vector_example_2(const int&, const int&);
matrix transformation_example_2(const time_type&);

prior_type prior_example_1(const int&);
feature_vector_type feature_vector_example_1(const int&, const int&);
matrix transformation_example_1(const time_type&);


prior_type prior_daily(const int&);
feature_vector_type feature_vector_daily(const int&, const int&);
matrix transformation_daily(const time_type&);






#endif
