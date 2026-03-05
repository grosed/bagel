
#ifndef ___BAGELR_H___
#define ___BAGELR_H___

#include "matrix.h" 
#include <list>
#include <vector>
#include <map>

struct bagelR
{
  double data;
  bagelR(const double&);
  double doit(const double&);
  std::list<int> taus();
  // int feature_vectors(Rcpp::List&);
  int feature_vectors(const std::vector<matrix>&);
  
};




#endif 
