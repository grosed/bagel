
#ifndef ___BAGELR_H___
#define ___BAGELR_H___


#include "matrix.h" 
#include <list>
#include <vector>
#include <map>
#include "bagel.h"
#include <memory>

struct bagelR
{
  double data;
   
  std::shared_ptr<bagel_type> sp_bagel; 


  std::map<int,matrix> M;
  const matrix& H(const int&,const int&);
  
  
  bagelR(const probability_type&,const probability_type&,const real_type&);
  // bagelR(const probability_type&);
  double doit(const double&);
  real_type update(const real_type&);
  std::list<int> taus();
  // int feature_vectors(Rcpp::List&);
  // int feature_vectors(const std::vector<matrix>&);
  int feature_vectors(const std::vector<int>&, const std::list<matrix>&);

  
};




#endif 
