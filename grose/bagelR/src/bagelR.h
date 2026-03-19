
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
  std::shared_ptr<bagel_type> sp_bagel; 
  std::map<int,matrix> M;
  std::map<int,prior_type> P;
  const matrix& feature_vector_from_R(const int&,const int&);
  const prior_type prior_from_R(const int&);
   
  bagelR(const probability_type&,const probability_type&,const real_type&);
  real_type update(const real_type&);
  std::list<int> get_taus();
  void set_feature_vectors(const std::vector<int>&, const std::list<matrix>&);
  double get_time();
  void set_priors(const std::vector<int>&,const std::list<matrix>&,const std::list<matrix>&);


};




#endif 
