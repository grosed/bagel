
#include <iostream>
#include "example_2.h"


int main()
{


prior_function_type  prior = prior_example_2;
 feature_vector_function_type feature_vector = feature_vector_example_2  ;
  
for(int t = 1; t <= 3; t++)
    {
      std::cout << "*****************************************" <<std::endl;
      //std::cout << prior(t).mu << std::endl;
      //std::cout << "------------------------------------" <<std::endl;
      //std::cout << prior(t).sigma << std::endl;
      //std::cout << "------------------------------------" <<std::endl;
      for(int tau = 0; tau <= t; tau++)
	{
	  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<std::endl;
	  std::cout << feature_vector(t,tau) << std::endl;
	}
    }
 
 return 0;
}




