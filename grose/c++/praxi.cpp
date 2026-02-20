
#include <iostream>
#include "example_2.h"

#include "bagel.h"

int main()
{


  time_type t = 1;
  tau_type tau = 0;
  probability_type p = 1.0;
  probability_type p0 = 0.9;
  real_type s = 1.0;
  
  particle_type particle_1(prior_example_2,feature_vector_example_2,tau,p0,p,s);

  particle_1 = theorem_2(particle_1,t);

  real_type y = 1.262954;

  
  particle_1 = theorem_3(particle_1,t,y);

   
  t = t + 1;
  y = -0.3262334;
  particle_type particle_2 =  theorem_2(particle_1,t);
  particle_1 = theorem_3(particle_1,t,y);
  particle_2 = theorem_3(particle_2,t,y);

 
  t = t + 1;
  y = 1.329799;
  particle_type particle_3 =  theorem_2(particle_1,t);
  particle_1 = theorem_3(particle_1,t,y);
  particle_2 = theorem_3(particle_2,t,y);
  particle_3 = theorem_3(particle_3,t,y);

  std::cout << particle_1.post.mu << std::endl;
  std::cout << particle_1.post.sigma << std::endl;
  std::cout << particle_2.post.mu << std::endl;
  std::cout << particle_2.post.sigma << std::endl;
  std::cout << particle_3.post.mu << std::endl;
  std::cout << particle_3.post.sigma << std::endl;
  

  

  




  
  

  
  /*
  
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

  */

    
 return 0;
}




