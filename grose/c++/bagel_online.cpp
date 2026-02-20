#include <iostream>
#include "example_2.h"
#include "bagel.h"
#include <iostream>
#include <string>
#include <list>


using namespace std;

std::string input_line;

int main()
{


  time_type t = 1;
  tau_type tau = 0;
  probability_type p = 1.0;
  probability_type p0 = 0.9;
  real_type s = 1.0;

  // initialisation
  std::list<particle_type> particles;
  particle_type initial_particle(prior_example_2,feature_vector_example_2,tau,p0,p,s);
  initial_particle = theorem_2(initial_particle,t);
  try
    {
    while(cin)
      {
	// get the next value from the data stream
        getline(cin, input_line);
	double y = std::stod(input_line);
	// std::cout << y << std::endl;
	if(t == 1)
	  {
	    initial_particle = theorem_3(initial_particle,t,y);
	  }
	else
	  {
	    // new particle
	    particles.push_back(theorem_2(initial_particle,t));
	    // update particles
	    initial_particle = theorem_3(initial_particle,t,y);
	    for(auto& p : particles)
	      {
		p = theorem_3(initial_particle,t,y);
	      }
	  }
	t = t + 1;
    };
    }
  catch(...)
    {
      std::cout << initial_particle.post.mu << std::endl;
      std::cout << initial_particle.post.sigma << std::endl; 
      return 0;
    }
    return 0;
}









