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
  
  try
    {
    while(cin)
      {
	// get the next value from the data stream
        getline(cin, input_line);
	double y = std::stod(input_line);
	if(t == 1)
	  {
	    particle_type initial_particle(prior_example_2,feature_vector_example_2,tau,p0,p,s);
	    initial_particle = theorem_2(initial_particle,t);
	    initial_particle = theorem_1(initial_particle,initial_particle,t);
	    initial_particle = theorem_4(initial_particle,t,y);
	    initial_particle = theorem_3(initial_particle,t,y);
	    initial_particle.weight = p0; 
	    particles.push_back(initial_particle);
	  }
	else
	  {
	    // new particle
	    particle_type initial_particle = particles.front();
	    std::cout << "here again " << initial_particle.weight << std::endl;
	    particles.push_back(theorem_2(initial_particle,t));
	    real_type sum_of_weights = 0.0;
	    for(auto& p : particles)
	      {
		p = theorem_1(initial_particle,p,t);
		p = theorem_4(p,t,y);
		p = theorem_3(p,t,y);
		std::cout << p.p0 << std::endl;
		std::cout << p.s << std::endl;
		sum_of_weights += p.weight;
	      }
	    std::cout << "here again " << initial_particle.weight << std::endl;
	    for(auto& p : particles)
	      {
	    	p.weight = p.weight/sum_of_weights;
	      }	    
	  }
	t = t + 1;
    };
    }
  catch(...)
    {
      particle_type& initial_particle = particles.front();
      std::cout << initial_particle.post.mu << std::endl;
      std::cout << initial_particle.post.sigma << std::endl;

      initial_particle = particles.back();
      std::cout << initial_particle.post.mu << std::endl;
      std::cout << initial_particle.post.sigma << std::endl;

      for(auto& p : particles)
	{
	  std::cout << p.weight << std::endl;
	}

      
      return 0;
    }
    return 0;
}









