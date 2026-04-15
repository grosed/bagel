#include <iostream>
#include "example_2.h"
#include "bagel.h"
#include <iostream>
#include <string>
#include <list>

#include "process_args.h"

using namespace std;

std::string input_line;

int main(int argc, char* argv[])
{

  std::tuple<double,double,double,int,double> command_line_args;
  try
    {
      command_line_args = process_args(argc,argv);
    }
  catch(const std::exception& e)
    {
      std::cerr << e.what() << std::endl;
    }
  
  prior_function_type prior_function = prior_example_2;
  feature_vector_function_type feature_vector_function = feature_vector_example_2;

  //probability_type p = 1.0;
  //probability_type p0 = 0.9;
  //real_type s = 1.0;


  
  probability_type p0 = std::get<0>(command_line_args);
  probability_type p = std::get<1>(command_line_args); 
  real_type s = std::get<2>(command_line_args);
  int n = std::get<3>(command_line_args);
  real_type t = std::get<4>(command_line_args);

  // dummy
  real_type lst_nu = 1.0;
  real_type lst_mu =0.0;

  
  bagel_type bagel(prior_function,feature_vector_function,p0,p,s,lst_nu,lst_mu,n);
  
  try
    {
    while(cin)
      {
	// get the next value from the data stream
        getline(cin, input_line);
	double y = std::stod(input_line);
	bagel = update(bagel,y);
	std::cout << weight_0_t(bagel) << std::endl;
	if(1.0 - weight_0_t(bagel) > t)
	  {
	    std::cout << std::endl;
	    auto particle_ratios = bagel.ratios();
	    for(auto& rs : particle_ratios)
	      {		
		for(auto& r : rs)
		  {
		    std::cout << r << " ";
		  }
		std::cout << std::endl;
	      }
	    break;
	  }
	
      };
        
      }
  catch(...)
    {
      std::cout << "end of data stream / invalid input" << std::endl;
      return 0;
    }
    return 0;
}









