
#include <iostream>


#include "examples.h"

#include "particle_type.h"
#include "bagel.h"
#include "process_args.h"
#include "noise_type.h"
#include "KL_divergence_type.h"

#include "model_type.h"



int main(int argc, char* argv[])
{

  // read command line arguments
  std::tuple<double,double,double,int,double> command_line_args;
  try
    {
      command_line_args = process_args(argc,argv);
    }
  catch(const std::exception& e)
    {
      std::cerr << e.what() << std::endl;
    }

  
  probability_type p0 = std::get<0>(command_line_args);
  probability_type p = std::get<1>(command_line_args); 
  real_type s = std::get<2>(command_line_args);
  int n = std::get<3>(command_line_args);
  real_type t = std::get<4>(command_line_args);

  
  // set up noise models
  known_variance kv;
  kv.sigma = 1.0;
  
  unknown_variance uv;
  uv.nu = 1.0;
  uv.iota = 1.0;

  // set up model
  model_type model;

  /*
  model.transformer_function =  transformation_example_1;
  model.prior_function =  prior_example_1;
  model.feature_vector_function =  feature_vector_example_1;
  */
  
  model.transformer_function =  transformation_daily;
  model.prior_function =  prior_daily;
  model.feature_vector_function =  feature_vector_daily;
  
  
  bagel_type<known_variance,KL_divergence_type::approximate> bagel(model,
								     kv, //kv,
								     p0,p,n);

  // read data from stdin
  std::string input_line;
  try
    {
      int i = 0;
      while(std::cin)
      {
	// get the next value from the data stream
        getline(std::cin, input_line);
	double y = std::stod(input_line);
	bagel = update(bagel,y);
	std::cout << weight_0_t(bagel) << std::endl;
	// test for threshold
	if(1.0 - weight_0_t(bagel) > t)
	  {
	    // report prediction
	    std::cout << std::endl;
	    std::cout << "----------------------------" << std::endl;
	    std::cout << "change detected at t = " << bagel.t << std::endl;
	    std::cout << "----------------------------" << std::endl;
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
