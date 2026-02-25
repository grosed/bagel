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
  
  prior_function_type prior_function = prior_example_2;
  feature_vector_function_type feature_vector_function = feature_vector_example_2;
  probability_type p = 1.0;
  probability_type p0 = 0.9;
  real_type s = 1.0;


  bagel_type bagel(prior_function,feature_vector_function,p0,p,s);
  
  try
    {
    while(cin)
      {
	// get the next value from the data stream
        getline(cin, input_line);
	double y = std::stod(input_line);
	bagel = update(bagel,y);
	std::cout << weight_0_t(bagel) << std::endl;
	
      };
        
      }
  catch(...)
    {
      std::cout << "end of data stream / invalid input" << std::endl;
      return 0;
    }
    return 0;
}









