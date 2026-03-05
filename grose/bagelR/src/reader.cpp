
#include <iostream>
#include <string>

using namespace std;

std::string input_line;

int main()
{
  try
    {
    while(cin)
      {
        getline(cin, input_line);
	double val = std::stod(input_line);
	std::cout << val << std::endl;
    };
    }
  catch(...)
    {
      return 0;
    }
}
