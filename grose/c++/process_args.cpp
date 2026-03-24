

#include "process_args.h" 


std::tuple<double,double,double,int> process_args(int argc, char** argv)
{
      double p0 = 0, p = 0, s = 0;
      int n = 0;
      bool hasp0 = false, hasp = false, hass = false, hasn = false;

      
      for (int i = 1; i < argc; ++i)
	{
	  std::string arg = argv[i];
	  
	  if (arg == "--p0" && i + 1 < argc)
	    {
	      p0 = std::stod(argv[++i]);
	      hasp0 = true;
	    }
	  else if (arg == "--p" && i + 1 < argc)
	    {
	      p = std::stod(argv[++i]);
	      hasp = true;
	    }
	  else if (arg == "--s" && i + 1 < argc)
	    {
	      s = std::stod(argv[++i]);
	      hass = true;
	    }
	  else if (arg == "--n" && i + 1 < argc)
	    {
	      n = std::stoi(argv[++i]);
	      hasn = true;
	    }
	  else
	    {
	      throw std::runtime_error("Unknown or incomplete argument: " + arg);
	    }
	}
      if (!hasp0 || !hasp || !hass || !hasn)
	{
	  throw std::runtime_error("Usage: " + std::string(argv[0]) + " --p0 <value> --p <value> --s <value> --n <value>");
	}				 
      
      return std::tuple<double,double,double,int>(p0,p,s,n);
}
