#ifndef ___PROCESS_ARGS_H___
#define ___PROCESS_ARGS_H___

#include <string>
#include <tuple>
#include <stdexcept>

std::tuple<double,double,double,int,double> process_args(int, char**);

#endif
