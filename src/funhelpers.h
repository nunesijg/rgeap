#ifndef FUNHELPERS_H
#define FUNHELPERS_H

#include <Rcpp.h>
using namespace Rcpp;

bool str_startswith (std::string const &fullString, std::string const &ending);
bool str_endswith (std::string const &fullString, std::string const &ending);
int str_size(const Rcpp::String& str);

#endif /* FUNHELPERS_H */
