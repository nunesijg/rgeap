#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#include <fstream>

//
//  Function helpers
//
//  Some function helpers to better performance
//


bool str_startswith (std::string const &fullString, std::string const &ending)
{
  return fullString.rfind(ending, 0) == 0;
}

bool str_endswith (std::string const &fullString, std::string const &ending)
{
  if (fullString.length() >= ending.length())
  {
    return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
  }
  else
  {
    return false;
  }
}

// [[Rcpp::export(substr.multi)]]
CharacterVector multiSubstring(std::string const &text,
                               IntegerVector& starts, IntegerVector& sizes)
{
  int count = starts.size();
  CharacterVector res(count);
  bool hasSizes = !(sizes.isNULL() || sizes.length() == 0);
  for (int i = 0; i < count; i++)
  {
    unsigned int start = starts[i] - 1; // R is one-based, C++ is zero-based
    if (start > text.length())
    {
      start = text.length();
    }
    unsigned int maxSize = text.length() - start;
    if (maxSize > text.length())
    {
      maxSize = text.length();
    }
    unsigned int size = hasSizes ? sizes[i % sizes.length()] : maxSize;
    if (size > maxSize)
    {
      size = maxSize;
    }
    char * substr = new char[size + 1];
    substr[size] = '\0';
    for (unsigned int j = 0; j < size; j++)
    {
      substr[j] = text[start + j];
    }
    res[i] = Rcpp::String(substr);
    delete [] substr;
  }
  return res;
}


// [[Rcpp::export(max.index)]]
uint64_t max_index(NumericVector& vals)
{
  if (vals.isNULL() || vals.size() == 0) return NA_REAL;
  NumericVector::iterator it = vals.begin();
  double bestval = *it;
  uint64_t counter = 1, bestind = 1;
  ++it;
  while (it != vals.end())
  {
    counter++;
    double currval = *it;
    if (currval > bestval)
    {
      bestval = currval;
      bestind = counter;
    }
    ++it;
  }
  return bestind;
}

// [[Rcpp::export(str.size)]]
int str_size(const Rcpp::String& str)
{
  return LENGTH(str.get_sexp());
}

// [[Rcpp::export(min.index)]]
uint64_t min_index(NumericVector& vals)
{
  if (vals.isNULL() || vals.size() == 0) return NA_REAL;
  NumericVector::iterator it = vals.begin();
  double bestval = *it;
  uint64_t counter = 1, bestind = 1;
  ++it;
  while (it != vals.end())
  {
    counter++;
    double currval = *it;
    if (currval < bestval)
    {
      bestval = currval;
      bestind = counter;
    }
    ++it;
  }
  return bestind;
}

