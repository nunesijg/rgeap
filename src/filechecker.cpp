#include <Rcpp.h>
#include <fstream>
#include "fileio.h"
using namespace Rcpp;
using namespace std;

//
//  File I/O Checker
//
//  Verifies the associated structure of input data files
//

// [[Rcpp::export(check.file.structure)]]
DataFrame checkFileStructure(const Rcpp::String& filename)
{
  NumericVector chinds;
  IntegerVector rows;
  IntegerVector cols;
  double chind = 0, lastlchind = 0, counter = 0;
  FILE *f = NULL;
  try
  {
    f = openTextFile(filename.get_cstring());
    if (f)
    {
      static char buffer[4096];
      int size, c = 1, r = 0,lastcols = 0;
      bool firstLine = true;
      while ((size = fread(buffer, 1, sizeof(buffer), f)) != 0)
      {
        for (int i=0; i<size; i++, counter++)
        {
          switch (buffer[i])
          {
          case '\t': ++c; continue;
          case '\n':
            if (c != lastcols)
            {
              if (firstLine)
                firstLine = false;
              else
              {
                chinds.push_back(chind);
                rows.push_back(r);
                cols.push_back(lastcols);
                chind = lastlchind;
              }
              lastcols = c;
              r = 0;
            }
            ++r;
            lastlchind = counter + 1;
            c = 1;
            continue;
          }
        }
      }
      if (counter > chind)
      {
        chinds.push_back(chind);
        rows.push_back(r);
        cols.push_back(lastcols);
      }
    }
    fcloseifopen(f);
  }
  catch(std::exception &ex)
  {
    fcloseifopen(f);
    forward_exception_to_r(ex);
  }
  catch(...)
  {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return DataFrame::create(_["charind"] = chinds, _["rows"] = rows, _["cols"]=cols);
}


// class FileDataStructure
// {
// public:
//   FileDataStructure(const Rcpp::String& filename)
//   {
//     _filename = filename;
//   }
// 
//   DataFrame checkFileStructure()
//   {
//     NumericVector v = {1,2};
//     return DataFrame::create(Named("v1") = v);
//   }
// 
// private:
//   Rcpp::String _filename;
//   //double _filename;
// };
// 
// RCPP_MODULE(mod_foo)
// {
//   class_<FileDataStructure>( "FileDataStructure" )
//   .constructor<Rcpp::String>()
//   .method( "checkFileStructure", &FileDataStructure::checkFileStructure )
//   ;
// }
