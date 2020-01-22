#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdarg.h>
#include "funhelpers.h"
using namespace Rcpp;
using namespace std;

//
//  File I/O Reader
//
//  Reads an already scanned file (faster than read.table)
//

bool fcloseifopen(FILE *f)
{
  if (f != NULL)
  {
    fclose(f);
    return true;
  }
  return false;
}

FILE * openTextFile(const char * fname)
{
  FILE *f;
  try
  {
    f = fopen(fname, "r");
    //errno_t err = fopen_s(&f, fname, "r");
    //if (err != 0)
    //{
    //  std::stringstream errstrs;
    //  errstrs << "Failed: could not open such file. Returned ";
    //  errstrs << err;
    //  throw std::runtime_error(errstrs.str());
    //}
    //else if (!f)
    if (f == NULL)
    {
      throw std::runtime_error(std::string("Failed: could not open such file"));
    }
    return f;
  }
  catch(std::exception &ex)
  {
    fcloseifopen(f);
    forward_exception_to_r(ex);
  }
  return NULL;
}

FILE * openBinFile(const char * fname)
{
  FILE *f;
  try
  {
    f = fopen(fname, "rb");
    if (f == NULL)
    {
      throw std::runtime_error(std::string("Failed: could not open such file"));
    }
    return f;
  }
  catch(std::exception &ex)
  {
    fcloseifopen(f);
    forward_exception_to_r(ex);
  }
  return NULL;
}

// [[Rcpp::export(track.matching.lines)]]
NumericVector trackMatchingLines(const Rcpp::String& filename,
                                CharacterVector& matchLines,
                                bool prefixOnly = false,
                                bool lineFeedOnly = true)
{
  NumericVector chinds;
  if (matchLines.isNULL() || matchLines.size() == 0) return chinds;
  CharacterVector::iterator it = matchLines.begin();
  FILE * f = NULL;
  try
  {
    f = openTextFile(filename.get_cstring());
    static char buffer[4096];
    int size, lineCounter = 0, cind = 0;
    bool skipthisline = false, reading = true, prefixMatch = false;
    Rcpp::String chStr = *it;
    const char * chPtr = chStr.get_cstring();
    int chSize = str_size(chStr);
    while (reading && ((size = fread(buffer, 1, sizeof(buffer), f)) != 0))
    {
      for (int i=0; i<size; i++)
      {
        switch (buffer[i])
        {
        case '\r':
        case '\n':
          if (lineFeedOnly && buffer[i] == '\r')
          {
            skipthisline = true;
            continue;
          }
          if (prefixMatch)
          {
            chinds.push_back(lineCounter, chStr);
            ++it;
            if (it == matchLines.end())
            {
              i = size;
              reading = false;
              continue;
            }
            else
            {
              chStr = *it;
              chPtr = chStr.get_cstring();
              chSize = str_size(chStr);
            }
          }
          cind = 0;
          skipthisline = false;
          prefixMatch = false;
          ++lineCounter;
          continue;
        case '\t':
          skipthisline = true;
          continue;
        default:
          if (prefixMatch && !prefixOnly)
          {
            prefixMatch = false;
          }
          if (skipthisline) continue;
          if (buffer[i] == chPtr[cind])
          {
            if (chSize == ++cind)
            {
              prefixMatch = true;
              skipthisline = true;
            }
          }
          else skipthisline = true;
          continue;
        }
      }
    }
    fcloseifopen(f);
    if (reading)
    {
      while (it != matchLines.end())
      {
        chStr = *it;
        chinds.push_back(NA_REAL, chStr);
        ++it;
      }
    }
  }
  catch(std::exception &ex)
  {
    fcloseifopen(f);
    forward_exception_to_r(ex);
  }
  return chinds;
}

// [[Rcpp::export(read.bin2buffer)]]
bool readBin2Mem(const Rcpp::String& filename, int skip, int binLen, intptr_t objPtr)
{
  FILE * f = NULL;
  try
  {
    f = openBinFile(filename.get_cstring());
    if (skip != 0) fseek(f, skip, 0);
    fread((void*)objPtr, 1, binLen, f);
    fcloseifopen(f);
  }
  catch(std::exception &ex)
  {
    fcloseifopen(f);
    forward_exception_to_r(ex);
    return false;
  }
  return true;
}



//
// // [[Rcpp::export(read.table.offset)]]
//DataFrame readTableOffset(const Rcpp::String& filename,
//                          const long offsetChars, const long amountChars,
//                          const int nrow, const int ncol)
//{
//  
//}

