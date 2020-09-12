#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "funhelpers.h"
#include "CharMatchIterator.h"
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

FILE * openTextFile(const char * fname, const bool binMode)
{
  FILE *f;
  try
  {
    if (binMode)
      f = fopen(fname, "rb");
    else
      f = fopen(fname, "r");
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

FILE * openTextFile(const char * fname)
{
  return openTextFile(fname, false);
}

bool setReadPos(FILE * f, long pos)
{
  fseek(f, 0, SEEK_END);
  long fileSize = ftell(f);
  if (pos < fileSize)
  {
    fseek(f, pos, 0);
    return true;
  }
  return false;
}

// [[Rcpp::export(track.matching.lines)]]
NumericVector trackMatchingLines(const Rcpp::String& filename,
                                CharacterVector& matchLines,
                                bool prefixOnly = false,
                                bool lineFeedOnly = true,
                                long binOffset = 0,
                                long maxBinRead=-1)
{
  NumericVector chinds;
  if (matchLines.isNULL() || matchLines.size() == 0) return chinds;
  CharacterVector::iterator it = matchLines.begin();
  FILE * f = NULL;
  try
  {
    f = openTextFile(filename.get_cstring(), true);
    static char buffer[4096];
    int size, lineCounter = 0, cind = 0;
    bool skipthisline = false, reading = true, prefixMatch = false, lastWasCr = false;
    const bool limitBinRead = maxBinRead >= 0;
    Rcpp::String chStr = *it;
    const char * chPtr = chStr.get_cstring();
    int chSize = str_size(chStr);
    long readCount = 0;
    if (binOffset > 0)
    {
      reading = setReadPos(f, binOffset);
    }
    while (reading && ((size = fread(buffer, 1, sizeof(buffer), f)) != 0))
    {
      for (int i=0; i<size; i++)
      {
        if (limitBinRead)
        {
          if (readCount > maxBinRead)
            break;
          readCount++; 
        }
        if (!reading) break;
        switch (buffer[i])
        {
        case '\r':
        case '\n':
          if (buffer[i] == '\r')
          {
            if (lineFeedOnly)
            {
              skipthisline = true;
              continue;
            }
            lastWasCr = true;
          }
          else if (lastWasCr)
          {
            lastWasCr = false;
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
          lastWasCr = false;
          continue;
        default:
          lastWasCr = false;
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


// [[Rcpp::export(track.matching.string.binpos)]]
NumericVector trackMatchingStringBinPos(const Rcpp::String& filename,
                                     StringVector& candidateStrings,
                                     long binOffset = 0,
                                     long maxBinRead = -1,
                                     bool lineStartOnly = true)
{
  NumericVector chind = NumericVector::create(NA_REAL);
  if (candidateStrings.isNULL() || candidateStrings.size() == 0) return chind;
  CharMatchIterator mit;
  mit.AppendStringVector(candidateStrings);
  FILE * f = NULL;
  try
  {
    f = openTextFile(filename.get_cstring(), true);
    static char buffer[4096];
    bool reading = true, skipthisline = false;
    int readSize;
    if (binOffset > 0)
    {
      reading = setReadPos(f, binOffset);
    }
    else binOffset = 0;
    size_t maxBinRead_size = 0;
    const bool limitBinRead = maxBinRead >= 0;
    if (limitBinRead)
      maxBinRead_size = static_cast<size_t>(maxBinRead);
    while (reading && ((readSize = fread(buffer, 1, sizeof(buffer), f)) != 0))
    {
      for (int i=0; i<readSize; i++)
      {
        if (limitBinRead && mit.totalCounter() >= maxBinRead_size)
        {
          reading = false;
        }
        if (!reading) break;
        char ch = buffer[i];
        switch (ch)
        {
        case '\r':
        case '\n':
          skipthisline = false;
          break;
        default:
          if (skipthisline) break;
          if (mit.Next(ch))
          {
            if (mit.IsFullyMatching())
            {
              chind[0] = mit.startMatchIndex() + binOffset;
              reading = false;
            }
          }
          else
          {
            mit.Reset();
            if (lineStartOnly)
              skipthisline = true;
          }
          continue;
        }
        mit.incrementTotalCounter();
      }
    }
    fcloseifopen(f);
  }
  catch(std::exception &ex)
  {
    fcloseifopen(f);
    forward_exception_to_r(ex);
  }
  return chind;
}

// [[Rcpp::export(track.next.line.binpos)]]
NumericVector trackNextLineBinPos(const Rcpp::String& filename,
                                    long binOffset = 0,
                                    int skipLines = 0)
{
  NumericVector chind = NumericVector::create(NA_REAL);
  FILE * f = NULL;
  try
  {
    f = openTextFile(filename.get_cstring(), true);
    static char buffer[4096];
    bool reading = true;
    int readSize;
    if (binOffset > 0)
    {
      reading = setReadPos(f, binOffset);
    }
    else binOffset = 0;
    long counter = binOffset;
    bool foundLine = false;
    if (skipLines < 0) skipLines = 0;
    while (reading && ((readSize = fread(buffer, 1, sizeof(buffer), f)) != 0))
    {
      for (int i=0; i<readSize; i++)
      {
        if (!reading) break;
        switch (buffer[i])
        {
        case '\r':
        case '\n':
          foundLine = true;
          break;
        default:
          if (foundLine)
          {
            if (skipLines == 0)
            {
              chind[0] = counter;
              reading = false;
              continue;
            }
            skipLines--;
            foundLine = false;
          }
          break;
        }
        counter++;
      }
    }
    fcloseifopen(f);
  }
  catch(std::exception &ex)
  {
    fcloseifopen(f);
    forward_exception_to_r(ex);
  }
  return chind;
  
}



// [[Rcpp::export(read.bin2buffer)]]
bool readBin2Mem(const Rcpp::String& filename, int skip, int binLen, intptr_t objPtr)
{
  FILE * f = NULL;
  try
  {
    f = openTextFile(filename.get_cstring(), true);
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

// [[Rcpp::export(clean.opened.streams)]]
bool cleanOpenedStreams()
{
  _fcloseall();
  return true;
}
