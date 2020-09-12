#include <Rcpp.h>
using namespace Rcpp;

//
//  Serialization functions
//
//  These functions are used to transfer data between C# and R in GEAP
//

static intptr_t staticRawPointer = (intptr_t)0;

// [[Rcpp::export(internal.ptr)]]
intptr_t internalPtr(RObject& obj)
{
  switch(obj.sexp_type())
  {
  case CLOSXP:
  case BUILTINSXP:
  case SPECIALSXP:
    stop("argument must not be a function");
    return 0;
  case NILSXP:
    stop("cannot trace NULL");
    return 0;
  case ENVSXP:
  case PROMSXP:
    stop("'internal.ptr' is not useful for promise and environment objects");
    return 0;
  case EXTPTRSXP:
  case WEAKREFSXP:
    stop("'internal.ptr' is not useful for weak reference or external pointer objects");
    return 0;
  case STRSXP:
    return (intptr_t)(void*)CHAR(STRING_ELT(obj, 0));
  }
  return (intptr_t)(DATAPTR_OR_NULL((SEXP)obj));
}

// Rcpp::String internalType(RObject& obj)
// {
//   static const std::string CMPPREF = "CMP:";
//   const int sxptype = obj.sexp_type();
//   if (sxptype == INTSXP && Rf_isFactor(obj))
//   {
//     return(CMPPREF + "factor");
//   }
//   switch(sxptype)
//   {
//   case INTSXP:
//   case REALSXP:
//   case STRSXP:
//     
//     break;
//   }
//   
//   return type2name(obj);
// }

// [[Rcpp::export(.get.mem.address)]]
intptr_t getMemAddress(SEXP& obj, bool outMemInfo = false)
{
  intptr_t retAddr;
  R_xlen_t size = 0;
  switch(TYPEOF(obj))
  {
  case NILSXP:
  {
    retAddr = 0;
    break;
  }
  case INTSXP:
  case REALSXP:
  case RAWSXP:
  case LGLSXP:
  {
    size = XLENGTH(obj);
    if (size == 0)
    {
      retAddr = 0;
    }
    else
    {
      retAddr = (intptr_t)(DATAPTR(obj));
    }
    break;
  }
  case CHARSXP:
  {
    size = XLENGTH(obj);
    Rcout << "??";
    retAddr = (intptr_t)(CHAR(obj));
    break;
  }
  case STRSXP:
  {
    size = XLENGTH(obj);
    if (size == 0)
    {
      retAddr = 0;
    }
    else
    {
      size = XLENGTH(STRING_ELT(obj, 0));
      retAddr = (intptr_t)(CHAR(STRING_ELT(obj, 0)));
    }
    break;
  }
  default:
  {
    stop("Cannot use this object type");
    return 0;
  }
  }
  if (outMemInfo)
  {
    Rcerr << "%MEM%\t" << "T:" << type2name(obj) << "\tA:" << (uintptr_t)retAddr << "\tL:" << size << "\n";
  }
  return retAddr;
}

// [[Rcpp::export]]
intptr_t rawPtr(RawVector& raw)
{
  staticRawPointer = (intptr_t)DATAPTR(raw);
  return staticRawPointer;
}

// [[Rcpp::export]]
intptr_t fixedRawPtr()
{
  return (intptr_t)(void*)(&staticRawPointer);
}

// [[Rcpp::export]]
CharacterVector rawToCharVec(RawVector& raw, bool firstOnly = false)
{
  uint32_t currSize = 0, i = 0, j = 0;
  std::vector<uint32_t> sizes;
  for (RawVector::iterator it = raw.begin(); it != raw.end(); ++it)
  {
    if ((*it) == '\0')
    {
      sizes.push_back(currSize);
      currSize = 0;
      if (firstOnly) break;
    }
    else
    {
      currSize++;
    }
  }
  if (currSize != 0) sizes.push_back(currSize);
  char * chptr = (char*)(RAW(raw));
  CharacterVector vec = CharacterVector(sizes.size());
  for (std::vector<uint32_t>::iterator it = sizes.begin(); it != sizes.end(); ++it, ++i)
  {
    uint32_t size = *it;
    vec[i] = Rf_mkCharLen(&chptr[j], size);
    j += size + 1;
  }
  return vec;
}

// [[Rcpp::export]]
RawVector charVecToRaw(StringVector& strvec)
{
  size_t totlen = 0;
  std::vector<size_t> sizes;
  uint32_t j = 0;
  for (R_xlen_t i = 0; i < strvec.size(); i++)
  {
    size_t size = LENGTH(strvec(i));
    totlen += (size + 1);
    sizes.push_back(size);
  }
  RawVector rawvec = RawVector(totlen);
  char * chptr = (char*)(RAW(rawvec));
  for (R_xlen_t i = 0; i < strvec.size(); i++)
  {
    memcpy(&chptr[j], CHAR(STRING_ELT(strvec, i)), sizes[i]);
    j += (sizes[i] + 1);
  }
  return rawvec;
}

// [[Rcpp::export]]
bool hasMinByteSize(SEXP& obj, double minSize)
{
  switch(TYPEOF(obj))
  {
  case NILSXP:
  {
    return false;
  }
  case INTSXP:
  case LGLSXP:
  {
    return (double)(sizeof(int) * XLENGTH(obj)) >= minSize;
  }
  case REALSXP:
  {
    return (double)(sizeof(double) * XLENGTH(obj)) >= minSize;
  }
  case RAWSXP:
  {
    return (double)(sizeof(char) * XLENGTH(obj)) >= minSize;
  }
  case STRSXP:
  {
    int size = (int)XLENGTH(obj);
    double totSize = 0;
    for (int i = 0; i < size; i++)
    {
      totSize += (double)XLENGTH(STRING_ELT(obj, i));
      if (totSize >= minSize) return true;
    }
    return false;
  }
  default:
  {
    return false;
  }
  }
}

// [[Rcpp::export]]
bool printPtrBytes(intptr_t ptr, unsigned int byteCount = 8)
{
  if (ptr == 0) return false;
  char * bins = (char*)ptr;
  for (unsigned int i = 0; i < byteCount; i++)
  {
    char byte = bins[i];
    for (int j = 0; j < 8; j++)
    {
      Rcout << ((byte & (1 << j)) >> j);
    }
    Rcout << ' ';
  }
  return true;
}

// [[Rcpp::export(multifactor.rebuild.base)]]
List multiFactorRebuild(IntegerVector& indexes, IntegerVector& sizes, StringVector& levels)
{
  R_xlen_t len = sizes.length(), allIndsLen = indexes.length(), lvlsLen = levels.length();
  List retLs = List(len);
  int k = 0;
  for (int i = 0; i < len; i++)
  {
    int currSz = sizes[i];
    StringVector strVec = StringVector(currSz);
    for (int j = 0; j < currSz && k < allIndsLen; j++, k++)
    {
      int ind = indexes[k];
      if (ind >= lvlsLen) 
      {
        stop("Index out of range");
        return List::create();
      }
      strVec[j] = levels[ind];
    }
    retLs[i] = strVec;
  }
  return retLs;
}

