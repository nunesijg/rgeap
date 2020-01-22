
#ifndef FILEIO_H
#define FILEIO_H

#include <Rcpp.h>
using namespace Rcpp;

bool fcloseifopen(FILE *f);
FILE * openTextFile(const char * fname);

#endif /* FILEIO_H */
