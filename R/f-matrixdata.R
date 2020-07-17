
#' @include base.R

##########################
# MATRIX MATH FUNCTIONS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

#' @name matrixdata
#' @title Matrix Data Helpers
#' 
#' @description Matrix functions for internal operations in GEAP.
#' They were meant to be exclusively used in GEAP but can also be applied in any other application.
#' 
#' @details
#' \code{column.classes} finds the classes correspondent to each column in a \code{data.frame}.
#' If a \code{matrix} is given, the same class is returned for all columns, since matrix values consist of a single class.
#' If a \code{list} is given, the classes will correspond to the list elements.
#' 
#' @return 
#' Named character vector with the column names (for \code{matrix} or \code{data.frame}) or
#' keys (for \code{list}) and their respective classes.
#' 
NULL


#' @rdname matrixdata
#' 
#' @export
#' 
column.classes <- function(m)
{
  if (is.matrix(m))
  {
    return(apply(m[NULL,,drop=F], 2, class))
  } else if (is.data.frame(m))
  {
    return(sapply(m[NULL,,drop=T], class))
  } else if (is.list(m))
  {
    return(sapply(m, function(x) class(x)[1]))
  }
  warning(sprintf("'m' is expected to be a matrix, list, or data.frame. %s given", class(m)[1]))
  logical(0)
}

.initialize.matrixdata <- function()
{
  loadpkgs('Biobase')
  .self.oneshot()
}

# Returns an integer vector with the index of the first occurrence of each value
.occur.index <- function(vc = NA_character_, sequential=T, keep.names=F)
{
  if ((length(vc) == 0) || all(is.na(vc))) return(NULL)
  if (!is.character(vc)) vc = as.character(vc)
  uvc = unique(vc)
  uvci = NULL
  if (sequential)
  {
    uvci = setNames(1L:length(uvc), uvc)
  } else {
    uvci = setNames((1:length(vc))[!duplicated(vc)], uvc)
  }
  vret = uvci[vc]
  if (keep.names) return(vret)
  return(unname(vret))
}

# Returns a list with the indexes of the first occurrence of each value
.occur.index.list <- function(vc = NA_character_)
{
  if ((length(vc) == 0) || all(is.na(vc))) return(NULL)
  if (!is.character(vc)) vc = as.character(vc)
  uvc = unique(vc)
  lsres = tapply(1:length(vc), as.factor(vc), list)[uvc]
  return(lsres)
}

# Returns a logical vector of the same length of vc, where TRUE values correspond to elements with duplicates
.occur.duplicates <- function(vc = NA_character_)
{
  visdups = setNames(duplicated(vc) | duplicated(vc, fromLast=T), vc)
  visdups = visdups[!duplicated(vc)]
  return(visdups)
}

# Groups the info by unique values from a specific column, and applies a function on elements from each group
.aggregate.dtinfo <- function(dtinfo, column=1, func, ...)
{
  if (is.null(dim(dtinfo))) stop("dtinfo parameter must be data.frame and matrix")
  if (ncol(dtinfo) <= 1) stop("dtinfo must have two or more columns")
  if (length(column) != 1 || is.na(column)) stop("column parameter must be a single atomic value")
  if (!is.function(func)) stop("func parameter must be a function")
  if (is.character(column))
  {
    cind = which(column == colnames(dtinfo))
    if (length(cind) == 0) stop(sprintf("No column named '%d'", column))
    column = cind[1]
  }
  infofs = as.character(dtinfo[,column])
  if (length(unique(infofs)) == length(infofs))
  {
    rownames(dtinfo) = as.character(dtinfo[,column])
    dtinfo = dtinfo[,-column, drop=F]
    return(dtinfo)
  }
  lsinfs = .occur.index.list(infofs)
  lsinfs = lsinfs[sapply(lsinfs, length) > 1]
  voccuns = !duplicated(infofs)
  dtinfowc = dtinfo[,-column,drop=F]
  dtinfotmp = dtinfo[voccuns,]
  rownames(dtinfotmp) = dtinfotmp[,column]
  dtinfotmp = dtinfotmp[,-column, drop=F]
  for (i in 1:length(lsinfs))
  {
    nm = names(lsinfs)[i]
    dtinfotmp[nm,] = apply(dtinfowc[lsinfs[[i]],,drop=F], 2, func, ...)
  }
  return(dtinfotmp)
}

# Groups the info by unique values from a specific column, and concatenates the group elements into a single data.frame
.aggregate.dtinfo.paste <- function(dtinfo, column=1, collapse='\n')
{
  .aggregate.dtinfo(dtinfo, column, paste0, collapse=collapse)
}


# Orders character vector based on order occurrence to another character vector
.order.byref <- function(x, refx, keep.names=F)
{
  inters = intersect(refx, x)
  if (length(inters) == 0) return(NULL)
  sorter = rep(0, length(x))
  selector = x %in% inters
  indexer = setNames(which(selector), x[selector])
  indexer = indexer[inters]
  indexer2 = which(!selector)
  if (keep.names)
  {
    return(c(indexer, setNames(indexer2, x[indexer2])))
  }
  return(c(unname(indexer), indexer2))
  
  #return(order(indexer))
  #ret = c(sort(x[indexer]) )
  #ret
}

# Sorts character vector based on order occurrence to another character vector
.sort.byref <- function(x, refx)
{
  return(x[.order.byref(x, refx)])
}

# [[geapexport bool ContainsZero(call vecOrMat)]]
#' @export
vecmat.has.zeroes <- function(v) any(v == 0)

# [[geapexport bool ContainsNA(call vecOrMat)]]
#' @export
vecmat.has.NA <- function(v) anyNA(v)

# [[geapexport bool ContainsInfinite(call vecOrMat)]]
#' @export
vecmat.has.infinite <- function(v) any(is.infinite(v))

# [[geapexport int RowCount(call matrix)]]
#' @export
mat.nrow <- function(m) nrow(m)

# [[geapexport string[] GetStrArray(call vecOrMat, int limit=0)]]
#' @export
vecmat.getstrarr <- function(v, limit=0)
{
  if (!is.character(v)) v = as.character(v)
  if (limit > 0)
  {
    return(head(v, n = limit))
  }
  return(v)
}

# [[geapexport string[] RowNames(call matrixOrESet, int limit=0)]]
#' @export
mat.rownames <- function(m, limit=0)
{
  .initialize.matrixdata()
  if(inherits(m, 'ExpressionSet'))
  {
    return(vecmat.getstrarr(rownames(exprs(m)), limit))
  }
  if(inherits(m, 'MAList'))
  {
    return(vecmat.getstrarr(rownames(m$M), limit))
  }
  vecmat.getstrarr(rownames(m), limit)
}

# [[geapexport string[] ColumnNames(call matrixOrESet, int limit=0)]]
#' @export
mat.colnames <- function(m, limit=0)
{
  if(inherits(m, 'ExpressionSet'))
  {
    return(vecmat.getstrarr(colnames(exprs(m)), limit))
  }
  if(inherits(m, 'MAList'))
  {
    return(vecmat.getstrarr(colnames(m$M), limit))
  }
  vecmat.getstrarr(colnames(m), limit)
}


#mat.merge <- function(...)
#{
#  .initialize.matrixdata()
#  argls = filter.args.bytype(types = 'matrix', ...)
#  if (length(argls) == 0) stop("At least one argument must be a matrix")
#  
#  
#}
