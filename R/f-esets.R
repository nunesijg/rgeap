#' @include base.R
#' @include f-fileio.R
#' @include f-matrixdata.R

##########################
# EXPRESSION SET FUNCTIONS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

.initialize.esets <- function()
{
  loadpkgs('Biobase')
  .self.oneshot()
}

# Insert the platform feature data to a ExpressionSet
# [[geapexport assign ESetInsertPlatformData(call eset, call platDataFrame)]]
#' @export
eset.insert.platdata <- function(eset, platdf)
{
  if (ncol(platdf) == 0) return(eset)
  platrnms = rownames(platdf)
  # TODO: platdf teria que ter uma coluna identificadora (ID), assim poderia-se dar merge nos row.names entre eset e platdf
  dt = exprs(eset)
  match.count = function(a, b) sum(!is.na(base::match(a, b)))
  rnms = rownames(dt)
  matchinds = integer(0)
  if (is.null(platrnms) || match.count(rnms, platrnms) == 0)
  {
    platrnms = platdf[,1]
    if (ncol(platdf) > 1 && match.count(rnms, platrnms) != length(rnms))
    {
      mcols = sapply(colnames(platdf), function(j) match.count(rnms, platdf[,j]))
      if (all(mcols) == 0) return(eset)
      bestcol = names(mcols[max.index(mcols)])
      platrnms = platdf[,bestcol]
    }
  }
  matchinds = base::match(rnms, platrnms)
  if (sum(!is.na(matchinds)) == 0) return(eset)
  platdf = platdf[matchinds,,drop=F]
  rownames(platdf) = rnms
  vmd = names(platdf)
  featureData(eset) = new("AnnotatedDataFrame", data = platdf, varMetadata = data.frame(vmd, row.names = vmd))
  eset
}



# Gets the columns of data.frame whose values inherit the type names
df.getcols.bytype <- function(df, types)
{
  colsel = sapply(df, inherits, what=types)
  return(df[,colsel, drop=F])
}

# Converts a matrix or data.frame to ExpressionSet (exprs slot for numeric and featureData for character) 
# [[geapexport assign DataFrame2ESet(call dfOrMat)]]
#' @export
df.to.eset <- function(df)
{
  .initialize.esets()
  if (inherits(df, 'ExpressionSet'))
  {
    return(df)
  }
  if (inherits(df, 'MAList'))
  {
    return(malist.to.eset(df))
  }
  if (is.matrix(df))
  {
    return(new("ExpressionSet", exprs=df))
  }
  valCols = df.getcols.bytype(df, c('numeric', 'integer'))
  #if (ncol(valCols) == 0) stop("Not enough columns with numeric values")
  valCols = as.matrix(valCols)
  eset = new("ExpressionSet", exprs = valCols)
  
  attrCols = df.getcols.bytype(df, c('character', 'factor', 'multifactor', 'list'))
  hasPlat = ncol(attrCols) != 0
  attributes(eset)$hasPlatform = hasPlat
  if (hasPlat)
  {
    eset = eset.insert.platdata(eset, attrCols)
    #vmd = names(attrCols)
    #featureData(eset) = new("AnnotatedDataFrame", data = attrCols, varMetadata = data.frame(vmd, row.names = vmd))
  }
  return(eset)
}

# Reads a table as data.frame and converts it to ExpressionSet (exprs slot for numeric and featureData for character) 
# [[geapexec assign ReadTXT2ESet(string fileName, int skip=0)]]
#' @export
txtfile.to.eset <- function(filename, nskip=0)
{
  .initialize.esets()
  dt = read.table(filename, header = T, sep='\t', skip=nskip, row.names = 1, blank.lines.skip = F, skipNul = F, check.names = F)
  eset = df.to.eset(dt)
  return(dt)
}

# Merges two or more matrices, data.frames or ExpressionSet's into a single ExpressionSet
# [[geapexec assign MergeESets(params call[] objs)]]
#' @export
esets.merge <- function(...)
{
  .initialize.esets()
  argls = .get.valid.args(...)
  argls$types = c('matrix', 'data.frame', 'ExpressionSet')
  argls = base::do.call(filter.args.bytype, argls)
  #argls = filter.args.bytype(types = c('matrix', 'data.frame', 'ExpressionSet'), ...)
  if (length(argls) == 0) stop("At least one argument must be a matrix, data.frame or ExpressionSet")
  if (length(argls) == 1) return(df.to.eset(argls[[1]]))
  rownms = NULL
  dtls = list()
  platls = list()
  platcolnms = character(0)
  matinds = NULL
  hasargnms = !is.null(names(argls))
  ncolfinal = 0
  for (i in 1:length(argls))
  {
    ci = as.character(i)
    dt = argls[[i]]
    ceset = df.to.eset(dt)
    dt = Biobase::exprs(ceset)
    if (eset.hasplatform(ceset))
    {
      dtplat = eset.getplatform(ceset)
      pcnms = base::setdiff(colnames(dtplat), platcolnms)
      if (length(pcnms) != 0)
      {
        platls[[ci]] = dtplat[,pcnms,drop=F]
        platcolnms = c(platcolnms, pcnms)
      }
    }
    #if (inherits(dt, "MAList")) dt = malist.to.eset(dt)
    #if (inherits(dt, "ExpressionSet")) dt = exprs(dt)
    #if (inherits(dt, "data.frame")) dt = as.matrix(df.getcols.bytype(df, c('integer', 'numeric')))
    dtls[[i]] = dt
    if (is.null(rownms))
    {
      rownms = rownames(dt)
      matinds = matrix(0L, nrow=nrow(dt), ncol=length(argls), dimnames = list(rownms, character(0)) ) #data.frame(V1=(1L:nrow(dt)), row.names = rownms)
      matinds[,1] = 1L:nrow(dt)
      ncolfinal = ncol(dt)
    } else {
      minds = match(rownms, rownames(dt))
      #rownms = intersect(rownms, rownames(dt))
      if (all(minds == 0)) stop(sprintf("Argument %s has no matching row names", (if (hasargnms) sprintf("%s (%d)", names(argls)[i], i) else as.character(i) ) ))
      matinds[,i] = minds
      ncolfinal = ncolfinal + ncol(dt)
    }
  }
  rsel = !apply(matinds, 1, anyNA)
  if (!any(rsel)) stop("Total matrix has no common intersections")
  matinds = matinds[rsel,]
  fmat = matrix(NA_real_, nrow=nrow(matinds), ncol=ncolfinal, dimnames = list(rownames(matinds), character(0)) )
  fdfplat = data.frame(row.names = rownames(matinds))
  #platmat = data.frame()
  currcol = 1L
  assign('ainds', matinds, envir=globalenv())
  assign('afmat', fmat, envir=globalenv())
  for (i in 1:length(dtls))
  {
    ci = as.character(i)
    if (ci %in% names(platls))
    {
      dtplat = platls[[ci]]
      fdfplat[, colnames(dtplat)] = dtplat[matinds[,i],,drop=F]
    }
    dt = dtls[[i]]
    if (ncol(dt) == 0) next
    fcol = currcol + ncol(dt) - 1L
    fmat[,(currcol:fcol)] = dt[matinds[,i],,drop=F]
    currcol = fcol + 1L
  }
  colnames(fmat) = make.unique(unlist(lapply(dtls, colnames)))
  eset = df.to.eset(fmat)
  if (ncol(fdfplat) != 0) eset = eset.insert.platdata(eset, fdfplat)
  eset
}

# Removes NA entries from ExpressionSet
# [[geapexport assign ESetRemoveInvalidValues(call eset, bool rmNA=true, bool rmZero=true, bool rmNegative=true)]]
#' @export
eset.rm.invalid.values <- function(eset, rm.na=T, rm.zero=T, rm.neg=T)
{
  .initialize.esets()
  m = eset
  iseset = inherits(eset, "ExpressionSet")
  if (iseset)
  {
    m = exprs(eset)
  }
  rminds = apply(m, 1, anyNA)
  if (rm.zero) rminds = rminds | apply(m, 1, function(v) any(v == 0, na.rm = T))
  if (rm.neg) rminds = rminds | apply(m, 1, function(v) any(v < 0, na.rm = T))
  if (!any(rminds)) return(eset)
  m = m[!rminds,,drop=F]
  if (iseset)
  {
    neweset = df.to.eset(m)
    if (eset.hasplatform(eset))
    {
      platf = eset.getplatform(eset)[!rminds,,drop=F]
      neweset = eset.insert.platdata(platf)
    }
    return(neweset)
  }
  m
}

# Gets the intensity values matrix from ExpressionSet, AffyBatch, EListRaw, EList, etc
# [[geapexport assign SetExprsMatrix(call eset)]]
#' @export
eset.exprs <- function(eset)
{
  .initialize.esets()
  if(inherits(eset, c("ExpressionSet", "AffyBatch"))) return(Biobase::exprs(eset))
  if(inherits(eset, c("EListRaw", "EList"))) return(eset$E)
  if(inherits(eset, "MAList")) return(Biobase::exprs(malist.to.eset(eset)))
  if(inherits(eset, "RGList")) return(eset$R)
  eset
}

# Gets the sample names from a data (ExpressionSet, matrix, etc)
# [[geapexport string[] GetSampleNames(call eset)]]
#' @export
eset.samplenames <- function(eset)
{
  .initialize.esets()
  if(inherits(eset, "ExpressionSet")) return(Biobase::sampleNames(eset))
  mat = eset.exprs(eset)
  smpnms = colnames(mat)
  smpnms
}

# Checks if a ExpressionSet has platform information
# [[geapexport bool ESetHasPlatform(call eset)]]
#' @export
eset.hasplatform <- function(eset)
{
  .initialize.esets()
  hasplat = ncol(featureData(eset)) != 0
  hasplat
}

# Gets the platform data.frame from a ExpressionSet, EList, MAList or data.frame
# [[geapexport assign GetPlatformData(call eset)]]
#' @export
eset.getplatform <- function(eset)
{
  .initialize.esets()
  if(inherits(eset, c("ExpressionSet", "AffyBatch"))) return(fData(eset))
  if(inherits(eset, c("EListRaw", "EList", "RGList", "MAList"))) return(eset$genes)
  plat = df.getcols.bytype(eset, c('character', 'factor', 'multifactor', 'list'))
  plat
}

