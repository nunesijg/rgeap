
#' @include base.R

##########################
# INPUT/OUTPUT FUNCTIONS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0


# Matches a regex pattern.
# Return the matching group captures, or character(0) if none
.regex.matches <- function(pattern, x)
{
  mt = gregexpr(pattern, x, perl=T)[[1]]
  mres = attributes(mt)
  mt = as.integer(mt)
  ismatch = mt[1][1] != -1
  vres = character(0)
  if (ismatch)
  {
    if ('capture.start' %in% names(mres))
    {
      vres = substr.multi(x, mres$capture.start[1,], mres$capture.length[1,])
      if (any(nchar(mres$capture.names) != 0))
      {
        vres = setNames(vres, mres$capture.names)
      }
    } else {
      vres = substr.multi(x, mt, mres$match.length)
    }
  }
  return(vres)
}

# Disposes an object being used
dispose <- function(x) UseMethod('dispose')
dispose.default <- function(x) F

# Closes a connection, if valid and opened
dispose.connection <- function(x)
{
  if (x %in% getAllConnections() && isOpen(x))
  {
    close(x)
    return(T)
  }
  return(F)
}

# Uses a resource (left side), executes expression (right side) and disposes
# the resource after the function call, if needed
`%using%` <- function(x, expr)
{
  xpar = substitute(x)
  xexpr = substitute(expr)
  locenv = new.env()
  locobj = NULL
  locname = NA_character_
  if (class(xpar) %in% c('(', '{', 'call'))
  {
    subj = as.character(xpar)
    vrgm = grep("^[a-zA-Z\\.\\_][\\w\\.\\_]*?$", subj, perl=TRUE, value = T)[1]
    if (length(vrgm) == 0 || is.na(vrgm))
    {
      locobj = local(expr = x, envir = locenv)
      if (length(ls(locenv)) != 0)
      {
        locname = ls(locenv)[1]
      }
    } else {
      locname = vrgm[1]
      assign(locname, NULL, envir = locenv)
    }
  } else if (class(xpar) %in% c('name'))
  {
    locname = as.character(xpar)
    assign(locname, NULL, envir = locenv)
  }
  on.exit({
    if (!is.na(locname))
    {
      locobj = locenv[[locname]]
    }
    dispose(locobj)
  })
  res = local(expr = expr, envir = locenv)
  invisible(res)
}

# Opens a new file connection to read
.open.filecon <- function(filename)
{
  if (!file.exists(filename)) stop(sprintf("File %s not found", basename(filename)))
  con = NULL
  if (length(grep("\\.gz$", filename, perl=T)) != 0)
  {
    con = gzfile(filename, open = 'rt')
  } else {
    con = file(filename, 'r')
  }
  if (is.null(con)) stop(sprintf("Could not open file: %s", basename(filename)))
  con
}

# Opens a new file connection to write
.write.filecon <- function(filename, open.mode='w')
{
  con = file(filename, open.mode)
  con
}

# Gets the text of a specific line from a file
.read.line.offset <- function(filename, skip=0, skipNul=T)
{
  strl = NA_character_
  (con = .open.filecon(filename)) %using% {
    if (skip > 0)
    {
      readLines(con, n = skip, ok = T, warn=F, skipNul = skipNul)
    }
    strl = readLines(con, n = 1, ok = T, warn = F, skipNul = skipNul)
  }
  return(strl)
}

# Gets a character vector with the lines from a file by specific indexes
.read.lines.byindex <- function(filename, indexes=1, one.based=T, skipNul = T)
{
  indexes = as.integer(indexes)
  if (is.null(indexes) || is.na(indexes)) stop("'indexes' must be a integer vector with one or more elements")
  vstrl = rep(NA_character_, length(indexes))
  if (!one.based) indexes = indexes + 1L
  skips = diff(c(0L, indexes)) - 1L
  (con = .open.filecon(filename)) %using% {
    for (i in 1:length(skips))
    {
      skip = skips[i]
      if (skip > 0)
      {
        readLines(con, n = skip, ok = T, warn=F, skipNul = skipNul)
      }
      vstrl[i] = readLines(con, n = 1, ok = T, warn = F, skipNul = skipNul)
    }
  }
  return(vstrl)
}

# Gets the first text line from a file
.read.firstline <- function(filename)
{
  .read.line.offset(filename, skip=0)
}

# [[geapexport void CloseAllStreams()]]
# Closes all connections being used
#' @export
conn.close.all <- function()
{
  closeAllConnections()
  while(grDevices::dev.cur() > 1)
  {
    tryCatch(grDevices::dev.off(), error=function(e) invisible(e))
  }
  invisible(0)
}

# [[geapexport assign ReadRDS(path fileName)]]
# Loads a RDS data file to a object
#' @export
read.rds.file <- function(filename)
{
  obj = readRDS(file = filename)
  obj
}

# [[geapexport void SaveRDS(path fileName, call obj)]]
# Saves an object to a RDS data file
#' @export
save.rds.file <- function(filename, obj)
{
  saveRDS(object = obj, file = filename, compress = 'bzip2')
  invisible(T)
}

# [[geapexport assign ReadTableAtPos(path fileName, long bytePos=0, int nrows=-1)]]
# Reads a table from a specific binary position. The nrows parameter does not include the row for column names
#' @export
read.table.atpos <- function(filename, bytepos=0, nrows=-1)
{
  dt = NULL
  (con = .open.filecon(filename)) %using% {
    seek(con, bytepos)
    dt = read.table(con,
                    header = T,
                    sep = '\t',
                    row.names = 1,
                    nrows = nrows,
                    blank.lines.skip = F,
                    skipNul = F,
                    check.names = F)
    
  }
  colnames(dt) = make.unique(colnames(dt), sep = '_')
  dt
}

# [[geapexport void WriteMatrix(path fileName, call obj, dots optArgs)]]
# Saves a matrix or data.frame to a text file
#' @export
write.matrix <- function(filename, obj, ...)
{
  m = obj
  if (!is.matrix(m) && !is.data.frame(m))
    m = eset.exprs(m)
  lsargs = list(..., format='txt', formatfn=utils::write.table, x=m, file=filename, quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
  m = lsargs$x
  if (identical(lsargs$row.names, TRUE) && identical(lsargs$col.names, TRUE) &&
      length(rownames(m)) != 0L && length(colnames(m)) != 0L)
  {
    if (identical(rownames(m), m[,1L,drop=TRUE]))
      lsargs$row.names = FALSE
    else if (!(toupper(colnames(m)[1]) %in% c("ID", "ID_REF", "SPOT_ID")))
    {
      m = data.frame(ID=rownames(m), m, check.names = FALSE)
      lsargs$x = m
      lsargs$row.names = FALSE
    }
  }
  writefn = switch(lsargs$format,
                   txt = utils::write.table,
                   csv = utils::write.csv,
                   get(lsargs$formatfn)
                   )
  writeargnms = switch(lsargs$format,
                       csv = setdiff(formalArgs(write.table),
                                     c("append", "col.names", "sep",
                                       "dec", "qmethod")),
                       txt = formalArgs(write.table),
                       formatArgs(writefn))
  lsargs = lsargs[intersect(writeargnms, names(lsargs))]
  suppressWarnings(do.call(writefn, lsargs))
  invisible(0)
}
