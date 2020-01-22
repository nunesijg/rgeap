
#' @include base.R
#' @include c-serialist.R

##########################
# Serialization
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

#' @title GEAP serial
#' 
#' @description Funtions to transfer data between GEAP executable and R. Standalone use is not advised and must be avoided.
#' 
#' @name serial
#' @rdname serial
#' 
NULL

#' @rdname serial
#' @export
setGeneric('serialprint', function(x, ...) standardGeneric('serialprint'))

#' @rdname serial
#' @export
setMethod('serialprint', 'character', function(x, ...)
{
  if (mayserialize(x, ...))
  {
    print(xserialize(x))
  } else if (length(x) != 0) {
    writeLines(paste0('"', gsub("\n", "\\\\n", gsub("\\\\", "\\\\\\\\", x)), '"'))
  }
  invisible(0)
})

#' @rdname serial
#' @export
setMethod('serialprint', 'integer', function(x, ...)
{
  if (mayserialize(x, ...))
  {
    print(xserialize(x))
  } else {
    writeLines(paste0(as.character(x), 'L'))
  }
  invisible(0)
})

#' @rdname serial
#' @export
setMethod('serialprint', 'numeric', function(x, ...)
{
  if (mayserialize(x, ...))
  {
    print(xserialize(x))
  } else {
    rets = as.character(x)
    rets[is.na(x)] = "NaN"
    writeLines(rets)
  }
  invisible(0)
})

#' @rdname serial
#' @export
setMethod('serialprint', 'factor', function(x, ...)
{
  if (mayserialize(x, ...))
  {
    print(xserialize(x))
  } else {
    writeLines(as.character(x))
  }
  invisible(0)
})

#' @rdname serial
#' @export
setMethod('serialprint', 'logical', function(x, ...)
{
  if (length(x) == 0)
  {
    writeLines(character(0))
  } else {
    if (mayserialize(x, ...))
    {
      print(xserialize(x))
    } else {
      writeLines(as.character(x))
    }
  }
  invisible(0)
})

setMethod('serialprint', 'serializable_list', function(x, ...)
{
  if (mayserialize(x, ...))
  {
    print(xserialize(x))
  } else {
    print(x)
  }
  invisible(0)
})

#' @rdname serial
#' @export
setMethod('serialprint', 'matrix', function(x, ...)
{
  if (mayserialize(x, ...))
  {
    print(xserialize(x))
  } else {
    writeLines(paste0(c(sprintf("%s.matrix", class(x[1,1])) , (if(is.null(colnames(x))) paste0('C', 1:ncol(x)) else colnames(x)) ), collapse = '\t'))
    lins = apply(x, 1, paste0, collapse='\t')
    rnms = (if (is.null(rownames(x))) paste0("R", 1:nrow(x)) else rownames(x))
    writeLines(paste(rnms, lins, sep='\t'))
  }
  
  invisible(0)
})

setMethod('serialprint', 'serialelem', function(x, ...)
{
  con = stderr()
  if (hasArg('con')) con = list(...)$con
  memtg = "<mem>"
  if ('serialfile' %in% names(attributes(x))) memtg = sprintf("<mem serialfile=\"%s\">", attr(x, 'serialfile'))
  outmsg = c(getOption('tag.memory', '!M!'), memtg, as.character(x), "</mem>")
  writeLines(text = outmsg, con = con)
  invisible(0)
})


#' @rdname serial
#' @export
setMethod('serialprint', 'NULL', function(x, ...)
{
  writeLines(character(0))
  invisible(0)
})

#' @rdname serial
#' @export
setMethod('serialprint', 'ANY', function(x, ...)
{
  print(x)
  invisible(0)
})

.isMaySerialize <- function(...) hasArg('may.serialize') && any(as.logical(list(...)$`may.serialize`), na.rm = T)

setGeneric('mayserialize', function(x, ...) standardGeneric('mayserialize'))

setMethod('mayserialize', 'serializable_atomic', function(x, ...)
{
  if (.isMaySerialize(...))
  {
    if (is.raw(x)) return(T)
    return(hasMinByteSize(x, 4096))
  }
  FALSE
})

setMethod('mayserialize', 'serializable_list', function(x, ...)
{
  .isMaySerialize(...)
})

setMethod('mayserialize', 'matrix', function(x, ...)
{
  .isMaySerialize(...) && hasMinByteSize(x, 2048)
})


setMethod('mayserialize', 'ANY', function(x, ...)
{
  return(FALSE)
})



# Adjusts memory allocation
.adjust.alloc.size <- function(size)
{
  if (size == 0) return(1L)
  if (size <= 4096) return(4096L)
  return(2L^as.integer(ceiling(log2(size))))
}

# Allocates memory to a RAW object, returning its memory address as numeric
#' @rdname serial
#' @export
alloc.call.pointer <- function(size=4096)
{
  size = .adjust.alloc.size(size)
  if (size != length(.command.buffer))
  {
    .parent.assign(.command.buffer, raw(size))
  }
  return(rawPtr(.command.buffer))
}

# Gets the fixed address of the pointer that points to RAW data
#' @rdname serial
#' @export
fixed.call.pointer <- function()
{
  return(fixedRawPtr())
}

# Callbacks the stores buffer. Not recommended to run outside GEAP executable
.buffer.callback <- function(penv)
{
  ops = options(warn = -1)
  on.exit(options(ops))
  retVal2 = NULL
  parseAndEval = function()
  {
    command = rawToCharVec(.command.buffer, firstOnly = T)
    commandExpr = parse(text=command)
    retVal3 = eval(commandExpr, penv)
    return(retVal3)
  }
  retVal2 = withCallingHandlers(parseAndEval(),
                                message = function(c) invokeRestart('muffleMessage'),
                                warning = function(w) invokeRestart('muffleWarning') )
  return(retVal2)
}

# Validates and runs a command encoded as base64 string
#' @rdname serial
#' @export
eval.base64 <- function(str64, finalmsg="!0!", tag.error=getOption("tag.error", "!E!"), envir=globalenv(), prints=F, ...)
{
  if (is.character(finalmsg)) on.exit(message(finalmsg))
  retval = tryCatch(eval(parse(text=rawToChar(openssl::base64_decode(str64))), envir = envir),
                    error=function(e) message(sprintf("%s\n%s", tag.error, as.character(e))))
  if (prints)
  {
    serialprint(retval, ...)
  }
  invisible(0)
  
}

# Method to run inside GEAP. Has no effect as a pure R function
#' @rdname serial
#' @export
callrun <- function(storeResult=F, prints=F, finalmsg="!0!", envir=NULL)
{
  retVal = NULL
  penv = if (is.environment(envir)) envir else parent.env()
  success = T
  tryCatch({
    capture.output( assign('retVal', value = .buffer.callback(penv)), type = 'output')
  }, error = function(e) { assign('success', value=F ); message(sprintf("%s%s", getOption('tag.error', '!E!'), e)) })
  if (storeResult)
  {
    .parent.assign(.command.last.return, retVal)
  }
  if (prints)
  {
    serialprint(retVal)
  }
  if (is.character(finalmsg)) message(finalmsg)
  invisible(success)
}

# Method to run inside GEAP and return a value. Has no effect as a pure R function
#' @rdname serial
#' @export
call2value <- function()
{
  retVal = NULL
  penv = parent.frame()
  tryCatch({
    capture.output( assign('retVal', value = .buffer.callback(penv)), type = 'output')
  }, error = function(e) message(sprintf("%s%s", getOption('tag.error', '!E!'), e)) )
  return(retVal)
}

# Returns a numeric with the object pointer
# [[geapexport IntPtr GetMemoryPtr(call obj)]]
#' @rdname serial
#' @export
getMemPtr <- function(obj)
{
  if (is(obj, 'serialelem.atomic'))
  {
    return(.get.mem.address(obj@value, F))
  }
  .get.mem.address(obj, F)
}

# Prints memory information
#' @rdname serial
printMemInfo <- function(obj, autoserialize=T)
{
  if (!autoserialize && !is(obj, 'serializable_atomic')) stop("Cannot return non-atomic memory without autoserialize option")
  if (autoserialize)
  {
    if (!is(obj, 'serializable_t')) stop("Only objects of type 'serializable_t' can be serialized")
    lastSer = xserialize(obj, TRUE)
    show(lastSer)
  } else {
    if (inherits(obj, 'serialelem'))
    {
      show(obj)
    } else {
      .get.mem.address(obj, T)
    }
  }
  invisible(0)
}

# Serializes a object, returns a serialelem object
# [[geapexport void XSerialize(call value, bool assignReturn = true)]]
#' @rdname serial
#' @export
xserialize <- function(value, assign.return = T)
{
  newobj = NA
  if (is(value, 'serializable_list'))
  {
    newobj = new('serialelem.list', value)
  } else if (is(value, 'serializable_atomic'))
  {
    newobj = new('serialelem.atomic', value)
  } else if (is(value, 'serialelem'))
  {
    newobj = value
  } else stop("'value' must be a serializable_t object")
  if (assign.return)
  {
    .parent.assign(.last.serialized, newobj)
  }
  if (get.geap.option("serialwritemode", "mem") == "file")
  {
    fpath = base::tempfile()
    (con = .write.filecon(fpath, 'wb')) %using% {
      attr(newobj, 'serialfile') = fpath
      write.serialelem(newobj, con)
    }
  }
  newobj
}

# Writes the serialelement inside a file
write.serialelem <- function(selem, con)
{
  if (inherits(selem, "serialelem.list"))
  {
    for (nm in names(selem))
    {
      write.serialelem(selem@value[[nm]], con)
    }
  } else if (inherits(selem, 'serialelem.atomic'))
  {
    selemval = selem@value
    if (!is.vector(selemval))
    {
      selemval = as.vector(selemval)
    }
    base::writeBin(selemval, con = con, useBytes = T)
  }
  invisible(0)
}

# Unserializes a serialized and returns its content
# [[geapgeneric void XUnserialize(call obj)]]
#' @rdname serial
#' @export
setGeneric('xunserialize', function(x) standardGeneric('xunserialize'))

#' @rdname serial
#' @export
setMethod('xunserialize', 'serialelem.list', function(x)
{
  x = switch (x@type,
    matrix = {
      retv = matrix(xunserialize(x@value$values), ncol=as.integer(x@attrs['ncol']), nrow=as.integer(x@attrs['nrow']))
      if ('col.names' %in% names(x@value)) colnames(retv) = x@value@`col.names`
      if ('row.names' %in% names(x@value)) rownames(retv) = x@value@`row.names`
      retv
    },
    data.frame = {
      retv = data.frame(row.names = xunserialize(x@value$`row.names`))
      for (nm in names(x@value))
      {
        if (nm != 'row.names')
        {
          retv[,nm] = xunserialize(x@value[[nm]])
        }
      }
      retv
    },
    list = lapply(x@value, xunserialize),
    factor = {
      lvls = xunserialize(x@value$levels)
      factor(x = lvls[xunserialize(x@value$indexes)], levels = lvls)
    }
  )
  return(x)
})

#' @rdname serial
#' @export
setMethod('xunserialize', 'serialelem.atomic', function(x)
{
  val = x@value
  if (x@type == 'character')
  {
    val = rawToCharVec(val)
  }
  return(val)
})

#' @rdname serial
#' @export
setMethod('xunserialize', 'ANY', function(x)
{
  return(x)
})

# Allocates bytes to store characters
# [[geapexport void AllocStrVector(long binsize, int vecLen, bool assignReturn = true)]]
#' @rdname serial
#' @export
xalloc.character <- function(binsize, veclength=1, assign.return=TRUE)
{
  bins = raw(binsize)
  attr(bins, 'count') = veclength
  attr(bins, 'type') = 'character'
  swmode = get.geap.option("serialwritemode", "mem")
  set.geap.option("serialwritemode", "mem")
  on.exit({set.geap.option('serialwritemode', swmode)})
  retval = xserialize(bins, assign.return)
  retval
}

# Allocates bytes to store characters for a matrix
# [[geapexport void AllocStrMatrix(long binSize, int ncol, int nrow, long binSizeColNames, long binSizeRowNames, bool assignReturn = true)]]
#' @rdname serial
#' @export
xalloc.character.matrix <- function(binsize, ncol, nrow, binsize.colnames=0, binsize.rownames=0, assign.return=TRUE)
{
  inpls = list()
  inpls$values = xalloc.character(binsize, ncol * nrow, F)
  if (binsize.colnames > 0)
  {
    inpls$col.names = xalloc.character(binsize.colnames, ncol, F)
  }
  if (binsize.rownames > 0)
  {
    inpls$row.names = xalloc.character(binsize.rownames, nrow, F)
  }
  attr(inpls, 'matrixtype') = 'character'
  attr(inpls, 'ncol') = ncol
  attr(inpls, 'nrow') = nrow
  retval = xserialize(inpls, assign.return)
  retval
}

# Allocates bytes to store characters for a matrix
# [[geapexport void AllocMatrix(call defval, int ncol, int nrow, long binSizeColNames, long binSizeRowNames, bool assignReturn = true)]]
#' @rdname serial
#' @export
xalloc.matrix <- function(defval, ncol, nrow, binsize.colnames=0, binsize.rownames=0, assign.return=TRUE)
{
  if (!is(defval, 'serializable_atomic')) stop("'defval' must be a serializable_atomic type")
  inpls = list()
  inpls$values = xserialize(rep(defval, ncol * nrow), F)
  if (binsize.colnames > 0)
  {
    inpls$col.names = xalloc.character(binsize.colnames, ncol, F)
  }
  if (binsize.rownames > 0)
  {
    inpls$row.names = xalloc.character(binsize.rownames, nrow, F)
  }
  attr(inpls, 'matrixtype') = class(defval)
  attr(inpls, 'ncol') = ncol
  attr(inpls, 'nrow') = nrow
  
  retval = xserialize(inpls, assign.return)
  retval
}


# Reads the bytes from a file to the target object. CAUTION: This may crash the entire R session
# [[geapexport void ReadBinFile2Obj(string fileName, int skip, int binLen, long objPtr)]]
#' @rdname serial
#' @export
read.bin.file2obj <- function(fname, skip, binLen, objPtr)
{
  if (!file.exists(fname)) stop(sprintf("file does not exist: '%s'", fname))
  read.bin2buffer(fname, skip, binLen, objPtr)
  #if (!(is(obj, 'serializable_atomic') || inherits(obj, 'serialelem.atomic') )) stop("'obj' must be a 'serializable_atomic' or 'serialelem.atomic'")
  #if (inherits(obj, 'serialelem.atomic'))
  #{
  #  read.bin2buffer(fname, skip, binLen, obj@value)
  #} else {
  #  read.bin2buffer(fname, skip, binLen, obj)
  #}
  invisible(0)
}

# Saves the objects (by name) from an environment to a rdata file
# [[geapexport void SaveSession(string filename, params string[] objnames)]]
#' @rdname serial
#' @export
save.session <- function(filename, objnames)
{
  penv = parent.frame()
  if (length(objnames) == 0)
  {
    objnames = ls(envir=penv)
  } else {
    objnames = intersect(objnames, ls(envir=penv))
  }
  save(list=objnames, file=filename, compress = 'bzip2', envir = penv)
  invisible(0)
}

# Loads the objects (by name) from a rdata file to the current environment 
# [[geapexport void LoadSession(string filename)]]
#' @rdname serial
#' @export
load.session <- function(filename)
{
  penv = parent.frame()
  load(file=filename, envir=penv, verbose=F)
  invisible(0)
}

# Calls an expression locally as if it was called in GEAP
# This function is used for debug purposes
metacall <- function(exprs, env=globalenv())
{
  etxt = paste0(deparse(substitute(exprs)), collapse = ' ')
  str64 = openssl::base64_encode(etxt)
  eval.base64(str64, envir = env, prints = T, finalmsg = NULL)
}
