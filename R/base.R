
#' @include rgeap-package.R
#' @include vardefs.R

##########################
# Essential functions
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

#' @title GEAP basic Run-Time
#' 
#' @description Basic functions to perform together with GEAP executable run-time.
#' 
#' @name geap-base
#' @rdname geapbase
#' 
NULL

# Use this in functions to be called once
.self.oneshot <- function()
{
  cl = sys.call(-1)
  f = get(as.character(cl[[1]]), mode="function", sys.frame(-2))
  cl = match.call(definition=f, call=cl)
  funcname = as.character(cl)
  penv = parent.env(parent.frame())
  lockedbind = bindingIsLocked(funcname, penv)
  if (lockedbind) unlockBinding(funcname, penv)
  assign(x = funcname, value = function(...){ invisible(0) }, pos = -1, envir = penv)
  if (lockedbind) lockBinding(funcname, penv)
  invisible(0)
}

# Assigns a value in the context of the parent function or environment
.parent.assign <- function(def, value, degree=1)
{
  defname = as.character(substitute(def))
  penv = parent.env(parent.frame(degree))
  if (bindingIsLocked(defname, penv)) unlockBinding(defname, penv)
  assign(x = defname, value = value, envir = penv)
  return(defname)
}

# Gets only the valid (existent) arguments inside ellipsis
.get.valid.args <- function(...)
{
  argls = list()
  for (i in 1:...length())
  {
    tryCatch({argls[[i]] = ...elt(i)}, error=function(x) invisible(0))
  }
  argls
}

# Example:
#loader.first <- function()
#{
#  .self.oneshot()
# ...  
#}

# All backslashes and multiple forward slashes become a single slash
# [[geapexport bool RemoveVars(params string[] varnames)]]
#' @rdname geapbase
#' @export
remove.vars <- function(varnms)
{
  exvars = sapply(varnms, exists)
  if (any(exvars))
  {
    exvars = names(exvars)[exvars]
    remove(list = exvars)
    return(T)
  }
  return(F)
}

# All backslashes and multiple forward slashes become a single slash
# [[geapexport string AdjustPath(string path)]]
#' @rdname geapbase
#' @export
adjust.path <- function(path) sub(pattern = '[\\\\/]+', replacement = '/', x = path, perl = T)

# Appends path to temporary directory
to.tempdir <- function(filename) adjust.path(file.path(temp.dir, filename))

# Defines a value for a local option
# [[geapexport void SetOption(string key, string value)]]
#' @export
set.geap.option <- function(key, value)
{
  argls = list()
  argls[[key]] = as.character(value)
  do.call(options, argls)
  invisible(T)
}

# Defines values for local options
# [[geapexport void SetOptions(params named_string[] args)]]
#' @export
set.geap.options <- function(...)
{
  argls = list(...)
  do.call(options, argls)
  invisible(T)
}

# Defines a value for a local option
# [[geapexport string GetOption(string key, string defValue="NULL")]]
#' @export
get.geap.option <- function(key, defVal)
{
  getOption(x = key, default = defVal)
}


# Prints the contents of a variable to be read by GEAP
#' @rdname geapbase
#' @export
print.result <- function(x)
{
  if(is.null(x)) writeLines('');
  if (is.vector(x)) writeLines(paste(x, sep = "\t", collapse = "\t"))
  else if (is.matrix(x) || is.data.frame(x))
  {
    writeLines(apply(x, 1, paste, sep="\t", collapse="\t"), sep="\n", useBytes = T)
  } else writeLines(x)
}

# Combines head() and tail() in a single function
head.tail <- function(mat, colname=NULL, maxRows = 10, decOrder = F)
{
  ord = order((if(is.null(colname)) 1:nrow(mat) else (if (colname == 'row.names') rownames(mat) else mat[,colname])), decreasing = decOrder)
  return(mat[unique(c(head(ord, n = ceiling(maxRows/2)),
                      tail(ord, n = floor(maxRows/2)))),,drop=F])
}

# Gets arguments from ellipsis of specified (or inherited from) types
filter.args.bytype <- function(types, ...)
{
  argls = list(...)
  argls = argls[sapply(argls, inherits, what=types)]
  argls
}

# Sends a status message to GEAP's GUI (C#)
.give.status <- function(percent = NA_integer_, message = NA_character_, engMsg = message)
{
  if (!is.na(percent))
  {
    if (percent < 0) percent = 0
    else if (percent > 100) percent = 100
    writeLines(paste0(tag.percent, as.character(round(percent))))
  }
  if(!is.na(message))
  {
    if (get.geap.option('geap.text.lang', 'en') == 'en')
    {
      writeLines(paste0(tag.status, engMsg))
    } else {
      writeLines(paste0(tag.status, message))
    }
  }
  invisible(0)
}

# Emits a error
.give.error <- function(message = "Error!"){
  writeLines(paste0(tag.error, message))
  stop("ERROR!")
}

# Checks if verbose is activated in arguments
.is.verbose <- function(...)
{
  na.exclude(c(as.logical(list(...)[['verbose']]), verbose.default))[1] # verbose.default : defined in vardefs.R
}

# Checks if the current OS is Windows
is.windows <- function()
{
  return(Sys.info()[['sysname']] %in% 'Windows')
}

# Loads one or multiple packages
# Use verbose=TRUE to enable library startup messages
# [[geapexec void LoadRequiredPackages(params string[] pkgNames)]]
#' @export
loadpkgs <- function(...)
{
  #argls = list(...)
  pkgs = unlist(filter.args.bytype(types = 'character', ...))
  verbose = .is.verbose(...)
  unloaded = NULL
  loadfun = NULL
  if (verbose)
  {
    loadfun = function(pkgname)
    {
      .give.status(message = sprintf("Carregando módulos: %s", pkgname), engMsg = sprintf("Loading modules: %s", pkgname))
      library(pkgname, verbose = T, character.only = T, logical.return = T)
    }
  } else {
    loadfun = function(pkgname)
    {
      library(pkgname, verbose = F, character.only = T, logical.return = T)
    }
  }
  loadedpkgs = search()
  for (pkgname in pkgs)
  {
    if (pkgname %in% loadedpkgs || sprintf("package:%s", pkgname) %in% loadedpkgs) next
    succ = suppressPackageStartupMessages(loadfun(pkgname))
    if (!succ)
    {
      unloaded[length(unloaded)+1] = pkgname
    }
  }
  if (length(unloaded) != 0) stop(sprintf("Could not load the following packages:\n%s", paste0(unloaded, collapse = ', ') ))
  invisible(T)
}

# Gets the names of all currently installed packages
# [[geapexport string[] GetAllPackages()]]
#' @export
all.packages <- function()
{
  pkgs = rownames(installed.packages(fields="Package"))
  pkgs
}

loaded.packages <- function()
{
  return(c(sub("^package:", "", search())))
}

# Checks if all the following packages are installed
# [[geapexport bool HasPackages(params string[] pkgnames)]]
#' @export
has.packages <- function(pkgnames)
{
  return(all(pkgnames %in% c(loaded.packages(), all.packages())))
}



# [[geapexport bool RObjectExists(string objName)]]
#' @export
robj.exists <- function(objName)
{
  exists(objName, where = parent.frame())
}

# [[geapexport string GetClass(call obj)]]
#' @export
robj.class <- function(obj)
{
  return(class(obj)[1])
}

# [[geapexport bool InheritsClass(call obj, params string[] whatClass)]]
#' @export
robj.inherits <- function(obj, what)
{
  return(inherits(obj, what))
}

# [[geapexport void GarbageCollect(bool full=true)]]
#' @export
garbage.collect <- function(full=T)
{
  gc(verbose = F, full = full)
  invisible(0)
}

