
#' @include base.R
#' @include f-strformat.R

##########################
# ASSERTION METHODS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0
#

# Environment declaration
.assert <- new.env()

# Definitions
local({
  match.types = function(object, types, strict=F, okNULL=F, okNA=F) 
  {
    argnm = as.character(substitute(object))
    if (!okNULL && is.null(object)) stop(sprintf("'%s' must not be NULL", argnm))
    if (!okNA && (is.list(object) || is.vector(object)) && is.na(object)) stop(sprintf("'%s' must not be NA", argnm))
    objcl = class(object)[1]
    if (strict)
    {
      if (objcl %in% types) return(T)
    } else {
      if (inherits(object, types)) return(T)
    }
    
    if (length(types) == 1)
    {
      pronom = .pronom.an(types)
      stop(sprintf("'%s' must be %s %s", argnm, pronom, types))
    }
    stop(sprintf("'%s' must be one of the folowing types: %s\n'%s' given", argnm, paste0(types, collapse = ', '), objcl ))
    invisible(F)
  }
  files.exist = function(filenames)
  {
    argnm = as.character(substitute(filenames))
    if (is.null(filenames) || any(is.na(filenames))) stop(sprintf("'%s' cannot be NULL or NA", argnm))
    vass = file.exists(filenames)
    if (any(!vass)) stop(sprintf("Could not find the following files: %s", paste0(basename(filenames[!vass]), collapse = ", " ) ))
    invisible(T)
  }
}, envir = .assert)
