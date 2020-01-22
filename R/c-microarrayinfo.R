
#' @include base.R

##########################
# MICROARRAYINFO CLASS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0
# 


#' @name MicroarrayInfo-class
#' @rdname MicroarrayInfo
#' @title Microarray File Information
#' 
#' @description Stores metadata related to microarray files.
#' Used internally by GEAP for optimizations, and has no intended use as a standalone object.
#' 
#' @slot filenames character vector of file names
#' @slot headers list of character vectors with heading information
#' @slot rownames row names in data matrix
#' @slot colnames column names in data matrix
#' @slot coltypes column classes described for each data column
#' @slot hasFullData logical. TRUE if full microarray data table is available
#' @slot dataRange named integer vector. 'begin' is the starting index of data table values, and 'end' is the final data index
#' @slot channels named integer, where names are the channel initials (e.g. 'g' and 'r'). Values are one-based indexedes
#' @slot params named character (optional). Additional parameters, when present
#' @slot meta list of extra attributes to be considered
#' 
#' @export
setClass('MicroarrayInfo',
         slots = c(filenames = 'character',
                   headers = 'list',
                   idcolumns = 'integer', # Names, if possible
                   rownames = 'character',
                   colnames = 'character',
                   coltypes = 'character',
                   hasFullData = 'logical',
                   dataRange = 'integer', # Named, if possible
                   channels ='integer', # Named, if possible
                   manufacturer = 'character',
                   params = 'character', # Named, if possible
                   meta = 'list')
)

setMethod('initialize', 'MicroarrayInfo',
          function(.Object, 
                   filenames = NA_character_,
                   headers = list(),
                   idcolumns = integer(0),
                   rownames = character(0),
                   colnames = character(0),
                   coltypes = character(0),
                   hasFullData = FALSE,
                   dataRange = c(begin=NA_integer_, end=NA_integer_),
                   channels = setNames(1L, 'g'),
                   manufacturer = "Unknown",
                   params = character(0),
                   meta = list(),
                   ...
                   )
          {
            idcolnames = idcolumns
            if (length(idcolnames) != 0)
            {
              if (is.integer(idcolnames) || is.numeric(idcolnames))
              {
                idcolnames = colnames[as.integer(idcolnames)]
              } else {
                idcolumns = setNames((1L:length(colnames)), colnames)[idcolumns]
              }
            }
            if (hasFullData && !all(idcolnames %in% colnames)) stop("'idcolumns' elements must be part of 'colnames' if full data is available")
            if (length(params) != 0 && (is.null(names(params)) || any(duplicated(names(params))) ) ) stop("'params', if present, must be a named character with unique fields")
            if (!all(c('begin', 'end') %in% names(dataRange))) stop("'dataRange' must be a named integer vector with 'begin' and 'end' named elements")
            .Object@filenames = filenames
            .Object@headers = headers
            .Object@rownames = rownames
            .Object@colnames = colnames
            .Object@idcolumns = idcolumns
            .Object@coltypes = coltypes
            .Object@hasFullData = hasFullData
            .Object@dataRange = dataRange
            .Object@channels = channels
            .Object@manufacturer = manufacturer
            .Object@params = params
            .Object@meta = meta
            validObject(.Object)
            return(.Object)
          })

setGeneric('channels', function(minfo) standardGeneric('channels'))
setMethod('channels', 'MicroarrayInfo', function(minfo) minfo@channels )
setGeneric('meta', function(object) standardGeneric('meta'))
setMethod('meta', 'MicroarrayInfo', function(object) object@meta )

setGeneric('parameters', function(object) standardGeneric('parameters'))
setMethod('parameters', 'MicroarrayInfo', function(object) object@params )

setMethod('dim', 'MicroarrayInfo', function(x) c(length(x@rownames), length(x@colnames) ) )
setMethod('dimnames', 'MicroarrayInfo', function(x) list(x@rownames, x@colnames) )
setMethod('names', 'MicroarrayInfo', function(x) x@filenames )
setMethod('names<-', 'MicroarrayInfo', function(x, value) x@filenames <- value)
setMethod('$', 'MicroarrayInfo', function(x, name) slot(x, name))
setMethod('range', c('MicroarrayInfo', 'logical'), function(x, na.rm=F) if (na.rm) x@dataRange[!is.na(x@dataRange)] else x@dataRange )

if (!isGeneric('channelNames', where = parent.frame())) setGeneric('channelNames', function(object, ...) standardGeneric('channelNames'))
setMethod('channelNames', 'MicroarrayInfo', function(object, ...) names(object@channels) )
setGeneric('idcolumns', function(object) standardGeneric('idcolumns'))
setMethod('idcolumns', 'MicroarrayInfo', function(object) object@idcolumns )

