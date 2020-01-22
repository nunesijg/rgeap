#' @include base.R

##########################
# MULTIFACTOR S3 CLASS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0
# 

# Rebuilds a multifactor type
# [[geapexport assign MultiFactorRebuild(call indexes, call sizes, call levels)]]
#' @rdname serial
#' @export
multifactor.rebuild <- function(indexes, sizes, lvls)
{
  retls = list()
  class(retls) = c('list', 'multifactor')
  if (length(sizes) == 0) return(retls)
  retls = multifactor.rebuild.base(as.integer(indexes), as.integer(sizes), lvls)
  retls
}

multifactor.build <- function(strvec, sep = ' /// ')
{
  if (is.multifactor(strvec)) return(strvec)
  if (multifactor.cancoerce(strvec, T))
  {
    multf = strvec
  } else {
    if (is.factor(strvec)) strvec = as.character(strvec)
    multf = strsplit(strvec, split = sep, fixed = T, perl = F)
  }
  class(multf) = c('list', 'multifactor')
  multf
}

multifactor.build.ifneeded <- function(vec, sep)
{
  strvec = vec
  if (is.factor(strvec)) strvec = as.character(vec)
  if (!is.character(strvec)) return(vec)
  sepinds = grepl(sep, vec, perl = F, fixed = T)
  if (!any(sepinds)) return(vec)
  multf = multifactor.build(strvec, sep)
  multf
}

multifactor.tolist <- function(multf)
{
  lvls = unique(unlist(multf))
  indxls = sapply(multf, match, lvls)
  lens = sapply(indxls, length)
  retls = list(levels=lvls, indexes=unlist(indxls), sizes=lens)
  class(retls) = c("list", "multifactor", "multifactor.serialized")
  retls
}

is.multifactor <- function(multf) is(multf, "multifactor")

multifactor.cancoerce <- function(mfls, fail=F)
{
  if (is.multifactor(mfls)) return(T)
  if (!is.list(mfls)) return(F)
  isallchars = !(all(sapply(mfls, is.character)) %in% c(F, NA))
  if (fail && !isallchars) stop("multifactor must be built from a list of characters")
  isallchars
}


# #' @export
# print.multifactor.indexes <- function(x, ...) print(paste0('[', paste0(x, collapse = '; '), ']'))


# Class for indexing multiple characters in each data.frame column
#setClass('indexlist', slots=c(values = 'integer'))
#setMethod('initialize', 'indexlist',
          #function(.Object, values)
          #{
            #if (!is.integer(values) || anyNA(values)) stop("values must be valid integer values")
            #.Object@values = values
            #validObject(.Object)
            #.Object
          #})

#
# Class for indexing multiple characters in each data.frame column
# setClass('multifactor',
#          slots = c(
#            indexlist = 'list',
#            levels = 'character'
#            )
#          )
# 
# setMethod('initialize', 'multifactor',
#           function(.Object, indexlist, levels)
#           {
#             if (!is.list(indexlist)) stop("indexlist must be a list of integers")
#             if (!all(sapply(indexlist, inherits, c('integer', 'indexlist')))) stop("multifactor only accepts integers")
#             maxind = length(levels)
#             if (!all(sapply(indexlist, function(i) !any(i <= 0L | i > maxind ))) ) stop("all integers in indexlist must be inside range of levels")
#             #indexlist = sapply(indexlist, function(e) { if (!inherits(e, 'indexlist')) new('indexlist', values=e) else e } )
#             .Object@indexlist = indexlist
#             .Object@levels = levels
#             validObject(.Object)
#             .Object
#           })
# 
# levels.multifactor <- function(x) x@levels

# #' @export
# as.data.frame.multifactor <- function(x, ...)
# {
#   if (length(x@indexlist) == 0) return(data.frame())
#   nms = names(x@indexlist)
#   if (is.null(nms)) nms = as.character(1:length(x@indexlist))
#   colnm = 'x'
#   if (sys.nframe() > 1 && grepl('as.data.frame', as.character(sys.call(-1)), fixed = T))
#   {
#     penv = parent.frame()
#     colnm = penv$vnames[[penv$i]]
#   }
#   #if (ls(parent.frame()) )
#   #print(parent.frame()$vnames)
#   dt=data.frame(row.names=nms)
#   dt = do.call('$<-', list(dt, colnm, x@indexlist))
#   attr(dt, 'levels') = x@levels
#   dt
# }

# #' @export
# setMethod('[[', c(x='multifactor', i='numeric', j='missing'), function(x, i) { x@levels[x@indexlist[[i]]] }  )
# setMethod('[', c(x='multifactor', i='numeric', j='missing'), function(x, i) { x[[i]] }  )
# 
# setMethod('as.character', c(x='multifactor'), function(x, ...) sapply(x@indexlist, function(l) x@levels[l] ) )
# 
# setMethod('as.list', c(x='multifactor'), function(x, ...) list(indexlist=(x@indexlist), levels=(x@levels)) )
# 
# setMethod('length', c(x='multifactor'), function(x) length(x@indexlist) )
# 
# setMethod('$<-', c(x='data.frame', value='multifactor'), function(x, name, value)
#   {
#   x[,name] <- value@indexlist
#   x
# })
# setMethod('[<-', c(x='data.frame', i='missing', j='character', value='multifactor'), function(x, j, value)
# {
#   x[,j] <- value@indexlist
#   x
# })

#setMethod('as.data.frame', c(x='multifactor', row.names='missing', ), function(x, ...) sapply(x@indexlist, function(l) x@levels[l] )   )



# #' @export
# multifactor <- function(x=list(), levels)
# {
#   if (is.null(x))
#   {
#     x = list()
#   }
#   if (missing(levels))
#   {
#     if (all(sapply(x, inherits, 'character')))
#     {
#       levels = unique(unlist(x))
#       x = sapply(x, match, levels)
#     } else stop("x must be an entire character list when 'levels' arguments is missing")
#   }
#   multf = new('multifactor', indexlist = x, levels = levels)
# }
