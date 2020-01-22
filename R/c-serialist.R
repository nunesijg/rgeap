
#' @include base.R
#' @include c-multifactor.R

##########################
# SERIALIST CLASS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0
# 

# Forward definition of serializable_t
setClassUnion('serializable_atomic', c('logical', 'integer', 'numeric', 'character', 'raw'))

# Subclass for element serializing
setClass('serialelem',
         slots = c(type = 'character',
                   count = 'integer'
                   )
         )

setClass('serialelem.atomic',
         slots = c(value = 'serializable_atomic',
                   binsize = 'numeric'
                   ),
         contains = 'serialelem'
         )
setClass('serialelem.list',
         slots = c(value = 'list',
                   attrs = 'character'
                   ),
         contains = 'serialelem'
         )

# Forward definition of serializable_t
setClassUnion('serializable_list', c('matrix', 'data.frame', 'factor', 'list', 'serialelem'))
setClassUnion('serializable_t', c('serializable_atomic', 'serializable_list'))

setMethod('initialize', 'serialelem.list',
          function(.Object, value)
          {
            if (!is(value, 'serializable_list')) stop("'value' must be a serializable_list object")
            rettype = NA_character_
            attrs = character(0)
            if (is.matrix(value))
            {
              tmpls = list()
              attrs['matrixtype'] = class(value[1,1])
              attrs['ncol'] = ncol(value)
              attrs['nrow'] = nrow(value)
              if (!is.null(rownames(value)))
              {
                tmpls$row.names = new('serialelem.atomic', rownames(value))
              }
              if (!is.null(colnames(value)))
              {
                tmpls$col.names = new('serialelem.atomic', colnames(value))
              }
              tmpls$values = new('serialelem.atomic', as.vector(value))
              rettype = class(value)
              value = tmpls
              
            } else if (is.data.frame(value))
            {
              tmpls = list()
              attrs['nrow'] = as.character(nrow(value))
              if (!is.null(rownames(value)))
              {
                tmpls$row.names = new('serialelem.atomic', rownames(value))
              }
              cols = if(is.null(colnames(value))) (1L:ncol(value)) else colnames(value)
              for (nm in cols)
              {
                valcol = value[, nm]
                if (is.list(valcol) && multifactor.cancoerce(valcol, F))
                {
                  valcol = multifactor.build(valcol)
                }
                if (.isSerializableDfColumn(valcol))
                {
                  if (is.multifactor(valcol))
                  {
                    tmpls[[as.character(nm)]] = new('serialelem.list', multifactor.tolist(valcol))
                  } else if (is.factor(valcol))
                  {
                    tmpls[[as.character(nm)]] = new('serialelem.list', valcol)
                  } else if (is(valcol, 'serializable_atomic' ))
                  {
                    tmpls[[as.character(nm)]] = new('serialelem.atomic', valcol)
                  }
                }
              }
              rettype = class(value)
              value = tmpls
            } else if (is.factor(value))
            {
              tmpls = list()
              tmpls$levels = new('serialelem.atomic', levels(value))
              tmpls$indexes = new('serialelem.atomic', as.integer(value))
              value = tmpls
              rettype = 'factor'
              
            } else if (is.list(value))
            {
              tmpls = list()
              rettype = 'list'
              if (length(class(value)) > 1) rettype = class(value)[2]
              if (length(value) != 0)
              {
                for (i in 1:length(value))
                {
                  if (inherits(value[[i]], 'serialelem'))
                  {
                    tmpls[[i]] = value[[i]]
                  } else if (is(value[[i]], 'serializable_list'))
                  {
                    tmpls[[i]] = new('serialelem.list', value[[i]])
                  } else if (is(value[[i]], 'serializable_atomic'))
                  {
                    tmpls[[i]] = new('serialelem.atomic', value[[i]])
                  }
                }
                lsnames = if (is.null(names(value))) rep('', length(value)) else names(value)
                names(tmpls) = lsnames
                if (('values' %in% lsnames) && !is.null(names(attributes(value))) &&
                    !anyNA(match(c('matrixtype', 'ncol', 'nrow'), names(attributes(value)) ) ))
                {# this is a matrix
                  rettype = 'matrix'
                  attrs['matrixtype'] = attr(value, 'matrixtype')
                  attrs['ncol'] = attr(value, 'ncol')
                  attrs['nrow'] = attr(value, 'nrow')
                } else if (length(value) > 1 && 'row.names' %in% lsnames && all(sapply(value, .isSerializableDfColumn)) && length(unique(sapply(value, length))) == 1L)
                {# this is a data.frame
                  rettype = 'data.frame'
                  attrs['nrow'] = length(value[['row.names']])
                }
              }
              value = tmpls
            } else stop("Unknown error during serialization??")
            .Object@value = value
            .Object@count = length(value)
            .Object@type = rettype
            .Object@attrs = attrs
            .Object
          })

setMethod('initialize', 'serialelem.atomic',
          function(.Object, value)
          {
            if (!is(value, 'serializable_atomic')) stop("'value' must be a serializable_atomic object")
            .Object@type = class(value)
            .Object@count = length(value)
            binsize = NA_real_
            if (is.character(value))
            {
              value = charVecToRaw(value)
            }
            if (is.raw(value)) {
              if (!is.null(attr(value, 'type')))
              {
                .Object@type = match.arg(attr(value, 'type'), getDirectSubclasses('serializable_atomic'))
                if (is.null(attr(value, 'count')))
                {
                  .Object@count = 1L
                } else {
                  .Object@count = as.integer(attr(value, 'count'))
                }
              }
              binsize = as.numeric(length(value))
            } else {
              veclen = length(value)
              binsize = switch (class(value),
                      integer = veclen * 4,
                      numeric = veclen * 8,
                      logical = veclen,
                      veclen)
            }
            .Object@value = value
            .Object@binsize = binsize
            validObject(.Object)
            return(.Object)
          })

setMethod('initialize', 'serialelem',
          function(.Object, value)
          {
            if (is(value, 'serializable_list'))
            {
              .Object = new('serialelem.list', value)
            } else if (is(value, 'serializable_atomic'))
            {
              .Object = new('serialelem.atomic', value)
            } else stop("'value' must be a serializable_t object")
            validObject(.Object)
            return(.Object)
          })

setMethod('show', 'serialelem',
          function(object)
          {
            serialprint(object)
            invisible(0)
          })


setMethod('as.character', 'serialelem.list',
          function(x, ...)
          {
            retval = c()
            attrtxt = ''
            if (length(x@attrs) != 0)
            {
              attrtxt = sprintf("%s ", paste(names(x@attrs), paste0("'", x@attrs , "'" ), sep = '=', collapse = ' '))
            }
            lsnm = if (hasArg('name')) list(...)[['name']] else ''
            retval[1] = sprintf("<list name='%s' type='%s' count='%d' %s>",
                             lsnm,
                             x@type,
                             x@count,
                             attrtxt)
            if (length(x@value) != 0)
            {
              lsnames = names(x@value)
              for (i in 1:length(x@value))
              {
                retval[(length(retval) + 1)] = as.character(x@value[[i]], name=lsnames[i] )
              }
            }
            retval[(length(retval) + 1)] = "</list>"
            retval = paste0(retval, collapse = '\n')
            retval
          })

setMethod('as.character', 'serialelem.atomic',
          function(x, ...)
          {
            args = list(...)
            nm = 'x'
            if ('name' %in% names(args))
            {
              nm = args[['name']]
            }
            if (nchar(nm) == 0 || any(grepl("[^\\w\\.]", nm, perl = T)))
            {
              nm = 'x'
            }
            retval = sprintf("<elem name='%s' type='%s' count='%d' ptr='%.0f' binsize='%.0f' />",
                             nm,
                             x@type,
                             x@count,
                             .get.mem.address(x@value),
                             x@binsize)
            retval
          })

setMethod('$', c('serialelem.list'), function(x, name) getElement(x@value, name))


setMethod('$<-', c(x='serialelem.list', value='serializable_t'),
          function(x, name, value) 
          {
            if (nchar(name) == 0) stop("Element name must be a non-empty character")
            x@value[[name]] = xserialize(value, F)
            x@count = length(x@value)
            x
          })

setMethod('names', c(x='serialelem.list'), function(x) names(x@value))

setMethod('[[', c(x='serialelem.list', i='character', j='missing'), function(x, i) x@value[[i]] )

setMethod('[[', c(x='serialelem.list', i='integer', j='missing'), function(x, i) x@value[[i]] )

setMethod('[[', c(x='serialelem.list', i='numeric', j='missing'), function(x, i) x@value[[as.integer(i)]] )

setMethod('[[<-', c(x='serialelem.list', i='integer', j='missing', value='serializable_t'),
          function(x, i, value)
          {
            x@value[[i]] = xserialize(value, F)
            x@count = length(x@value)
            if (is.null(names(x@value))) names(x@value) = rep('', length(x@value))
            x
          } )
setMethod('[[<-', c(x='serialelem.list', i='numeric', j='missing', value='serializable_t'),
          function(x, i, value)
          {
            x@value[[as.integer(i)]] = xserialize(value, F)
            x@count = length(x@value)
            if (is.null(names(x@value))) names(x@value) = rep('', length(x@value))
            x
          })

setMethod('[[<-', c(x='serialelem.list', i='character', j='missing', value='serializable_t'),
          function(x, i, value)
          {
            x@value[[i]] = xserialize(value, F)
            x@count = length(x@value)
            x
          } )


setMethod('length', 'serialelem', function(x) x@count)

getDirectSubclasses <- function(classname)
{
  cdef = getClassDef(classname)
  subcls = names(cdef@subclasses)[sapply(cdef@subclasses, function(x) x@distance == 1 )]
  subcls
}

.isSerializableDfColumn <- function(column) any(class(column) %in% c('numeric', 'integer', 'logical', 'character', 'factor', 'multifactor'))
