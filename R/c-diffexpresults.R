
#' @include base.R

##########################
# DIFEXPRESULTS CLASS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0
# 

#' @name DifExpResults-class
#' @rdname DifExpResults
#' @aliases DifExpResults
#' @title Differential Expression Comparison Results
#' 
#' @description Container for output data from differential expression analyses.
#' The contents stored in a \code{DifExpResults} represent a comparison result of type "Experimental vs Control".
#' 
#' @slot resultsData \code{data.frame} of results data, including \emph{logFC} and \emph{P.Value} columns. See details below.
#' @slot plotData \code{matrix} with data to be plotted in GEAP executable
#' @slot adjMethod character, parameter used to correct p-values.
#' @slot expGroup character vector, sample names from experimental group
#' @slot ctrlGroup character vector, sample names from control group
#' 
#' @seealso \link{diffexpr.compare}
#' 
#' @export
setClass('DifExpResults',
         slots = c(resultsData = 'data.frame',
                   plotData = 'matrix',
                   adjMethod='character',
                   expGroup = 'character',
                   ctrlGroup = 'character')
         )

setMethod('initialize', 'DifExpResults',
          function(.Object, resultsData=data.frame(), plotData=data.frame(), adjMethod='none', expGroup=NULL, ctrlGroup=NULL)
          {
            rownames(resultsData) = as.character(resultsData[,1])
            .Object@resultsData = resultsData
            if (inherits(plotData[,1], c('character', 'factor')))
            {
              rownames(plotData) = as.character(plotData[,1])
              plotData = plotData[,-1]
            }
            .Object@plotData = as.matrix(plotData)
            .Object@adjMethod = adjMethod
            .Object@expGroup = expGroup
            .Object@ctrlGroup = ctrlGroup
            return(.Object)
          })

setMethod('show', 'DifExpResults',
          function(object)
          {
            cat('Differential Expression Results'); cat('\n')
            cat(paste(sprintf('Columns: (%d)', ncol(resultsData(object)) ))); cat('\n')
            cat(paste(sprintf('Rows: (%d)', nrow(resultsData(object))) )); cat('\n')
          })

# [[geapgeneric void DifExpResultsData(call resObjName)]]
#' @export
setGeneric('resultsData', function(object) standardGeneric('resultsData'))
setMethod('resultsData', 'DifExpResults', function(object) object@resultsData )

# [[geapgeneric void DifExpPlotData(call resObjName)]]
#' @export
setGeneric('plotData', function(object) standardGeneric('plotData'))
setMethod('plotData', 'DifExpResults', function(object) object@plotData )

setMethod('dim', 'DifExpResults', function(x) dim(x@resultsData) )
setMethod('dimnames', 'DifExpResults', function(x) dimnames(x@resultsData) )
setMethod('names', 'DifExpResults', function(x) colnames(x@resultsData) )
setMethod('$', 'DifExpResults', function(x, name) getElement(resultsData(x), name))

# [[geapgeneric robj_RDataFrame DifExpTopTable(call resObjName, params named_call[] optArgs)]]
#' @export
setGeneric('top.table', function(x, ...) standardGeneric('top.table'))
setMethod('top.table', 'DifExpResults', function(x, ...)
          {
            resdt = resultsData(x)
            args = list(...)
            bycol = if (hasArg('splitcol')) args[['splitcol']] else intersect(c('Gene Symbol', 'Symbol', 'ID'), colnames(x))
            ordercol = if (hasArg('orderby')) args[['orderby']] else 'logFC'
            pvalue = if (hasArg('p.value')) args[['p.value']] else 1
            splitsep = if (hasArg('sep')) args[['sep']] else NULL
            issep = is.character(splitsep) && length(splitsep) == 1
            resdt = resdt[resdt$adj.P.Val < pvalue,,drop=F]
            if (hasArg('limit') && is.numeric(args[['limit']]))
            {
              resdt = head.tail(resdt, ordercol, args[['limit']])
            } else {
              resdt = head.tail(resdt, ordercol, nrow(resdt))
            }
            if (nrow(resdt) == 0) return(resdt)
            resdt
          })


