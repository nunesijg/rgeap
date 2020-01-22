
#' @include base.R
#' @include c-diffexpresults.R

##########################
# Differential Expression
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0
# 

.initialize.diffexpr <- function()
{
  loadpkgs('Biobase', 'limma')
  .self.oneshot()
}

#' @rdname diffexpr
#' @name diffexpr.compare
#' @title Differential Expression Analysis
#' 
#' @description Compares the differential expression between experiment and control samples.
#' 
#' @seealso \linkS4class{DifExpResults}
NULL

# [[geapexport assign DiffExprCompare(call seriesData, string[] expGroup, string[] ctrlGroup, bool doLog2=false, string adjMethod="BH" ) ]]
#' @rdname diffexpr
#' 
#' @param seriesData a \code{matrix} or \code{data.frame} or \code{ExpressionSet}
#' @param expGroup a character vector for experimental group. Elements must be contained in column names from \code{seriesData}
#' @param ctrlGroup a character vector for control group. Elements must be contained in column names from \code{seriesData}
#' @param doLog2 a logical whether to \eqn{log_2}{log2} the values
#' @param adjMethod a character parameter to adjust p.values. See \code{\link[stats]{p.adjust}.methods}
#' 
#' @return a \code{\linkS4class{DifExpResults}} object is returned, containing data to be processed by GEAP
#' @export
diffexpr.compare <- function(seriesData, expGroup, ctrlGroup, doLog2 = F, adjMethod = 'BH')
{
  .initialize.diffexpr()
  compMatrix = NULL
  if (is.data.frame(seriesData))
  {
    compMatrix = as.matrix(seriesData[,column.classes(seriesData) == 'numeric'])
  } else if (is.matrix(seriesData))
  {
    compMatrix = seriesData
  } else if (inherits(seriesData, 'ExpressionSet'))
  {
    compMatrix = exprs(seriesData)
  } else {
    stop("Data must be a matrix, data.frame or ExpressionSet")
  }
  if (is.null(compMatrix) || !is.numeric(compMatrix))
  {
    stop("Data must consist of numeric values")
  }
  groupcols = unique(c(expGroup, ctrlGroup))
  compMatrix = compMatrix[,groupcols]
  if (doLog2) compMatrix = log2(compMatrix)
  rowCount = nrow(compMatrix)
  plotDisplayVals = data.frame(ID = rownames(compMatrix), exp = NA_real_[1:rowCount], ctrl = NA_real_[1:rowCount],
                                logFC = NA_real_[1:rowCount], nlogpval = NA_real_[1:rowCount], row.names = rownames(compMatrix))
  # Comparison with 3 or more samples, producing logFC and p-value columns
  if(ncol(compMatrix) > 2)
  {
    ##### Matrix composed only by experiment and control factors
    mdesign = ifelse(groupcols %in% expGroup, 1, 2)
    # Groups' design based on matrix design
    designGroups = model.matrix(~ 0+factor(mdesign))
    # Group names (we only have two)
    colnames(designGroups) = c("g1","g2")
    # Sample names must correspond to the groups
    rownames(designGroups) = groupcols
    # Groups design is remade to allow overlaps between samples
    designGroups[,1:2] = 0
    designGroups[rownames(designGroups) %in% ctrlGroup, 1] = 1
    designGroups[rownames(designGroups) %in% expGroup, 2] = 1
    # Builds a linear model from groups design
    groupReplicates = lmFit(compMatrix, designGroups)
    # Makes a contrast matrix from the linear model
    contrastMatrix = makeContrasts(g2-g1, levels=designGroups)
    # Groups contrasts model
    groupContrast = contrasts.fit(groupReplicates, contrastMatrix)
    # Gets the differentially expressed genes
    difExpData=eBayes(groupContrast)
    # Gets the data.frame of results
    results<-topTable(difExpData,coef=1,number=1000000, p.value = 1, sort.by = "none", adjust.method = adjMethod)
    results<-cbind(ID = rownames(results), results)
    if (doLog2)
    {
      plotDisplayVals$exp = 2^(groupReplicates$coefficients[,"g2"])
      plotDisplayVals$ctrl = 2^(groupReplicates$coefficients[,"g1"])
    } else {
      plotDisplayVals$exp = groupReplicates$coefficients[,"g2"]
      plotDisplayVals$ctrl = groupReplicates$coefficients[,"g1"]
    }
    plotDisplayVals$logFC = results$logFC
    plotDisplayVals$nlogpval = -log10(results$adj.P.Val)
  } else if (ncol(compMatrix) == 2)
  { # When only two samples are compared
    results = data.frame(ID = rownames(compMatrix), logFC = rep(NA, rowCount), AveExpr = rep(NA, rowCount), row.names = rownames(compMatrix))
    plotDisplayVals$exp = compMatrix[, 1]
    plotDisplayVals$ctrl = compMatrix[, 2]
    for (i in 1:rowCount)
    {
      results[i,"logFC"] = log2(compMatrix[i, 1]) - log2(compMatrix[i, 2])
      results[i, "AveExpr"] = (compMatrix[i, 1] + compMatrix[i, 2])/2
    }
  }
  #if (adjMethod == 'none' && 'adj.P.Val' %in% names(results))
  #  results = results[, !names(results) %in% c('adj.P.Val')]
  resultList = new('DifExpResults', resultsData=results, plotData=plotDisplayVals,
                   expGroup=expGroup, ctrlGroup=ctrlGroup, adjMethod=adjMethod )
  #resultList[["displayvalues"]] = plotDisplayVals # Slots: exp, ctrl, logFC, nlogpval
  #resultList[["results"]] = results # Matriz de resultados que conhecemos
  #resultList[["difexpressas"]] <- difExpressas # Deprecated
  return(resultList)
}



