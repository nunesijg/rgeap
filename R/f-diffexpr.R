
#' @include base.R
#' @include c-diffexpresults.R
#' @include f-esets.R

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

# [[geapexport assign DiffExprCompare(call seriesData, string[] expGroup, string[] ctrlGroup, string ebayesMethod, dots optargs ) ]]
#' @rdname diffexpr
#' 
#' @param seriesData a \code{matrix} or \code{data.frame} or \code{ExpressionSet}
#' @param expGroup a character vector for experimental group. Elements must be contained in column names from \code{seriesData}
#' @param ctrlGroup a character vector for control group. Elements must be contained in column names from \code{seriesData}
#' @param ebayesMethod a character parameter for Empirical Bayes method. Valid options are 'ebayes' and 'treat'
#' @param do.log2 a logical whether to \eqn{log_2}{log2} the values
#' @param adjust.method a character parameter to adjust p-values. See \code{\link[stats]{p.adjust}.methods}
#' 
#' @return a \code{\linkS4class{DifExpResults}} object is returned, containing data to be processed by GEAP
#' @export
diffexpr.compare <- function(seriesData, expGroup, ctrlGroup, ebayesMethod = 'ebayes', ...)
{
  .initialize.diffexpr()
  ebayesMethod = match.arg(tolower(ebayesMethod), c('ebayes', 'treat'))
  compMatrix = NULL
  if (is.data.frame(seriesData))
  {
    compMatrix = as.matrix(seriesData[,column.classes(seriesData) == 'numeric'])
  } else if (is.matrix(seriesData))
  {
    compMatrix = seriesData
  } else if (inherits(seriesData, 'ExpressionSet'))
  {
    compMatrix = eset.exprs(seriesData)
  } else {
    compMatrix = tryCatch(eset.exprs(seriesData), error=function(e) NULL)
    if (is.null(compMatrix))
    stop("Data must be a matrix, data.frame or ExpressionSet")
  }
  if (is.null(compMatrix) || !is.numeric(compMatrix))
  {
    stop("Data must consist of numeric values")
  }
  if (is.factor(expGroup)) expGroup = as.character(expGroup)
  if (is.numeric(expGroup)) expGroup = colnames(compMatrix)[expGroup]
  if (is.factor(ctrlGroup)) ctrlGroup = as.character(ctrlGroup)
  if (is.numeric(ctrlGroup)) ctrlGroup = colnames(compMatrix)[ctrlGroup]
  
  doLog2 = ...arg(do.log2, F)
  
  if (ebayesMethod == 'ebayes')
  {
    adjMethod = ...arg(adjust.method, 'BH')
    ebayesCall = substitute(eBayes(groupContrast))
    tableCall = substitute(topTable(difExpData,coef=1, number=1000000, p.value = 1, sort.by = "none", adjust.method = adjMethod))
  }
  else
  {
    adjMethod = 'TREAT'
    lfc = ...arg(lfc, log2(1.2))
    ebayesCall = substitute(treat(groupContrast, lfc=lfc))
    tableCall = substitute(topTreat(difExpData,coef=1, number=1000000, p.value = 1, sort.by = "none"))
  }
  
  groupcols = unique(c(expGroup, ctrlGroup))
  compMatrix = compMatrix[,groupcols]
  compMatrix = eset.rm.invalid.values(compMatrix) # Removes the invalid values, since they will return errors
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
    difExpData = eval(ebayesCall)
    # Gets the data.frame of results
    results = eval(tableCall)
    results = cbind(ID = rownames(results), results)
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
  resultList = new('DifExpResults', resultsData=results, plotData=plotDisplayVals,
                   expGroup=expGroup, ctrlGroup=ctrlGroup, adjMethod=adjMethod )
  return(resultList)
}



