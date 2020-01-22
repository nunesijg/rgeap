
#' @include base.R

##########################
# Results Helpers
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

# Returns the lines of results with greater and lesser logFC, ordered by logFC value
summary.results.single <- function(resultsMatrix, platformMatrix = NULL, symbolColName = "", maxRows = 10)
{
  sumMatrix = head.tail(resultsMatrix[,c("ID", "logFC")], "logFC", decOrder = T, maxRows = maxRows)
  if (is.null(platformMatrix) || nchar(symbolColName) == 0)
  {
    return(sumMatrix)
  }
  sumMatrix = cbind(platformMatrix[rownames(sumMatrix), symbolColName], sumMatrix)
  colnames(sumMatrix)[1] = paste0(symbolColName, "/ID")
  sumMatrix = as.matrix(sumMatrix)
  for(i in 1:nrow(sumMatrix))
  {
    if (nchar(as.character(sumMatrix[i,1])) == 0)
      sumMatrix[i,1] = sumMatrix[i,2]
  }
  return(sumMatrix[,-2])
}

.split.multigenes <- function(genes, split=' /// ')
{
  genes = as.character(genes)
  genes = unique(unlist( strsplit(genes, split, fixed=T) ) )
  genes
}
