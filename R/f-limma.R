
#' @include base.R
#' @include f-fileio.R
#' @include f-esets.R

##########################
# LIMMA METHODS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0
# 

.initialize.limma <- function()
{
  loadpkgs('Biobase', 'limma')
  .self.oneshot()
}

# Creates a new EListRaw object
# [[geapexport assign CreateEListRaw(dots args)]]
#' @export
create.elistraw <- function(...)
{
  .initialize.limma()
  rgls = new("EListRaw", list(...))
  rgls
}

# Creates a new RGList object
# [[geapexport assign CreateRGList(dots args)]]
#' @export
create.rglist <- function(...)
{
  .initialize.limma()
  rgls = new("RGList", list(...))
  rgls
}

# Treats array data (background and normalization), returns an ExpressionSet
# [[geapexec assign TreatLimma(call ergObjName, string bgMethod = "auto", string bgNormExpMethod = "saddle", double bgOffset=0, string normMethod = "quantile", string normCyclicMethod = "fast" )]]
#' @export
treat.limma <- function(ergls, bg.method = 'auto', bg.normexp.method = 'saddle', bg.offset=0, norm.method = 'quantile', norm.cyclic.method = 'fast')
{
  .initialize.limma()
  wasmat = F
  if (inherits(ergls, 'ExpressionSet'))
  {
    ergls = eset.exprs(ergls)
    wasmat = T
  }
  ergclass = class(ergls)[1]
  if (ergclass == 'matrix')
  {
    ergls = new("EListRaw", list(E = ergls))
    ergclass = 'EListRaw'
    wasmat = T
  }
  bg.method = match.arg(bg.method[1], c('auto', 'none', 'subtract', 'half', 'minimum', 'movingmin', 'edwards', 'normexp'))
  bg.normexp.method = match.arg(bg.normexp.method[1], c('saddle', 'mle', 'rma', 'rma75'))
  normoptions = c('none', 'scale', 'quantile', 'cyclicloess', 'vsn')
  if (ergclass == 'RGList') normoptions = c(normoptions, 'Aquantile', 'Gquantile', 'Rquantile', 'Tquantile') # TODO: Tquantile requer definição dos 'targets'
  norm.method = match.arg(norm.method[1], normoptions)
  norm.cyclic.method = match.arg(norm.cyclic.method[1], c('fast', 'affy', 'pairs'))
  if (bg.method != 'none')
  {
    .give.status(message = sprintf("Aplicando correção de fundo no objeto %s...", ergclass), engMsg = sprintf("Applying background correction on %s values...", ergclass))
    lsres = limma::backgroundCorrect(ergls, method = bg.method, normexp.method = bg.normexp.method, offset = bg.offset, verbose = F)
  }
  else
  {
    minval = min(ergls$E) # Workaround to avoid negative values
    if (minval <= 0)
    {
      warning(sprintf("Negative values found.\nAll values were offset by %.2f", (-minval + 1e-10)))
      ergls$E = ergls$E + (-minval + 1e-10)
    }
    lsres = ergls
  }
  .give.status(message = sprintf("Normalizando valores do objeto %s...", ergclass), engMsg = sprintf("Normalizing %s values...", ergclass))
  lsres = switch(norm.method,
                 vsn = limma::normalizeVSN(lsres),
                 limma::normalizeBetweenArrays(lsres, method=norm.method, cyclic.method = norm.cyclic.method)
                 )
  if (!wasmat)
  {
    lsres = limma::avereps(lsres, ID=lsres$genes[,1])
  }
  if (inherits(lsres, 'EList'))
  {
    validSel = !apply(lsres$E, 1, anyNA)
    eset = df.to.eset(lsres$E[validSel,,drop=F])
    if (!wasmat) eset = eset.insert.platdata(eset, lsres$genes[validSel,,drop=F])
    return(eset)
  }
  if (inherits(lsres, 'MAList'))
  {
    eset = malist.to.eset(lsres)
    return(eset)
  }
  lsres
}

# Given a MAList, returns an ExpressionSet with separated (R and G) columns
malist.to.eset <- function(mals)
{
  .initialize.limma()
  if (!inherits(mals, 'MAList')) stop("'mals' must be a MAList")
  mat = limma::exprs.MA(mals)
  dimnames(mat) = list(rownames(mals$M),
                       paste0(colnames(mals$M)[as.integer(sapply(1:ncol(mals$M), rep, 2))], c("_G", "_R")) )
  eset = df.to.eset(mat)
  if (!is.null(mals$genes) && ncol(mals$genes) != 0)
  {
    eset = eset.insert.platdata(eset, mals$genes)
  }
  eset
}
