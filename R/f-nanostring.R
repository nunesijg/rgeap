
##########################
# NANOSTRING METHODS
# ########################
# Nunes et al, 2020
# Last updated version: 0.3.2
# 


.initialize.nanostring <- function()
{
  loadpkgs('NanoStringNorm')
  .self.oneshot()
}

# [[geapexec assign ReadNanoStringRCC(path[] fileNames)]]
#' @export
read.nanostring.RCC <- function(rccfiles)
{
  .initialize.nanostring()
  dirs = unique(dirname(rccfiles))
  basenms = basename(rccfiles)
  rawset = NULL
  for (dirnm in dirs)
  {
    exclude = list.files(dirnm, full.names = FALSE)
    exclude = exclude[!(exclude %in% basenms)]
    rset = NanoStringNorm::read.markup.RCC(dirnm, rcc.pattern = NULL, exclude = exclude)
    cat("\n")
    if (is.null(rawset))
    {
      rawset = rset
    }
    else
    {
      newcols = setdiff(colnames(rset$x), colnames(rawset$x))
      if (length(newcols) == 0L) next
      rawset$x[, newcols] = rset$x[, newcols]
      rawset$header[, newcols] = rset$header[, newcols]
    }
  }
  rawset
}

# Treats RAW data (NanoString)
# Parameters:
# - bgCorrect: none, mean, mean.2sd, max
# - sampleContent: none, housekeeping.sum, housekeeping.geo.mean, total.sum, low.cv.geo.mean, top.mean, top.geo.mean
# - otherNorm : none, quantile, zscore
# [[geapexec assign TreatNanoString(call rawSet)]]
#' @export
treat.nanostring <- function(rawset, bgCorrect='none', sampleContent='none', otherNorm='none', is.log=FALSE, ...)
{
  .initialize.nanostring()
  normset = NanoStringNorm::NanoStringNorm(rawset, Background = bgCorrect,
                                           SampleContent = sampleContent, OtherNorm = otherNorm,
                                           is.log=is.log, verbose=FALSE, ...)
  eset = df.to.eset(normset$normalized.data)
  eset
}