
#' @include base.R
#' @include f-fileio.R
#' @include f-limma.R

##########################
# GENEPIX METHODS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0
# 

.initialize.genepix <- function()
{
  loadpkgs('Biobase', 'limma')
  .self.oneshot()
}

# Reads multiple files, returns a EListRaw (single channel) or RGList (dual channel)
# [[geapexec assign ReadGenePixGPR(string[] filenames, string method, bool singleChannel = true) ]]
#' @export
read.genepix <- function(filenames, method='mean', green.only=F)
{
  if (length(filenames) == 0) stop("At least one file name must be provided")
  .initialize.genepix()
  method = match.arg(method[1], choices = c('median', 'mean'))
  srcmethod = sprintf("genepix.%s", method)
  gprmrg = NULL
  filecount = length(filenames)
  for (i in 1:filecount)
  {
    fnm = filenames[i]
    .give.status(percent= as.integer(100 * (i - 1) / filecount),
                 message=sprintf("Lendo %s (%d de %d)", basename(fnm), i, filecount),
                 engMsg=sprintf("Reading %s (%d of %d)", basename(fnm), i, filecount))
    gprhead = limma::readGPRHeader(fnm)
    gprcols = NULL
    regms = .regex.matches("(\\d+)[^\t]*\t(\\d+)[^\t]*$", gprhead$Wavelengths)
    if (length(regms) == 2)
    {
      capmethod = c(median='Median', mean='Mean')[method]
      gprcols = list(R = sprintf("F%s %s", regms[1], capmethod), G = sprintf("F%s %s", regms[2], capmethod), 
                     Rb = sprintf("B%s %s", regms[1], capmethod), Gb = sprintf("B%s %s", regms[2], capmethod))
    }
    gprls = limma::read.maimages(fnm, source = srcmethod, green.only = green.only, columns = gprcols, verbose = F)
    if (is.null(gprmrg))
    {
      gprmrg = gprls
    } else {
      gprmrg = cbind(gprmrg, gprls)
    }
  }
  basenms = make.names(sub('\\.[^\\.]*?$', '', basename(filenames)), T)
  gprmrg$targets = basenms
  for (lslot in c('R', 'G', 'Rb', 'Gb'))
  {
    if (lslot %in% names(gprmrg))
    {
      colnames(gprmrg[[lslot]]) = basenms
    }
  }
  if ('ID' %in% names(gprmrg$genes))
  {
    gprmrg$genes = gprmrg$genes[, c('ID', setdiff(names(gprmrg$genes), 'ID'))]
    vidcol = gprmrg$genes$ID
    naids = is.na(vidcol)
    if (any(naids))
    {
      gprmrg = gprmrg[!naids, ]
    }
  }
  gprmrg
}
