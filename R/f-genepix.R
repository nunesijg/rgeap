
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
# [[geapexec assign ReadGenePixGPR(path[] filenames, string method, bool singleChannel = true) ]]
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
    waveinfo = NULL
    wavepatt = NULL
    wavepatterns = c(
      Wavelengths = "(\\d+)[^\t]*\t(\\d+)[^\t]*$",
      RatioFormulation = "W1/W2\\s*?\\((\\d+?)(?:\\s*?nm)?/(\\d+?)(?:\\s*?nm)?\\)"
    )
    wavepatterns['RatioFormulations'] = wavepatterns['RatioFormulation']
    for (winf in names(wavepatterns))
    {
      if (winf %in% names(gprhead))
      {
        waveinfo = gprhead[[winf]]
        wavepatt = wavepatterns[winf]
        break
      }
    }
    if (is.null(waveinfo)) stop("This GenePix version is too old or is not supported. Contact the developer or edit the file to indicate the Wavelengths entry in the headers according to the columns (ex.: Wavelengths=635\t532)")
    regms = .regex.matches(wavepatt, waveinfo)
    capmethod = c(median='Median', mean='Mean')[method]
    if (length(regms) == 2)
    {
      gprcols = list(R = sprintf("F%s %s", regms[1], capmethod), G = sprintf("F%s %s", regms[2], capmethod), 
                     Rb = sprintf("B%s %s", regms[1], capmethod), Gb = sprintf("B%s %s", regms[2], capmethod))
    }
    else
    {
      wcolpatterns = c(R="(\\d+)\t[a-zA-Z/]+", G="[a-zA-Z/]+\t(\\d+)")
      for (colnm in names(wcolpatterns))
      {
        colpatt = wcolpatterns[colnm]
        if (grepl(colpatt, waveinfo, perl = T))
        {
          wlen = sub(colpatt, '\\1', waveinfo, perl = T)
          gprcols = list()
          gprcols[colnm] = sprintf("F%s %s", wlen, capmethod)
          gprcols[sprintf("%sb", colnm)] = sprintf("B%s %s", wlen, capmethod)
          break
        }
      }
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
  for (lslot in c('R', 'G', 'Rb', 'Gb', 'E', 'Eb'))
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
