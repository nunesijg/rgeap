
#' @include base.R
#' @include f-fileio.R
#' @include c-microarrayinfo.R

##########################
# AGILENT METHODS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0
# 

.initialize.agilent <- function()
{
  loadpkgs('Biobase', 'limma')
  .self.oneshot()
}

# Reads an Agilent TXT file, returns a MicroarrayInfo with information regarding the sample
read.agilent.header <- function(filename)
{
  .initialize.agilent()
  headindexes = track.matching.lines(filename, c('FEPARAMS', 'STATS', 'TYPE', 'FEATURES'), prefixOnly = T)
  skip = headindexes[['FEATURES']]
  if (is.na(skip)) stop("This is not a valid Agilent sample file")
  vheaders = setNames(.read.lines.byindex(filename, indexes=headindexes, one.based=F), names(headindexes))
  headls = strsplit(vheaders, '\t', fixed=T)
  headls = lapply(headls, function(e) e[-1] )
  coltypes = setNames(c('integer', 'character', 'numeric', 'logical'), c('integer', 'text', 'float', 'boolean'))[headls$TYPE]
  colnms = headls$FEATURES
  possiblecols = c("ProbeName", "GeneName", "SystematicName",
                   "Row", "Col", "Start", "Sequence",
                   "SwissProt", "GenBank", "Primate",
                   "GenPept", "ProbeUID", "ControlType")
  #colids = which(setNames(colnms, colnms) %in% possiblecols)
  colids = .order.byref(colnms, possiblecols, keep.names=T)
  colids = colids[names(colids) %in% possiblecols]
  channels = c()
  chnprefs = c('g', 'r')
  for (pref in chnprefs)
  {
    vnm = sprintf("%sMedianSignal", pref)
    if (vnm %in% colnms)
    {
      channels[pref] = as.integer(length(channels) + 1L)
    }
  }
  headls$TYPE = NULL
  minfo = new('MicroarrayInfo',
              filenames=filename,
              headers=headls,
              idcolumns=colids,
              colnames=colnms,
              coltypes=coltypes,
              hasFullData=T,
              dataRange=c(begin=as.integer(headindexes[['TYPE']]), end=NA_integer_),
              manufacturer='Agilent',
              channels=channels
              )
  return(minfo)
  #vheaders = 
}

# Reads multiple files, returns a EListRaw (single channel) or RGList (dual channel)
# [[geapexec assign ReadAgilentTXT(path[] filenames, string method, bool singleChannel = true) ]]
#' @export
read.agilent <- function(filenames, method='median', green.only=F)
{
  .initialize.agilent()
  method = match.arg(method[1], choices = c('median', 'mean'))
  .assert$files.exist(filenames)
  hagil = read.agilent.header(filenames[1])
  chnms = sort(channelNames(hagil))
  dualch = !green.only && max(channels(hagil)) == 2
  colsvals = NA_character_
  chsuffix = c(median='MedianSignal', mean='MeanSignal')[method]
  retclass = 'EListRaw'
  if (dualch)
  {
    retclass = 'RGList'
    colsvals = c(G = sprintf("%s%s", chnms[1], chsuffix),
                  Gb = sprintf("%sBG%s", chnms[1], chsuffix),
                  R = sprintf("%s%s", chnms[2], chsuffix),
                  Rb = sprintf("%sBG%s", chnms[2], chsuffix))
  } else {
    colsvals = c(E = sprintf("%s%s", chnms[1], chsuffix),
                  Eb = sprintf("%sBG%s", chnms[1], chsuffix))
  }
  skip = range(hagil)[['begin']] + 1L
  colidnms = names(idcolumns(hagil))
  genedt = read.columns(filenames[1], colidnms, skip = skip, quote = '')[,colidnms]
  basenms = make.names(sub('\\.[^\\.]*?$', '', basename(filenames)), T)
  datals = lapply(colsvals, function(i) matrix(NA_real_, nrow = nrow(genedt), ncol = length(basenms),
                                               dimnames = list(1:nrow(genedt), basenms ) ))
  filecount = length(filenames)
  for (i in 1:filecount)
  {
    fagil = filenames[i]
    .give.status(percent= as.integer(100 * (i - 1) / filecount),
                 message=sprintf("Lendo %s (%d de %d)", basename(fagil), i, filecount),
                 engMsg=sprintf("Reading %s (%d of %d)", basename(fagil), i, filecount))
    dtvals = read.columns(fagil, colsvals, skip=skip, quote = '')
    for (nmslot in names(colsvals))
    {
      nmnew = basenms[i]
      nmcol = colsvals[nmslot]
      datals[[nmslot]][,nmnew] = dtvals[,nmcol]
    }
  }
  datals$targets = filenames
  datals$genes = genedt
  datals$source = sprintf("agilent.%s", method)
  #datals$Class = retclass
  retls = new(retclass, datals) #do.call('new', datals)
  retls
}
