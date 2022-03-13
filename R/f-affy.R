
#' @include base.R

##########################
# AFFY
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

.initialize.affy <- function()
{
  loadpkgs('affy', 'affyio', 'affxparser')
  affyenv = environment(expresso)
  unlockBinding("updateMe", affyenv)
  affyenv$updateMe <- function (object, ...)
  {
    .local <- function (object) 
    {
      increment <- get("increment", object@internals)
      i <- get("i", object@internals) + increment
      milestones.i <- get("milestones.i", object@internals)
      milestones <- get("milestones", object@internals)
      touched <- FALSE
      if (i == 1) cat('\n')
      while (milestones.i <= length(milestones) && i >= milestones[milestones.i])
      {
        #cat("#")
        perc <- round(100 * i / tail(milestones, n=1))
        .give.status(percent = perc,
                   message = paste0("Aplicando tratamento estatístico (Método Expresso): ", perc, '%'),
                   engMsg = paste0("Applying statistical treatment (Expresso Method): ", perc, '%'))
        milestones.i <- milestones.i + increment
        touched <- TRUE
      }
      if (touched) {
        assign("milestones.i", milestones.i, envir = object@internals)
        if (.Platform$OS.type == "windows") 
          flush.console()
      }
      assign(x="i",value = i, envir = object@internals)
    }
    .local(object, ...)
  }
  .self.oneshot()
}

.initialize.affy.customcdf <- function()
{
  loadpkgs('makecdfenv')
  .self.oneshot()
}

.initialize.affy.treatgcrma <- function()
{
  loadpkgs('gcrma')
  .self.oneshot()
}

.initialize.affy.treatplier <- function()
{
  loadpkgs('plier')
  .self.oneshot()
}


.initialize.oligo <- function()
{
  loadpkgs('oligoClasses', 'oligo')
  .self.oneshot()
}

# Checks if a CEL file is only supported by oligo
is.oligo.only.celfile <- function(celFileName)
{
  cdfnm = getCdfPkgName(celFileName)
  return(length(grep("gene2[01]st|ex[1-2][0-1]st|hta20|mta10", cdfnm)) == 1)
}

# Reads CEL files using affxparser
read.affxparser.cel.files <- function(samplefiles)
{
  .initialize.affy()
  mexprs = NULL
  cdfnm = NULL
  colnms = sub(pattern = ".CEL$", replacement = '', x = basename(samplefiles))
  nr = 0L
  nc = 0L
  i = 0L
  for (smpfile in samplefiles)
  {
    i = i + 1L
    isinit = is.null(mexprs)
    smpdata = affxparser::readCel(smpfile, readHeader = isinit, verbose = FALSE)
    if (isinit)
    {
      mexprs = matrix(NA_real_, nrow=length(smpdata$intensities), ncol=length(samplefiles),
                      dimnames=list(1L:length(smpdata$intensities), colnms))
      cdfnm = smpdata$header$chiptype
      nc = smpdata$header$cols
      nr = smpdata$header$rows
    }
    mexprs[,i] = smpdata$intensities
  }
  assayenv = new.env()
  assayenv$exprs = mexprs
  new('AffyBatch', assayData=assayenv,
      cdfName=cdfnm,
      nrow=nr,
      ncol=nc,
      annotation=affy::cleancdfname(cdfnm))
}

# [[geapexec assign AffyReadCelFiles(path[] fileNames)]]
#' @export
read.affy.cel.files <- function(samplefiles) 
{
  .initialize.affy()
  if (is.oligo.only.celfile(samplefiles[1]))
  {
    affyRaw = read.affxparser.cel.files(samplefiles)
  }
  else
  {
    affyRaw = affy::ReadAffy(filenames = samplefiles, widget = F, verbose = F)
    colnames(exprs(affyRaw)) = sub(pattern = ".CEL$", replacement = '', x = colnames(exprs(affyRaw)), ignore.case = T)
    rownames(affyRaw@phenoData@data) = sub(pattern = ".CEL$", replacement = '', x = rownames(affyRaw@phenoData@data), ignore.case = T)
    rownames(affyRaw@protocolData@data) = sub(pattern = ".CEL$", replacement = '', x = rownames(affyRaw@protocolData@data), ignore.case = T)
  }
  return(affyRaw)
}

# [[geapexec assign OligoReadCelFiles(path[] fileNames, dots optArgs)]]
#' @export
read.oligo.cel.files <- function(samplefiles, ...) 
{
  .initialize.oligo()
  smpnms = sub(pattern = ".CEL$", replacement = '', x = basename(samplefiles), ignore.case = T)
  affyRaw = oligo::read.celfiles(samplefiles, sampleNames = smpnms, verbose = FALSE, ...)
  return(affyRaw)
}


# [[geapexec assign AffyTreatAffyBatch(call affyBatch, string treatOption, string bgCorrect, string normMethod, string pmCorrect, string summaryMethod)]]
#' @export
treat.affy <- function(affyRaw, treatmentPackage, bgCorrect, normMethod, pmcorrect, summary)
{
  .initialize.affy()
  doNormalize = normMethod != "none"
  doBgcorrect = bgCorrect != "none"
  if (treatmentPackage == "expresso")
  {
    .give.status(message = "Aplicando tratamento estatístico (Expresso)", engMsg = "Applying statistical treatment (Expresso)")
    eset = affy::expresso(affyRaw,
                          bg.correct = doBgcorrect,
                          bgcorrect.method = bgCorrect,
                          normalize = doNormalize,
                          normalize.method = normMethod,
                          pmcorrect.method = pmcorrect,
                          summary.method = summary)
    
  } else if(treatmentPackage == "gcrma")
  {
    .initialize.affy.treatgcrma()
    .give.status(message = "Aplicando tratamento estatístico (GC-RMA)", engMsg = "Applying statistical treatment (GC-RMA)")
    eset = gcrma::gcrma(affyRaw,
                        normalize = doNormalize,
                        optical.correct = doBgcorrect,
                        GSB.adjust = FALSE)
    
  } else if(treatmentPackage == "plier")
  {
    .initialize.affy.treatplier()
    .give.status(message = "Aplicando tratamento estatístico (Plier)", engMsg = "Applying statistical treatment (Plier)")
    eset = plier::justPlier(affyRaw,
                            normalize = doNormalize,
                            norm.type = pmcorrect)
  }
  return(eset)
}

# [[geapexec assign OligoTreatFSet(call affyBatch, string bgCorrect, string normMethod, string summaryMethod)]]
#' @export
treat.oligo <- function(affyRaw, bgCorrect, normMethod, summaryMethod)
{
  .initialize.oligo()
  normMethod = sub("^quantiles", "quantile", normMethod)
  bgCorrect = match.arg(bgCorrect, c(oligo::backgroundCorrectionMethods(), 'none'))
  normMethod = match.arg(normMethod, c(oligo::normalizationMethods(), 'none'))
  summaryMethod = match.arg(summaryMethod, oligo::summarizationMethods())
  .give.status(message = "Aplicando tratamento estatístico (Oligo)", engMsg = "Applying statistical treatment (Oligo)")
  pmi = oligo::pmindex(affyRaw)
  pnVec = oligo::probeNames(affyRaw)
  tbls = DBI::dbListTables(db(affyRaw))
  if (oligo::manufacturer(affyRaw) == "Affymetrix" && "bgfeature" %in% tbls)
  {
    sql = paste("SELECT man_fsetid, fid", "FROM bgfeature", 
                 "INNER JOIN featureSet", "USING(fsetid)")
    tmpQcPm = DBI::dbGetQuery(db(affyRaw), sql)
    pmi = c(pmi, tmpQcPm[["fid"]])
    pnVec = c(pnVec, tmpQcPm[["man_fsetid"]])
  }
  idx = order(pnVec)
  pnVec = pnVec[idx]
  pmi = pmi[idx]
  theClass = class(exprs(affyRaw))
  if ("matrix" %in% theClass)
  {
    afb = affyRaw
    if (bgCorrect != 'none')
    {
      .give.status(message='Aplicando correção de fundo...', engMsg = 'Applying background correction...')
      if (bgCorrect == 'LESN') # Correcting a bug where R_stretch_down was not found in oligo package. Used affyPLM instead
      {
        library('affyPLM')
        pmMat = oligo::pm(afb)
        oligo::pm(afb) = matrix(.C("R_stretch_down", as.double(as.vector(pmMat)),
                                   0.25, dim(pmMat)[1], dim(pmMat)[2], 
                                   5L, 32.0, PACKAGE = "affyPLM")[[1]], 
                                dim(pmMat)[1], dim(pmMat)[2])
      }
      else
      {
        if (bgCorrect == 'mas' && is.null(mmSequence(afb)))
          stop("Cannot use the MAS background correction method due to missing mismatch (MM) table. Use RMA or LESN instead")
        afb = oligo::backgroundCorrect(afb, method = bgCorrect, verbose = FALSE)
      }
    }
    if (normMethod != 'none')
    {
      .give.status(message='Normalizando...', engMsg = 'Normalizing...')
      afb = oligo::normalize(afb, method = normMethod, verbose = FALSE)
    }
    pms = exprs(afb)[pmi, ]
    colnames(pms) =sampleNames(affyRaw)
    .give.status(message='Aplicando sumarização...', engMsg = 'Applying summarization...')
    exprs = oligo::summarize(pms, probes = pnVec, method = summaryMethod, verbose = FALSE)
    if (summaryMethod == 'plm') exprs = exprs$Estimates
    colnames(exprs) = sampleNames(affyRaw)
  }
  else
  {
    stop("basicRMA not implemented for '", theClass, 
         "' objects.")
  }
  out = new("ExpressionSet")
  slot(out, "assayData") = Biobase::assayDataNew(exprs = exprs)
  slot(out, "phenoData") = phenoData(affyRaw)
  slot(out, "featureData") = Biobase:::annotatedDataFrameFromMatrix(exprs, 
                                                                    byrow = TRUE)
  slot(out, "protocolData") = protocolData(affyRaw)
  slot(out, "annotation") = slot(affyRaw, "annotation")
  if (validObject(out))
  {
    return(out)
  }
  else
  {
    stop("Resulting object is invalid.")
  }
}


# Reads a CEL file using the affxparser header checking method
read.cel.header.core <- function(celFile)
{
  if (!file.exists(celFile)) stop("File not found: ", celFile)
  .initialize.affy()
  iscelfn = affxparser::isCelFile
  affxenv = environment(iscelfn)
  locenv = new.env(parent=affxenv)
  locenv$filename = celFile
  on.exit({if (inherits(locenv$con, 'connection') && isOpen(locenv$con)) close(locenv$con) })
  suppressWarnings(eval(functionBody(iscelfn), envir=locenv))
  if (!locenv$isCelFile) stop("Not a valid CEL file: ", celFile)
  header = locenv$header
  if (identical(header$version, 3L))
  {
    header = affyio::read.celfile.header(celFile, info = 'full', verbose = FALSE)
    header$version = 3L
  }
  else
  {
    header = affxparser::readCelHeader(celFile)
    names(header)[which(grepl('^parameters$', names(header)))] = 'algorithmparameters'
    if (is.character(header$header) && nchar(header$header) != 0L)
    {
      headsplit = grep('=', strsplit(header$header, split='\n')[[1]], value = TRUE)
      headfields = sub('=.*', '', headsplit)
      headvals = sub('.+?=(.*)', '\\1', headsplit)
      indsms = match(tolower(headfields), names(header))
      names(header)[na.omit(indsms)] = headfields[!is.na(indsms)]
      header[headfields[is.na(indsms)]] = headvals[is.na(indsms)]
    }
    defnms = c(cols='Cols', rows='Rows', total='Total',
               algorithm='Algorithm', algorithmparameters='AlgorithmParameters',
               chiptype='cdfName', datheader='DatHeader', cellmargin='CellMargin',
               noutliers='nOutliers', nmasked='nMasked')
    inds.repnms = match(names(defnms), names(header))
    sel.repnms = !is.na(inds.repnms)
    names(header)[na.omit(inds.repnms)] = defnms[sel.repnms]
    header$filename = NULL
    header$header = NULL
  }
  titlepatt = "^\\[\\d+?\\.\\.\\d+?\\]\\s+?([\\w\\-\\(\\)]+?)\\:CLS=.+"
  if (!is.null(header$DatHeader) && grepl(titlepatt, header$DatHeader, perl = TRUE))
    header$Title = sub(titlepatt, "\\1", header$DatHeader, perl=TRUE)
  
  header
}


# Reads the header information of a CEL file and returns a string array with the information in the format 'key=value'
# [[geapexport string[] ReadCELHeaderInfo(string celFilePath)]]
#' @export
read.cel.headerinfo <- function(celFile)
{
  .initialize.affy()
  celheader = read.cel.header.core(celFile)
  if (length(celheader) == 0) stop(sprintf("This CEL file (%s) has no headers", basename(celFile)))
  celheader$DatHeader = NULL
  vdims = sapply(celheader[sapply(celheader, is.numeric)], paste0, collapse=', ')
  vinfo = c(Platform=celheader$cdfName, ScanDate=celheader$ScanDate, vdims, Algorithm=celheader$Algorithm)
  headermetals = celheader[sapply(celheader, is.character)]
  headermetals$cdfName = NULL
  headermetals$AlgorithmParameters = NULL
  vinfo = c(vinfo, unlist(headermetals))
  if (is.character(celheader$AlgorithmParameters))
  {
    parls = strsplit(strsplit(celheader$AlgorithmParameters, ';', fixed=TRUE)[[1]], ':', fixed=TRUE)
    vinfo = c(vinfo, setNames(sapply(parls, function(p) c(p, '')[2L]), sapply(parls, `[`, i=1L)) )
  }
  names(vinfo) = sub(' ', '_', names(vinfo))
  vinfo = vinfo[unique(names(vinfo))]
  vinfo = paste(names(vinfo), vinfo, sep='=')
  vinfo
}

# Creates a probe-level (CDF) package from a CDF File. Returns the output directory of the created package
# [[geapexec string CDF2Package(string cdfFilePath, string pkgName, string orgName)]]
#' @export
cdf.makepkg <- function(cdfFile, pkgName, orgName)
{
  .initialize.affy.customcdf()
  dirnm = dirname(cdfFile)
  on.exit({clean.opened.streams()})
  .give.status(message = sprintf("Compilando pacote CDF (%s)", pkgName),
               engMsg = sprintf("Compiling CDF package (%s)", pkgName))
  makecdfenv::make.cdf.package(basename(cdfFile),
                               packagename = pkgName,
                               cdf.path = dirnm,
                               package.path = dirnm,
                               species = sub(' ', '_', orgName),
                               unlink = TRUE,
                               verbose = FALSE)
  outpath = file.path(dirnm, pkgName)
  outpath
}
