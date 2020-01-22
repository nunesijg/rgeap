
#' @include base.R
#' @include f-fileio.R

##########################
# ILLUMINA FUNTIONS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0


.initialize.illumina <- function()
{
  loadpkgs('illuminaio', 'beadarray', 'limma')
  .self.oneshot()
}

# Reads an Illumina data matrix file, returns an ExpressionSetIllumina
# [[geapexport assign ReadIlluminaMatrix(string fileName)]]
#' @export
read.illumina.matrix <- function(filename)
{#TODO: There is ...?
  .initialize.illumina()
  skip = 0L
  colshds = NA_character_
  (con = .open.filecon(filename)) %using% {
    repeat {
      cline = readLines(con, n=1L)
      if (!startsWith(cline, "#"))
      {
        colshds = strsplit(cline, '\t', fixed = T)[[1]]
        break
      }
      skip = skip + 1L
    }
  }
  colid = colshds[1]
  colsls = list()
  if( any(grepl("AVG_Signal$", colshds, perl=TRUE)) ) colsls$exprs = "AVG_Signal"
  if( any(grepl("BEAD_STDERR$", colshds, perl=TRUE)) ) colsls$se.exprs = "BEAD_STDERR"
  if( any(grepl("Avg_NBEADS$", colshds, perl=TRUE)) ) colsls$nObservations = "Avg_NBEADS"
  if( any(grepl("Detection Pval$", colshds, perl=TRUE)) ) colsls$Detection = "Detection Pval"
  dtret = readBeadSummaryData(filename, skip=skip, ProbeID = colid, columns = colsls)
  dtret
}

# Reads one or more Illumina data matrix files, returns an ExpressionSetIllumina
# [[geapexport assign ReadIlluminaMatrices(params string[] fileNames)]]
#' @export
read.illumina.matrices <- function(filenames)
{
  eset = NULL
  for (fname in filenames)
  {
    if (is.null(eset))
    {
      eset = read.illumina.matrix(fname)
    } else {
      eset = esets.merge(eset, read.illumina.matrix(fname))
    }
  }
  eset
}

# Normalizes an ExpressionSetIllumina
# [[geapexport assign TreatIlluminaMatrix(call objName, string normMethod = "quantile", string transform = "none")]]
#' @export
treat.illumina.matrix <- function(elumi, norm.method = 'quantile', transform = 'none')
{
  .initialize.illumina()
  norm.method = match.arg(norm.method[1], c('quantile', 'qspline', 'vsn', 'rankInvariant', 'median', 'none'))
  transform = match.arg(transform[1], c('none', 'log2', 'neqc', 'rsn', 'vst'))
  elumi = normaliseIllumina(elumi, norm.method, transform)
  elumi
}

# Reads a BGX file, returns a data.frame
# [[geapexport assign ReadIlluminaBGX(string fileName)]]
#' @export
read.illumina.bgx <- function(filename)
{
  .initialize.illumina()
  bgx <- readBGX(filename)
  bgx <- bgx$probes
  rownames(bgx) <- bgx$Probe_Id
  return(bgx)
}

# [[geapexport robj_RList CheckIlluminaIDAT(string fileName, int limitIds=1000)]]
#' @export
check.illumina.idat <- function(filename, limitids=1000L)
{
  .initialize.illumina()
  idatdt = readIDAT(filename)
  retls = list(barcode = idatdt$Barcode,
               arrayids = head(idatdt$Quants$IllumicodeBinData, limitids) )
  retls
}


# Exemplos: GSE80080 e GSE75550
# [[geapexec assign ReadIlluminaIDAT(string[] idatFiles, string bgxfile, string[] annotations, int tolerance=0)]]
#' @export
read.illumina.idat <- function(idatfiles, bgxfile, annotation='Symbol', tolerance = 0L)
{
  .initialize.illumina()
  .give.status(message="Lendo arquivos IDAT", engMsg="Reading IDAT files")
  retls = limma::read.idat(idatfiles, bgxfile, annotation = annotation, tolerance = tolerance, verbose = F)
  if (inherits(retls, c('EListRaw', 'EList'))) colnames(retls$E) = make.unique(basename(colnames(retls$E)), '_')
  retls
  # OK <- requireNamespace("illuminaio", quietly = TRUE)
  # if (!OK) 
  #   stop("illuminaio package required but is not installed (or can't be loaded)")
  # idatafiles <- as.character(idatfiles)
  # nsamples <- length(idatfiles)
  # elist <- new("EListRaw")
  # elist$source <- "illumina"
  # elist$targets <- data.frame(IDATfile = idatfiles, stringsAsFactors = FALSE)
  # .give.status(message="Lendo arquivo BGX", engMsg="Reading manifest file (BGX)")
  # bgx = NULL
  # if (length(bgxfile) != 0)
  # {
  #   bgx = illuminaio::readBGX(bgxfile)
  # }
  # nregprobes <- nrow(bgx$probes)
  # nctrlprobes <- nrow(bgx$control)
  # nprobes <- nregprobes + nctrlprobes
  # elist$genes <- rbind(bgx$probes[, c("Probe_Id", "Array_Address_Id")], 
  #                      bgx$controls[, c("Probe_Id", "Array_Address_Id")])
  # elist$genes$Status <- "regular"
  # elist$genes$Status[(nregprobes + 1):nprobes] <- bgx$controls[, 
  #                                                              "Reporter_Group_Name"]
  # if (!is.null(annotation) && !is.na(annotation)) {
  #   annotation <- as.character(annotation)
  #   annotation <- intersect(names(bgx$probes), annotation)
  # }
  # if (length(annotation)) {
  #   ac <- annotation %in% names(bgx$controls)
  #   for (i in 1:length(annotation)) {
  #     elist$genes[[annotation[i]]] <- NA_character_
  #     elist$genes[[annotation[i]]][1:nregprobes] <- bgx$probes[[annotation[i]]]
  #     if (ac[i]) 
  #       elist$genes[[annotation[i]]][(nregprobes + 1L):nprobes] <- bgx$controls[[annotation[i]]]
  #   }
  # }
  # elist$E <- matrix(0, nprobes, nsamples)
  # elist$E[] <- NA
  # colnames(elist$E) <- removeExt(idatfiles)
  # rownames(elist$E) <- elist$genes[, "Array_Address_Id"]
  # elist$other$STDEV <- elist$other$NumBeads <- elist$E
  # for (j in 1:nsamples) {
  #   if (verbose) 
  #     cat("\t", idatfiles[j], "... ")
  #   tmp <- illuminaio::readIDAT(idatfiles[j])
  #   if (verbose) 
  #     cat("Done\n")
  #   if ("IllumicodeBinData" %in% colnames(tmp$Quants)) {
  #     ind <- match(elist$genes$Array_Address_Id, tmp$Quants$IllumicodeBinData)
  #   }
  #   else {
  #     ind <- match(elist$genes$Array_Address_Id, rownames(tmp$Quants))
  #   }
  #   if (anyNA(ind)) {
  #     nna <- sum(is.na(ind))
  #     if (nna > tolerance) 
  #       stop("Can't match all ids in manifest with those in idat file ", 
  #            idatfiles[i], "\n", sum(is.na(ind)), " missing - please check that you have the right files, or consider setting 'tolerance'=", 
  #            sum(is.na(ind)))
  #     i <- which(!is.na(ind))
  #     ind <- ind[i]
  #     if ("MeanBinData" %in% colnames(tmp$Quants) && "DevBinData" %in% 
  #         colnames(tmp$Quants) && "NumGoodBeadsBinData" %in% 
  #         colnames(tmp$Quants)) {
  #       elist$E[i, j] <- tmp$Quants$MeanBinData[ind]
  #       elist$other$STDEV[i, j] <- tmp$Quants$DevBinData[ind]
  #       elist$other$NumBeads[i, j] <- tmp$Quants$NumGoodBeadsBinData[ind]
  #     }
  #     else {
  #       elist$E[i, j] <- tmp$Quants[ind, "Mean"]
  #       elist$other$STDEV[i, j] <- tmp$Quants[ind, "SD"]
  #       elist$other$NumBeads[i, j] <- tmp$Quants[ind, 
  #                                                "NBeads"]
  #     }
  #   }
  #   else {
  #     if ("MeanBinData" %in% colnames(tmp$Quants) && "DevBinData" %in% 
  #         colnames(tmp$Quants) && "NumGoodBeadsBinData" %in% 
  #         colnames(tmp$Quants)) {
  #       elist$E[, j] <- tmp$Quants$MeanBinData[ind]
  #       elist$other$STDEV[, j] <- tmp$Quants$DevBinData[ind]
  #       elist$other$NumBeads[, j] <- tmp$Quants$NumGoodBeadsBinData[ind]
  #     }
  #     else {
  #       elist$E[, j] <- tmp$Quants[ind, "Mean"]
  #       elist$other$STDEV[, j] <- tmp$Quants[ind, "SD"]
  #       elist$other$NumBeads[, j] <- tmp$Quants[ind, 
  #                                               "NBeads"]
  #     }
  #   }
  # }
  # if (verbose) 
  #   cat("Finished reading data.\n")
  # return(elist)
}


