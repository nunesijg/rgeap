
#' @include base.R
#' @include f-esets.R
#' @include f-matrixdata.R

##########################
# GEO DATA AND MATRICES
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

.initialize.geomats <- function()
{
  loadpkgs('Biobase')
  .self.oneshot()
}

# Reads a GSEXXX_series_matrix.txt, returns an ExpressionSet
# [[geapexport void ReadGEOSeriesMatrix(string filename, bool includeMetaInfo=true) ]]
#' @export
read.series.matrix <- function(filename, include.metainfo=T)
{
  .initialize.geomats()
  con = .open.filecon(filename)
  on.exit({ if(!is.null(con)) close(con) })
  infolabel = c(hseries='!Series_', hsample='!Sample_',
                   mbegin='!series_matrix_table_begin',
                   mend='!series_matrix_table_end')
  fstruct = track.matching.lines(filename, infolabel, T)
  if (any(is.na(fstruct))) stop(sprintf("This file has no line starting with %s",
                                        paste0(names(fstruct)[is.na(fstruct)], collapse = ' or ') ))
  names(fstruct) = names(infolabel)
  fstruct = as.list(fstruct)
  # Reading the headers
  dtseries = read.table(con, sep = '\t', header = F,
                        nrows = fstruct$hsample - fstruct$hseries - 1, allowEscapes = T,
                        blank.lines.skip = TRUE, skipNul = TRUE,
                        comment.char = '', stringsAsFactors = F)
  readLines(con, 1)
  dtsmps = read.table(con, sep = '\t', header = F,
                      nrows = fstruct$mbegin - fstruct$hsample, allowEscapes = T,
                      blank.lines.skip = TRUE, skipNul = TRUE,
                      comment.char = '', stringsAsFactors = F)
  readLines(con, 1)
  dtvals = read.table(con, sep = '\t', header = 1,
                      nrows = fstruct$mend - fstruct$mbegin - 2, allowEscapes = T,
                      blank.lines.skip = TRUE, skipNul = TRUE,
                      na.strings =  c("NA", "null", "NULL", "Null"),
                      comment.char = '', stringsAsFactors = F)
  close(con)
  con = NULL
  rownames(dtvals) = make.unique(dtvals[,1])
  dtvals = dtvals[,-1]
  # Organizing the data headers
  dtseries[1,1] = sub('^ï»¿', '', dtseries[1,1])
  dtseries = .aggregate.dtinfo.paste(dtseries)
  names(dtseries) = dtseries[grep('geo_accession$', rownames(dtseries), ignore.case = T, perl = T),][1]
  rownames(dtseries) = make.unique(sub('^!Series_', '', rownames(dtseries)), sep='_')
  if (include.metainfo)
  {
    dtsmps = .aggregate.dtinfo.paste(dtsmps)
    rownames(dtsmps) = make.unique(sub('^!Sample_', '', rownames(dtsmps)), sep='_')
    colnames(dtsmps) = colnames(dtvals)
    dtsmps = as.data.frame(t(dtsmps))
  } else {
    dtsmps = data.frame(row.names=colnames(dtvals))
  }
  gplname = dtseries[grep('platform_id$', rownames(dtseries), ignore.case = T, perl = T),][1]
  fd = new("AnnotatedDataFrame", data = data.frame(row.names = rownames(dtvals)))
  pd = as(dtsmps, "AnnotatedDataFrame")
  eset = new("ExpressionSet", phenoData = pd, 
              annotation = gplname, featureData = fd, exprs = as.matrix(dtvals))
  if (include.metainfo)
  {
    notes(eset) = dtseries
  }
  return(eset)
}

# Reads a GSEXXX_series_matrix.txt, returns an ExpressionSet
# [[geapexport assign ReadGEOSeriesMatrices(string[] filenames, bool includeMetaInfo=false) ]]
#' @export
read.series.matrices <- function(filenames, include.metainfo=F)
{
  eset = NULL
  for (fname in filenames)
  {
    if (is.null(eset))
    {
      eset = read.series.matrix(fname, include.metainfo)
    } else {
      eset = esets.merge(eset, read.series.matrix(fname, include.metainfo))
    }
  }
  eset
}

# Reads a GEO file in SOFT format, returns a data.frame
# [[geapexport assign ReadSOFT(string filename, bool multifactors = false, string multfsep = null) ]]
#' @export
read.soft.data <- function(filename, multifactors=F, multfsep=' /// ')
{
  if (!file.exists(filename)) stop(sprintf("File %s not found", basename(filename)))
  .initialize.geomats()
  multifactors = multifactors && !is.null(multfsep) && !is.na(multfsep) && nchar(multfsep) != 0
  fline = .read.firstline(filename)
  type = NA_character_
  quote = ''
  if(startsWith(fline, "!Series"))
  {
    type = "series_matrix"
    quote = '"'
  } else {
    vinfo = .regex.matches(pattern = "^\\^(?<type>SERIES|SAMPLE|PLATFORM|ANNOTATION) \\= (?<id>\\w+?)$",
                           x = fline)
    if (length(vinfo) != 0)
    {
      type = tolower(vinfo['type'])
    }
  }
  nskip = 0
  nrows = -1
  ## VERIFICAR ESTRUTURA DOS DADOS AQUI
  if (!is.na(type))
  {
    vdelims = c(sprintf("!%s_table_begin", type), sprintf("!%s_table_end", type))
    vinds = track.matching.lines(filename, vdelims, prefixOnly = F, lineFeedOnly = F)
    if (!any(is.na(vinds)))
    {
      nskip = vinds[1] + 1
      nrows = (vinds[2] - vinds[1]) - 2
    } else {
      type = NA_character_
    }
  }
  if (is.na(type))
  {
    stop("This SOFT is in invalid format!")
  }
  dt = data.frame()
  (con = .open.filecon(filename)) %using% {
    dt = read.table(filename, header = T, sep = '\t', quote = quote, skip = nskip, row.names = 1,
                    nrows = nrows, blank.lines.skip = F, skipNul = F, check.names = F)
    if (multifactors)
    {
      colsel = column.classes(dt) %in% c("character", "factor")
      if (any(colsel))
      {
        for (colnm in colnames(dt)[colsel])
        {
          multf = multifactor.build.ifneeded(dt[, colnm], multfsep)
          dt = do.call('$<-', list(dt, colnm, multf))
        }
      }
    }
  }
  return(dt)
}
