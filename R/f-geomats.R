
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
# [[geapexport assign ReadGEOSeriesMatrix(path filename, bool includeMetaInfo=true) ]]
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
  rownames(dtvals) = make.unique(as.character(dtvals[,1]))
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
# [[geapexport assign ReadGEOSeriesMatrices(path[] filenames, bool includeMetaInfo=false) ]]
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
# [[geapexport assign ReadSOFT(path filename, bool multifactors = false, string multfsep = null) ]]
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
  # Below the file structure must be checked
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
    dt = read.table(con, header = T, sep = '\t', quote = quote, skip = nskip,
                    nrows = nrows, blank.lines.skip = F, skipNul = F, check.names = F)
    if (anyNA(dt[,1L]))
    {
      dt = dt[-which(is.na(dt[,1L])),]
    }
    row.names(dt) = dt[,1L]
    dt = dt[,-1L, drop=FALSE]
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
  colnames(dt) = make.unique(colnames(dt), sep = '_')
  return(dt)
}

# Reads a GEO file in SOFT format, returns a data.frame. Returns a data.frame with the table contents
# [[geapexport assign ReadMultiSOFT(path filename) ]]
#' @export
read.multisoft.file <- function(filename)
{
  dtindexes = read.multisoft.table.binpos(filename)
  eset = NULL
  for (r in rownames(dtindexes))
  {
    lcount = dtindexes[r, 'linecount']
    if (is.na(lcount) || lcount <= 0) next
    ttype = dtindexes[r,'type']
    tbindex = dtindexes[r, 'tablebindex']
    dt = read.table.atpos(filename, tbindex, lcount)
    if (ttype %in% 'sample')
    {
      dt = dt[,1L,drop=FALSE]
      colnames(dt) = r
    }
    if (is.null(eset))
    {
      eset = df.to.eset(dt)
    }
    else
    {
      eset = esets.merge(eset, dt)
    }
  }
  eset
}

read.multisoft.table.binpos <- function(filename)
{
  vheaderprefs = c('SERIES', 'SAMPLE', 'PLATFORM', 'ANNOTATION')
  vheaderpatt = sprintf("^%s = ", vheaderprefs)
  indBin = track.matching.string.binpos(filename, vheaderpatt)
  if (is.na(indBin)) stop("Invalid SOFT format")
  offsets = list()
  headerregpatt = paste0(sprintf("(?:%s)", sub('\\^', '\\\\^', vheaderpatt)), collapse = '|')
  vgeotypepatt = sprintf("\\^(%s) = .+", paste0(vheaderprefs, collapse = '|'))
  vgeotypes = character(0)
  while (!is.na(indBin))
  {
    geoname = NULL
    (con = .open.filecon(filename)) %using% {
      seek(con, indBin)
      sline = readLines(con, n=1L)
      geoname = sub(headerregpatt, '', sline, perl = T)
      vgeotypes[geoname] = tolower(sub(vgeotypepatt, '\\1', sline, perl=T))
    }
    offsets[geoname] = indBin
    indBin = track.matching.string.binpos(filename, vheaderpatt, binOffset = indBin + 1)
  }
  offsets = data.frame(row.names = names(offsets), hindex=unlist(offsets))
  offsets$size = c(offsets$hindex[-1], file.size(filename)) - offsets$hindex
  offsets$type = vgeotypes
  offsets$tablebindex = rep(NA_real_, nrow(offsets))
  offsets$linecount = rep(NA_integer_, nrow(offsets))
  
  for (i in 1:nrow(offsets))
  {
    hindex = offsets[i, 'hindex']
    maxbin = offsets[i, 'size']
    dttype = offsets[i, 'type']
    tabpattbegin = sprintf("!%s_table_begin", dttype)
    tabpattend = sprintf("!%s_table_end", dttype)
    tbindex = track.matching.string.binpos(filename, tabpattbegin, binOffset = hindex, maxBinRead = maxbin)
    tlinecount = diff(track.matching.lines(filename, c(tabpattbegin, tabpattend), binOffset = tbindex, maxBinRead = maxbin)) - 2
    if (!is.na(tlinecount))
    {
      offsets[i, 'linecount'] = tlinecount
      tbindex = track.next.line.binpos(filename, tbindex)
      offsets[i, 'tablebindex'] = tbindex
    }
  }
  
  offsets
}

# Summarizes the ID column from a table, removing or merging duplicated values based on duplicates from other columns
# [[geapexport assign TableSummarizeRowNamesByCol(call table, string idCol, string[] dupRefCols, string dupSep) ]]
#' @export
table.summarize.rownames.bycol <- function(dt, idcol, duprefcols=NULL, dupsep=' /// ')
{
  if (is.character(idcol))
  {
    if (!(idcol %in% colnames(dt))) stop(sprintf("Missing ID column: %s", idcol))
    idcol = which(idcol == colnames(dt))[1]
  }
  vcolid = dt[, idcol]
  if (anyNA(vcolid))
  {
    # Removes the rows where the ID column is NA
    dt = dt[!is.na(vcolid), ]
    vcolid = dt[, idcol]
  }
  if (anyDuplicated(vcolid))
  {
    seldups = duplicated(vcolid)
    duprefcols = as.vector(na.omit(duprefcols))
    if (!is.null(duprefcols) && length(duprefcols != 0))
    {
      # Finds and removes the values which were also duplicated in reference columns
      duprefcols = duprefcols[duprefcols %in% colnames(dt)]
      for (rfcol in duprefcols)
      {
        vcolref = dt[,rfcol]
        selrefdups = duplicated(vcolref)
        seldups = seldups & selrefdups
      }
      dt = dt[!seldups,,drop=FALSE]
      vcolid = dt[, idcol]
      seldups = duplicated(vcolid)
    }
    if (any(seldups))
    {
      # Aggregates the remaining duplicated values
      vuid = as.character(unique(vcolid))
      ind.dups.id = match(vcolid, vcolid, nomatch = -1L)
      facs.id = factor(vcolid, levels=vuid)
      df.res = data.frame(row.names = vuid)
      for (colnm in colnames(dt)[-idcol])
      {
        vcol = dt[,colnm]
        ind.col = match(vcol, vcol)
        lsuniques = tapply(ind.col, facs.id, unique, simplify = FALSE)
        df.res[,colnm] = sapply(lsuniques, function(e) { paste0(vcol[e], collapse=dupsep) })
      }
      if (is.matrix(dt))
        df.res = as.matrix(df.res)
      dt = df.res
      return(dt)
    }
  }
  rownames(dt) = dt[, idcol]
  dt = dt[, -idcol]
  dt
}

