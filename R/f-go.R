
#' @include base.R
#' @include f-fileio.R

##########################
# GENE ONTOLOGY FUNCTIONS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

.initialize.go <- function()
{
  loadpkgs('stats', 'data.table', 'devtools', 'topGO', 'DBI')
  .self.oneshot()
}

# Reads a data.frame containing the GO annotation for an organism
read.go.annot.file <- function(fname, accession.col, symbol.col, go.col, ...)
{
  args.read.table = list(header=FALSE,
                         quote='',
                         comment.char='!',
                         sep='\t',
                         colClasses='character',
                         blank.lines.skip=TRUE,
                         skipNul=TRUE)
  argls = list(...)
  for (rtarg in intersect(formalArgs(read.table), names(argls)))
    args.read.table[[rtarg]] = argls[[rtarg]]
  args.read.table$file = fname
  dt = do.call(read.table, args.read.table)
  indexes = list(
    accession=accession.col,
    symbol=symbol.col,
    go=go.col
  )
  dt = data.frame(accession=dt[,indexes$accession, drop=TRUE],
                  symbol=dt[,indexes$symbol,drop=TRUE],
                  go=dt[,indexes$go,drop=TRUE])
  dt
}

# Reads a GAF file. TAIR files are also supported
# [[geapexport assign ReadGAF(path fileName)]]
#' @export
read.gaf <- function(fname)
{
  dt = read.go.annot.file(fname,
                          accession.col = 2L,
                          symbol.col = 3L,
                          go.col = 5L)
  dt
}


# Performs a straightforwad gene ontology analysis, returning a data.frame with the most significant gene ontologies
# [[geapexec assign GOAnalysisCustom(call go2genelist, call selGenes, call allGenes, double pValCutoff, int nodeSize, string testType, string algType, string pAdjustMethod, string ontology, string plotOutFile=null)]]
#' @export
go.analysis.custom <- function(go2genelist, selGenes, allGenes=NULL, pval.cutoff=0.05, nodeSize=5, testType='fisher', algType='classic', pAdjust='fdr', ontology='BP', plotOutFile='')
{
  .initialize.go()
  testType = testType[1]
  algType = algType[1]
  if (is.logical(pAdjust))
  {
    pAdjust = ifelse(pAdjust, 'fdr', 'none')
  } else if (!(is.character(pAdjust) && (pAdjust[1] %in% stats::p.adjust.methods)))
  {
    stop("pAdjust must be a single character parameter. See options in p.adjust.methods")
  }
  if (!(testType %in% whichTests())) stop("testType must be a valid character option from whichTests()")
  if (!(algType %in% whichAlgorithms())) stop("algType must be a valid character option from whichAlgorithms()")
  if (is.null(allGenes)) allGenes = unique(na.exclude(unlist(go2genelist)))
  allGenes = as.character(allGenes)
  selGenes = as.character(na.exclude(selGenes))
  nmatchgens = sum(na.exclude(selGenes %in% allGenes))
  if (nmatchgens == 0L)
    stop("None of the selected genes has been found in the reference list. ",
         "Another column should be selected as the gene/probe identifier.")
  geneList = factor(as.integer(allGenes %in% selGenes))
  names(geneList) = allGenes
  # GARGALO 1
  godata = new("topGOdata", ontology=ontology,
               allGenes=geneList,
               nodeSize = nodeSize,
               annot = annFUN.GO2genes,
               GO2genes = go2genelist,
               description="Gene ontology results")
  # GARGALO 2
  .give.status(message=sprintf("Rodando teste (%s) com %s dentre %s genes", testType, nmatchgens, length(allGenes)),
               engMsg=sprintf("Running test (%s) with %s of %s genes", testType, nmatchgens, length(allGenes)))
  resTest = topGO::runTest(godata, algorithm = algType, statistic = testType)
  nSigGenes = sum(score(resTest) < pval.cutoff)
  gotable = topGO::GenTable(godata, `P.Value`=resTest, topNodes = nSigGenes, orderBy="P.Value" )
  .give.status(message="Inserindo genes na tabela", engMsg="Appending genes to table")
  gotable = go.append.result.genes(godata, gotable)
  if (!is.null(plotOutFile) && is.character(plotOutFile) && nchar(plotOutFile) > 1  && dir.exists(dirname(plotOutFile)))
  {
    .give.status(message="Criando gráficos dos resultados", engMsg="Plotting results")
    gr.vector.create(plotOutFile, height=10, width=10, pointsize = 16)
    on.exit({grDevices::graphics.off()})
    showSigOfNodes(godata, score(resTest), firstSigNodes = 5)
    while (!is.null(dev.list())) dev.off()
    closeAllConnections()
  }
  gotable
}

# Appends the gene names to GO results table. The new column (list of characters) is appended to the table
# [[geapexport assign GOAppendResultGenes(call godata, call gotable, string[] terms, bool significantOnly, int geneCutoff)]]
#' @export
go.append.result.genes <- function(godata, gotable, whichTerms=NULL, significantOnly=T, geneCutoff=50)
{
  .initialize.go()
  if (length(whichTerms) == 0)
  {
    whichTerms = gotable$GO.ID
  } else {
    gotable = gotable[gotable$GO.ID %in% whichTerms,]
  }
  genls = topGO::genesInTerm(godata, whichTerms)
  retls = sapply(gotable$GO.ID, function(x) character(0), simplify = F)
  
  for (goterm in whichTerms)
  {
    termgenes = genls[[goterm]]
    pval = sort(topGO::geneScore(godata, termgenes, use.names = TRUE))
    if (significantOnly)
    {
      termgenes = names(pval[pval != 1])
    } else {
      termgenes = names(pval)
    }
    length(termgenes) = min(length(termgenes), geneCutoff)
    if (length(termgenes) != 0)
    {
      retls[[goterm]] = termgenes
    }
  }
  class(retls) = c('list', 'multifactor')
  gotable$Genes = retls
  gotable
}

# Produces a list of gene ontologies to genes based on a annotation database
# If pkgpath is provided and is a valid package directory path, pkgname is ignored
# geneids options are 'entrez', 'genbank', 'alias', 'ensembl', 'symbol', 'genename', 'unigene'
# Ontology can be 'BP', 'MF' or 'CC'
# [[geapexport assign Org2GeneList(path pkgPath, string pkgName, string geneids, string ontology)]]
#' @export
org2genelist <- function(pkgpath=NULL, pkgname=NULL, geneids='symbol', ontology='BP')
{
  .initialize.go()
  loaded = F
  if (is.character(pkgpath) && nchar(pkgpath) != 0L && dir.exists(pkgpath))
  {
    pkgenv = devtools::load_all(pkgpath)
    pkgname = pkgenv$env$.packageName
    loaded = T
  } else if (is.character(pkgname))
  {
    if (pkgname %in% loaded.packages())
    {
      loaded = T
    } else if (pkgname %in% all.packages())
    {
      library(pkgname, character.only = T)  
      loaded = T
    }
  }
  if (!loaded) stop("Expected either a valid directory from pkgpath or installed package name in pkgname")
  on.exit({detach(sprintf('package:%s', pkgname), unload = T, character.only = T)})
  defargs = c("entrez", "genbank", "alias", "ensembl",
              "symbol", "genename", "unigene", "uniprot", "refseq")
  
  tableName = c("genes", "accessions", "alias", "ensembl", 
                "gene_info", "gene_info", "unigene", "uniprot", "refseq")
  
  keyName = c("gene_id", "accessions", "alias_symbol", "ensembl_id", 
              "symbol", "gene_name", "unigene_id", "uniprot_id", "accession")
  names(tableName) = names(keyName) = defargs
  mapping = paste(sub(".db$", "", pkgname), ".db", sep = "")
  connfname = paste(sub(".db$", "", mapping), "dbconn", sep = '_')
  if (!exists(connfname))
    stop("Package '", pkgname, "' was loaded, but the underlying functions are not available. ",
         "This problem may occur due to a compilation error.\n",
         "Hence, the required package may need to be locally compiled from the source. ",
         "If this is a temporary installation, one can also accept to install it in ",
         "the local library.")
  # Opens the connection to local package database
  con = get(connfname)()
  .give.status(message="Mapeando ontologias do organismo", engMsg="Mapping the organism ontologies")
  tablenms = dbGetQuery(con, "SELECT name FROM sqlite_master WHERE type='table'")[,1L]
  geneids_key = tolower(geneids)
  genetable = if (geneids_key %in% names(tableName)) tableName[geneids_key] else geneids_key
  geneID = keyName[geneids_key]
  if (!(genetable %in% tablenms))
  {
    alternatives = switch(tolower(geneids),
                          alias = c(gene2alias = 'alias'),
                          symbol = c(gene2systematic = 'gene_name'),
                          genename = c(gene2systematic = 'systematic_name'),
                          NULL)
    if (!is.null(alternatives))
    {
      for (tbalt in names(alternatives))
      {
        if (!(tbalt %in% tablenms)) next
        tbalt.colnms = dbGetQuery(con, sprintf("PRAGMA table_info('%s')", tbalt))$name
        colalt = alternatives[tbalt]
        if (!(colalt %in% tbalt.colnms)) next
        genetable = tbalt
        geneID = colalt
      }
    }
    if (!(genetable %in% tablenms))
      stop(paste(sprintf("This organism database (%s package) has no table named '%s', required for '%s' entries.",
                         pkgname, genetable, geneids),
                         "Another identifier should be used instead."))
  }
  if (is.na(geneID))
  {
    tb.colnms = dbGetQuery(con, sprintf("PRAGMA table_info('%s')", genetable))$name
    geneID = setdiff(tb.colnms, "_id")[1]
  }
  tb.onto = sprintf("go_%s", tolower(ontology))
  if (!(tb.onto %in% tablenms))
    stop(sprintf("This organism database (%s) has no %s ontology data", pkgname, ontology))
  
  qsql = paste(sprintf("SELECT DISTINCT %s, go_id FROM", geneID),
               sprintf("%s INNER JOIN %s USING (_id)", genetable, tb.onto))
  go2genels = dbGetQuery(con, qsql)
  go2genels = split(go2genels[[geneID]], go2genels[['go_id']])
  go2genels
}

# Producest a list of gene ontologies to genes based on a GAF file
# geneids options are 'accession' and 'symbol'
# [[geapexport assign Gaf2GeneList(string filePath, string geneids)]]
#' @export
gaf2genelist <- function(filepath, geneids='symbol')
{
  dtgaf = read.gaf(filepath)
  if (!(is.character(geneids) && geneids %in% c('symbol', 'accession')))
  {
    stop("geneids must be a character ('accession' or 'symbol')")
  }
  dtgaf = data.frame(go=dtgaf$go, gene=dtgaf[,geneids[1]])
  go2genels = tapply(as.character(dtgaf$gene), dtgaf$go, c, simplify = F)
  go2genels = lapply(go2genels, c)
  go2genels
}
