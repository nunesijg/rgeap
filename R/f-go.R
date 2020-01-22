
#' @include base.R
#' @include f-fileio.R

##########################
# GENE ONTOLOGY FUNCTIONS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

.initialize.go <- function()
{
  library("stats")
  library("data.table")
  library("topGO")
  library("devtools")
  .self.oneshot()
}

# Reads a GAF file
# [[geapexport assign ReadGAF(string fileName)]]
#' @export
read.gaf <- function(fnamegaf)
{
  dtgaf = read.table(fnamegaf, header = F, quote = '', comment.char = '!', sep = '\t', colClasses = 'character', blank.lines.skip = T, skipNul = T)
  dtgaf = data.frame(accession=dtgaf[,2], symbol=dtgaf[,3], go=dtgaf[,5])
  dtgaf
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
  if (is.null(allGenes)) allGenes = unique(unlist(go2genelist))
  allGenes = as.character(allGenes)
  selGenes = as.character(selGenes)
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
  resTest = runTest(godata, algorithm = algType, statistic = testType)
  nSigGenes = sum(score(resTest) < pval.cutoff)
  gotable = GenTable(godata, `P.Value`=resTest, topNodes = nSigGenes, orderBy="P.Value" )
  .give.status(message="Inserindo genes na tabela", engMsg="Appending genes to table")
  gotable = go.append.result.genes(godata, gotable)
  if (!is.null(plotOutFile) && is.character(plotOutFile) && nchar(plotOutFile) > 1  && dir.exists(dirname(plotOutFile)))
  {
    .give.status(message="Criando gráficos dos resultados", engMsg="Plotting results")
    win.metafile(plotOutFile, height = 10, width=10, family='sans')
    showSigOfNodes(godata, score(resTest), firstSigNodes = 5)
    while (!is.null(dev.list())) dev.off()
    closeAllConnections()
  }
  gotable
}

# Appends the gene names to GO results table. The new column has results like GENE1|GENE2|GENE3
#go.append.result.genes <- function(gotable, go2genelist, selGenes)
#{
#  filtgos = go2genelist[names(go2genelist) %in% gotable$GO.ID]
#  filtgos = lapply(filtgos, function(v) v[v %in% selGenes])
#  filtgos = filtgos[sapply(filtgos, length) != 0]
#  gosgenes = sapply(filtgos, paste0, collapse='|')
#  gotable$genes = rep("", nrow(gotable))
#  gotable$genes[na.omit(match(names(gosgenes), gotable$GO.ID))] = gosgenes
#  gotable
#}

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
# [[geapexport assign Org2GeneList(string pkgPath, string pkgName, string geneids, string ontology)]]
#' @export
org2genelist <- function(pkgpath=NULL, pkgname=NULL, geneids='symbol', ontology='BP')
{
  .initialize.go()
  loaded = F
  if (is.character(pkgpath) && dir.exists(pkgpath))
  {
    pkgenv = load_all(pkgpath)
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
  on.exit(detach(sprintf('package:%s', pkgname), unload = T, character.only = T))
  go2genels = annFUN.org(ontology, mapping = pkgname, ID = geneids)
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
