
#' @include base.R

##########################
# AFFY
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

.initialize.affy <- function()
{
  loadpkgs('affy', 'affyio')
  affyenv <- environment(expresso)
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

# [[geapexec assign AffyReadCelFiles(string[] fileNames)]]
#' @export
read.cel.files <- function(samplefiles) 
{
  .initialize.affy()
  affyRaw <- ReadAffy(filenames = samplefiles, widget = F, verbose = F)
  colnames(exprs(affyRaw)) <- sub(pattern = ".CEL$", replacement = '', x = colnames(exprs(affyRaw)), ignore.case = T)
  rownames(affyRaw@phenoData@data) <- sub(pattern = ".CEL$", replacement = '', x = rownames(affyRaw@phenoData@data), ignore.case = T)
  rownames(affyRaw@protocolData@data) <- sub(pattern = ".CEL$", replacement = '', x = rownames(affyRaw@protocolData@data), ignore.case = T)
  return(affyRaw)
}

# [[geapexec assign AffyTreatAffyBatch(call affyBatch, string treatOption, string bgCorrect, string normMethod, string pmCorrect, string summaryMethod)]]
#' @export
treat.affy <- function(affyRaw, treatmentPackage, bgCorrect, normMethod, pmcorrect, summary)
{
  .initialize.affy()
  if (treatmentPackage == "expresso")
  {
    #library("simpleaffy")
    doBgcorrect <- (bgCorrect != "none")
    eset <- expresso(affyRaw, bg.correct = doBgcorrect, bgcorrect.method = bgCorrect,
                     normalize.method = normMethod, pmcorrect.method = pmcorrect, summary.method = summary)
    
  } else if(treatmentPackage == "gcrma")
  {
    library("gcrma")
    eset <- gcrma(affyRaw,
                  normalize = (if(normMethod != "none") TRUE else FALSE),
                  optical.correct =  (if(bgCorrect != "none") TRUE else FALSE))
    
  } else if(treatmentPackage == "plier")
  {
    library("plier")
    eset <- justPlier(affyRaw,
                      normalize = (if(normMethod != "none") TRUE else FALSE),
                      norm.type = pmcorrect)
  }
  return(eset)
}
