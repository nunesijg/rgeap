
#' @include base.R

##########################
# Quality Metrics
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

.qmetricsbase <- function(...) invisible(0)

.initialize.qmetrics <- function()
{
  loadpkgs('Biobase', 'arrayQualityMetrics')
  aqmenv = environment(arrayQualityMetrics)
  ################
  # CUSTOM HEATMAP
  aqm.heatmap2 = function (x, rgb.near = c(0,0,0), rgb.far = c(1,1,1), ...)
  {
    colorRange = rgb(seq(rgb.near[1], rgb.far[1], l = 256), seq(rgb.near[2], rgb.far[2], l = 256), 
                     seq(rgb.near[3], rgb.far[3], l = 256))
    m = genefilter::dist2(x$M)
    out = aqmenv$outliers(m, method = "sum")
    out@description = c("sum of distances to other arrays <i>S<sub>a</sub></i>", 
                        "data-driven")
    dend = as.dendrogram(hclust(as.dist(m), method = "single"))
    ord = order.dendrogram(dend)
    colnames(m) = rownames(m) = paste0(ifelse(seq_len(x$numArrays) %in% 
                                                out@which, "* ", ""), seq_len(x$numArrays))
    haveDend = (ncol(m) <= aqmenv$arrayQualityMetricsGlobalParameters$maxNumberOfArraysForDrawingDendrogram)
    if (haveDend) {
      theLegend = list(right = list(fun = latticeExtra::dendrogramGrob, args = list(x = dend, side = "right")))
      fillOrd = seq_len(x$numArrays)
    }
    else {
      theLegend = NULL
      fillOrd = ord
    }
    maxNrColors = 0
    ng = length(x$intgroup)
    if (ng > 0) {
      palettes = c("Set1", "Set2", "Set3", "Accent", "Dark2", 
                   "Paired", "Pastel1", "Pastel2")
      stopifnot(all(palettes %in% rownames(RColorBrewer::brewer.pal.info)))
      palettes = rep(palettes, ceiling(ng/length(palettes)))
      key = rects = vector(mode = "list", length = ng)
      names(rects) = rep("rect", ng)
      for (i in seq_len(ng)) {
        colors = RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[palettes[i], 
                                            "maxcolors"], palettes[i])
        fac = factor(x$pData[[x$intgroup[i]]])
        fac = aqmenv$maximumLevels(fac, n = length(colors))
        colors = colors[seq_len(nlevels(fac))]
        key[[i]] = list(rect = list(col = colors), text = list(levels(fac)))
        rects[[i]] = list(col = "transparent", fill = colors[as.integer(fac)[fillOrd]], 
                          type = "rectangle")
        if (length(colors) > maxNrColors) 
          maxNrColors = length(colors)
      }
      key = unlist(key, recursive = FALSE)
      key$rep = FALSE
      thekey = lattice::draw.key(key = key)
      if (haveDend) {
        theLegend$right$args = append(theLegend$right$args, 
                                      list(size.add = 1, add = rects))
      }
      else {
        lay = grid::grid.layout(nrow = 1, ncol = ng, heights = unit(1, 
                                                              "null"), widths = unit(rep(1, length = ng), 
                                                                                     rep("lines", ng)), respect = FALSE)
        g = grid::frameGrob(layout = lay)
        dy = 1/x$numArrays
        y = seq_len(x$numArrays) * dy
        for (i in seq_len(ng)) {
          g = grid::placeGrob(g, grid::rectGrob(y = y, height = dy, 
                                    vjust = 1, gp = do.call(grid::gpar, rects[[i]])), 
                        row = 1, col = i)
        }
        idem = function(x) x
        theLegend = list(right = list(fun = idem, args = list(x = g)))
      }
    }
    else {
      thekey = NULL
    }
    hfig = lattice::levelplot(m[ord, ord], scales = list(x = list(rot = 90)), 
                     legend = theLegend, colorkey = list(space = "left"), 
                     xlab = "", ylab = "", col.regions = colorRange, main = thekey)
    nout = length(out@which)
    legend = paste0("Heatmap false color.")
    new("aqmReportModule", plot = hfig, section = "Between array comparison", 
        title = "Distances between arrays", id = "hm", legend = legend, 
        size = c(w = 5 + x$numArrays * 0.075, h = 3 + x$numArrays * 
                   0.075 + maxNrColors * 0.2), colors = x$arrayColors, 
        outliers = out)
  }
  ##############
  # QMETRICSBASE
  qmetricsf = function (expressionset, outdir,
                              do.logtransform = FALSE, intgroup = character(0), 
                              reporttitle = "Quality Metrics Results", ...) 
  {
    if (!dir.exists(outdir)) dir.create(outdir, showWarnings = F, recursive = T)
    m = list()
    old.seed = setRNG::setRNG(kind = "default", seed = 28051968, 
                      normal.kind = "default")
    on.exit({setRNG::setRNG(old.seed); grDevices::graphics.off()})
    if (is.matrix(expressionset)) expressionset = df.to.eset(expressionset)
    x = aqmenv$prepdata(expressionset, intgroup = intgroup, do.logtransform = do.logtransform)
    .give.status(percent=16, message = "Gerando HeatMaps...", engMsg = "Building HeatMaps...")
    m$heatmap = aqm.heatmap2(x, ...)
    m$heatmap_out = aqmenv$aqm.outliers(m$heatmap)
    .give.status(percent=32, message = "Gerando PCA...", engMsg = "Building PCA...")
    m$pca = aqmenv$aqm.pca(x, ...)
    .give.status(percent=46, message = "Gerando Box Plot...", engMsg = "Building BoxPlot...")
    m$boxplot = aqmenv$aqm.boxplot(x, ...)
    m$boxplot_out = aqmenv$aqm.outliers(m$boxplot)
    .give.status(percent=54, message = "Calculando estimativas de densidade...", engMsg = "Computing density estimates...")
    m$density = aqmenv$aqm.density(x, ...)
    .give.status(percent=68, message = "Calculando desvios-padrão de intensidade...", engMsg = "Computing intensity standard deviations...")
    m$meansd = aqmenv$aqm.meansd(x, ...)
    .give.status(percent=80, message = "Calculando distribuição de intensidades...", engMsg = "Computing intensity distributions...")
    m$probesmap = aqmenv$aqm.probesmap(x, ...)
    .give.status(percent=90, message = "Gerando dados MA...", engMsg = "Building MA plot data...")
    m$maplot = aqmenv$aqm.maplot(x, ...)
    m$maplot_out = aqmenv$aqm.outliers(m$maplot)
    .give.status(percent=96, message = "Gerando gráficos espaciais...", engMsg = "Building spatial plots...")
    m = append(m, aqmenv$aqm.spatial(x, ...))
    .give.status(percent=99, message = "Renderizando gráficos...", engMsg = "Rendering plots...")
    dpi = aqmenv$arrayQualityMetricsGlobalParameters$dpi
    grDevices::graphics.off()
    smpnames = colnames(x$M)
    outdf = data.frame(row.names=smpnames)
    outimgfiles = character(0)
    svgMode = !is.windows() || getOption('force.svg', T)
    ext = if (svgMode) ".svg" else ".wmf"
    for(i in 1:length(m))
    {
      module = m[[i]]
      if (is.null(module)) next
      mname = module@id
      outimgf = sub('[\\\\/]{2,}', '/', paste0(outdir, '/', mname, ext))
      outimgfiles[mname] = outimgf
      h = module@size['h'] * dpi / 75
      w = module@size['w'] * dpi / 75
      suppressWarnings(expr = {
        if (is.windows() && !getOption('force.svg', F))
        {
          win.metafile(filename = outimgf, height = h, width = w, family = "sans")
        } else {
          svg(filename = outimgf, height = h, width = w, family = "sans")
        }
        aqmenv$makePlot(module)
      })
      grDevices::graphics.off()
      #while (!is.null(dev.list())) dev.off()
      modstats = module@outliers@statistic
      if (length(modstats) != 0)
      {
        modnames = names(modstats)
        if (is.null(modnames)) modnames = smpnames
        outdf[modnames, mname] = ifelse(modstats > module@outliers@threshold, 'x', '')
      }
    }
    return(list(modules=m, outlierdf=outdf, imgfiles=outimgfiles))
  }
  .parent.assign('.qmetricsbase', qmetricsf)
  .self.oneshot()
}

# Analyzes a expression set
# Returns a list with slots:
# - modules : list of 'aqmReportModule' objects
# - outlierdf : data.frame of outlier, with the found outliers marked with 'x' (columns are hm, box and ma, rows are sample names)
# - imgfiles : named character vector with output files (hm, out hm, pca, box, out box, dens, msd, ma, out ma)
# [[geapexec robj_RList QMetricsAnalyze(call eset, string outdir, bool doLog, string[] intgroup, Color rgbNear, Color rgbFar)]]
#' @export
qmetrics.analyze <- function(eset, outdir, do.log, intgroup=NULL, rgb.near=NULL, rgb.far=NULL)
{
  .initialize.qmetrics()
  if (is.null(rgb.near)) rgb.near = c(1,1,0)
  if (is.null(rgb.far)) rgb.far = c(0,0,1)
  aqmres = .qmetricsbase(eset, outdir, do.log=do.log, intgroup=intgroup, rgb.near=rgb.near, rgb.far=rgb.far)
  aqmres$modules = NULL
  aqmres
}


