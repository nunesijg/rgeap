
#' @include base.R

##########################
# GRAPHICS FUNCTIONS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

.initialize.graphics <- function()
{
  loadpkgs('Cairo', 'png')
  .self.oneshot()
}

gr.plot.raw.png <- function(expr)
{
  .initialize.graphics()
  toeval = substitute(expr)
  on.exit(graphics.off())
  Cairo::Cairo(file='/dev/null')
  base::eval.parent(expr = expr)
  im = Cairo:::.image(dev.cur())
  r = Cairo:::.ptr.to.raw(im$ref, 0, im$width * im$height * 4)
  dim(r) = c(4, im$width, im$height)
  r[c(1,3),,] = r[c(3,1),,]
  p = writePNG(r, raw())
  p
}

# Gets the extension used in vector images
gr.vector.ext <- function(fnames=character(0L))
{
  svgMode = !is.windows() || getOption('force.svg', T)
  ext = if (svgMode) ".svg" else ".wmf"
  if (length(fnames) == 0L) return(ext)
  sub('[\\\\/]{2,}', '/', paste0(fnames, ext))
}

# Creates a vector image using the current configuration (either SVG or WMF)
gr.vector.create <- function(filename, width = 7, height = 7, pointsize=12, family="sans")
{
  svgMode = !endsWith(tolower(filename), ".wmf") &&
    (endsWith(tolower(filename), ".svg") || !is.windows() || getOption('force.svg', TRUE))
  if (svgMode)
  {
    svg(filename = filename, height = height, width = width, pointsize=pointsize, family = family)
  } else {
    win.metafile(filename = filename, height = height, width = width, pointsize=pointsize, family = family)
  }
}
