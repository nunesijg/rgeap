
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
