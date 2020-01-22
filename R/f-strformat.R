
#' @include base.R

##########################
# STRING FORMAT METHODS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0
#

.pronom.an <- function(word)
{
  pronom = if (tolower(substring(word[1], 1, 1)) %in% c('a', 'i', 'u', 'e', 'o')) 'an' else 'a'
  pronom
}
