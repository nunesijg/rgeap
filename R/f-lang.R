
#' @include base.R

##########################
# Language Functions
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0

# [[geapexport void SetTextLanguage(string twoLetterLang)]]
#' @export
set.geap.text.language <- function(lang)
{
  lang = match.arg(lang, c("en", "pt"))
  set.geap.option('geap.text.lang', lang)
  invisible(lang)
}
