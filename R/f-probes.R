
#' @include base.R

##########################
# PROBE PROCESSING METHODS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0
#

.initialize.probes <- function()
{
  loadpkgs('devtools', 'Biobase')
  .self.oneshot()
}

.initialize.probes.affy <- function()
{
  loadpkgs('devtools', 'Biobase', 'affy')
  .self.oneshot()
}

.initialize.probes.oligo <- function()
{
  loadpkgs('devtools', 'Biobase', 'oligo')
  .self.oneshot()
}

# Loads a local package
# [[geapinternal string LoadLocalPackage(path pkgDir)]]
#' @export
loadLocalPackage <- function(pkgdir)
{
  .initialize.probes()
  pkgstatus = load_all(path = pkgdir, recompile = F, reset = T, export_all = F, quiet = T)
  return(pkgstatus$env$.packageName)
}

# Unloads a package, returing TRUE if it was previously loaded or FALSE if it's already detached
# [[geapinternal bool UnloadPackage(string pkgName)]]
#' @export
unloadPackage <- function(pkgName)
{
  .initialize.probes()
  pkgnmsp = sprintf("package:%s", pkgName)
  if (pkgnmsp %in% search())
  {
    detach(pkgnmsp, character.only = T, unload = T, force = T)
    return(T)
  }
  return(F)
}

# Gets the internal platform name from a CEL file
# [[geapexport string GetAffyPlatformName(path celFile)]]
#' @export
getAffyPlatformName <- function(celfile)
{
  .initialize.probes.affy()
  plat = affy::whatcdf(celfile)
  plat = tolower(gsub("[\\s\\_\\-]", "", plat, perl = T))
  plat
}

# Gets the CDF package name from a CEL file
# [[geapexport string GetCDFPkgName(path celFile)]]
#' @export
getCdfPkgName <- function(celfile)
{
  .initialize.probes.affy()
  plat = getAffyPlatformName(celfile)
  if (!endsWith(plat, "cdf")) plat = sprintf("%scdf", plat)
  plat
}

# Gets the first probes from a CDF package name (which must be loaded beforehand)
# [[geapexport string[] GetCDFProbeNames(string pkgname, int limit=0)]]
#' @export
getCdfProbeNames <- function(pkgname, amount=0)
{
  if (!exists(pkgname)) stop(sprintf("The required package ('%s') does not exist or is not currently loaded", pkgname))
  if (amount > 0)
  {
    return(head(names(get(pkgname)), n = amount))
  }
  return(names(get(pkgname)))
}

# Gets the first probes from a 'probe' package name (which must be loaded beforehand)
# [[geapexport string[] GetAnnProbeNames(string pkgname, int limit=0)]]
#' @export
getAnnProbeNames <- function(pkgname, amount=0)
{
  if (!exists(pkgname)) stop(sprintf("The required package ('%s') does not exist or is not currently loaded", pkgname))
  if (amount > 0)
  {
    return(head(unique(get(pkgname)$'Probe.Set.Name'), n = amount))
  }
  return(unique(get(pkgname)$'Probe.Set.Name'))
}

# Gets the first probes from a 'pd' package name (which must be loaded beforehand)
# [[geapexport string[] GetPdProbeNames(string pkgname, int limit=0)]]
#' @export
getPdProbeNames <- function(pkgname, amount=0)
{
  if (!exists(pkgname)) stop(sprintf("The required package ('%s') does not exist or is not currently loaded", pkgname))
  pkg = get(pkgname)
  sql = "select distinct man_fsetid from featureSet order by featureSet.fsetid ASC"
  if (amount > 0)
    sql = sprintf("%s limit %d", sql, amount)
  conn = oligo::db(pkg)
  on.exit({dbDisconnect(conn)})
  res = as.character(dbGetQuery(conn, sql)[,1])
  res
}


