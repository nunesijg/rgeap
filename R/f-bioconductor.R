
#' @include base.R

##########################
# INPUT/OUTPUT FUNCTIONS
# ########################
# Nunes et al, 2019
# Last updated version: 0.3.0


.initialize.bioconductor <- function()
{
  loadpkgs('BiocManager')
  repos = getOption('repos')
  if ('CRAN' %in% names(repos) && repos[['CRAN']] == '@CRAN@') repos['CRAN'] = "https://cloud.r-project.org/"
  repos = c(repos, BiocManager::repositories())
  repos = repos[unique(names(repos))]
  options(repos = repos)
  .self.oneshot()
}

# [[geapexec void InstallPackages(params string[] pkgList)]]
# Installs packages from available repositories (including CRAN and Bioconductor)
#' @export
install.pkgs <- function(pkglist)
{
  .initialize.bioconductor()
  i = 0
  tot = length(pkglist)
  pkgpaths = utils::installed.packages(fields = 'LibPath')[,'LibPath']
  for (pkgnm in pkglist)
  {
    instdir = .libPaths()[1]
    pkglocalname = pkgnm
    repos = getOption("repos")
    if (grepl("[\\\\/]", pkgnm) && file.exists(pkgnm))
    {
      pkglocalname = gsub("^([^_]+?)_?(?:\\d+?\\.)+?tar\\.gz", '\\1', basename(pkgnm), perl=TRUE)
      repos = NULL
    } else if (pkgnm %in% names(pkgpaths))
    {
      instdir = pkgpaths[[pkgnm]]
    }
    .give.status(percent=as.integer(100 * i / tot),
                 message=sprintf("Instalando pacote: %s", pkglocalname),
                 engMsg=sprintf("Installing package: %s", pkglocalname))
    utils::install.packages(pkgnm, repos=repos, quiet = TRUE, verbose = FALSE, lib = instdir, type='both')
    if (!(pkglocalname %in% getAllInstalledPackages()))
    {
      .give.status(percent=as.integer(100 * i / tot),
                   message=sprintf("Instalando pacote (2ª tentativa): %s", pkglocalname),
                   engMsg=sprintf("Installing package (2nd attempt): %s", pkglocalname))
      utils::install.packages(pkglocalname, repos=repos, quiet = TRUE, verbose = FALSE, lib = instdir, type='both')
    }
    i = i + 1
  }
  invisible(0)
}

# [[geapexport robj_RDataFrame OldPackagesTable()]]
# Gets a data.frame of packages to be updated
#' @export
get.old.pkgs <- function()
{
  .initialize.bioconductor()
  dt = as.data.frame(utils::old.packages())
  if (ncol(dt) == 0) dt = data.frame(Package=NA_character_, LibPath=NA_character_, Repository=NA_character_)[-1,,drop=F]
  dt
}

# [[geapexport string GetBiocManagerVersion()]]
# Gets the current version of BiocManager as string
#' @export
get.biocmanager.version <- function()
{
  .initialize.bioconductor()
  vers = as.character(BiocManager::version())
  vers
}

# [[geapexport string[] GetInstalledPackages()]]
# Gets a character vector listing all currently installed packages
#' @export
getAllInstalledPackages <- function()
{
  pkgs = names(utils::installed.packages(fields = 'Package', noCache = TRUE)[,'Package'])
  pkgs
}

# [[geapexport string[] GetLoadedPackages(bool includeBase = false)]]
# Gets a character vector listing all currently loaded packages
#' @export
getLoadedPackages <- function(includeBase=FALSE)
{
  if (includeBase)
    pkgs = unique(sub("^package:", "", grep("^package:", search(), value = T, perl = T), perl = T))
  else
    pkgs = names(sessionInfo()$otherPkgs)
  pkgs
}

# [[geapexport int UnloadLocalLoadedPackages()]]
# Unloads only the locally loaded packages (e.g., by devtools). Returns the number of detached packages
#' @export
unloadLocalLoadedPackages <- function()
{
  nsnames = names(sessionInfo()$otherPkgs)
  localpkgs = nsnames[vapply(nsnames,
                             function(ns) '.__DEVTOOLS__' %in% names(getNamespace(ns)), FALSE, USE.NAMES = TRUE)]
  localpkgs = setdiff(localpkgs, packageName())
  nUnloaded = 0L
  for (pkgnm in localpkgs)
  {
    pkgns = sprintf('package:%s', pkgnm)
    success = FALSE
    tryCatch({
      detach(pkgns, unload = TRUE, character.only = TRUE, force = TRUE)
      success = TRUE
    }, error=function(e) warning(e))
    if (success)
      nUnloaded = nUnloaded + 1L
  }
  nUnloaded
}

# [[geapexport int UnloadAllLoadedPackages()]]
# Unloads the previously loaded packages except rgeap and base packages. Returns the number of detached packages
#' @export
unloadAllLoadedPackages <- function()
{
  pkgs = setdiff(getLoadedPackages(), packageName())
  nUnloaded = 0L
  for (pkg in pkgs)
  {
    pkgns = sprintf('package:%s', pkg)
    success = FALSE
    tryCatch({
      detach(pkgns, unload = TRUE, character.only = TRUE, force = TRUE)
      success = TRUE
    }, error=function(e) warning(e))
    if (success)
      nUnloaded = nUnloaded + 1L
  }
  nUnloaded
}
