
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
    if (grepl("[\\\\/]", pkgnm) && file.exists(pkgnm))
    {
      pkglocalname = sub("^([^_]+?)_?(?:\\d+?\\.)+?tar\\.gz", '\\1', basename(pkgnm), perl=T)
    } else if (pkgnm %in% names(pkgpaths))
    {
      instdir = pkgpaths[[pkgnm]]
    }
    .give.status(percent=as.integer(100 * i / tot),
                 message=sprintf("Instalando pacote: %s", pkglocalname),
                 engMsg=sprintf("Installing package: %s", pkglocalname))
    utils::install.packages(pkgnm, quiet = T, verbose = F, lib = instdir, type='both')
    if (!(pkglocalname %in% getAllInstalledPackages()))
    {
      .give.status(percent=as.integer(100 * i / tot),
                   message=sprintf("Instalando pacote (2ª tentativa): %s", pkglocalname),
                   engMsg=sprintf("Installing package (2nd attempt): %s", pkglocalname))
      utils::install.packages(pkglocalname, quiet = T, verbose = F, lib = instdir, type='both')
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
  pkgs = names(utils::installed.packages(fields = 'Package')[,'Package'])
  pkgs
}

# [[geapexport string[] GetLoadedPackages()]]
# Gets a character vector listing all currently loaded packages
#' @export
getLoadedPackages <- function()
{
  pkgs = unique(sub("^package:", "", grep("^package:", search(), value = T, perl = T), perl = T))
  pkgs
}
