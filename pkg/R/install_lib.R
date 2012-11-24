libs=c("seqinr","phangorn","optimx","expm","multicore","minqa","mgcv","numDeriv","ppso")
type=getOption("pkgType")                           
CheckInstallPackage <- function(packages, repos="http://cran.r-project.org",
                                depend=c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances"), ...) {
  installed=as.data.frame(installed.packages())
  for(p in packages) {
    if(is.na(charmatch(p, installed[,1]))) { 
      install.packages(p, repos=repos, dependencies=depend, ...) 
    }
  }
} 
CheckInstallPackage(packages=libs)