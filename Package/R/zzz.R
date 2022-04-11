# It's a file where one usually puts actions on load of the package. 
# It is tradition/convention that it's called zzz.R and could be called anything.R

# 'loadModule' - imports all modules that were defined in src to expose them in the package:MSGARCH environment
# all object that are defined within the respective modules are then available to the package
loadModule("TVMSgarch", TRUE)
loadModule("MSgarch", TRUE)
loadModule("sARCH", TRUE)
loadModule("sGARCH", TRUE)
loadModule("eGARCH", TRUE)
loadModule("gjrGARCH", TRUE)
loadModule("tGARCH", TRUE)
loadModule("Normal", TRUE)
loadModule("Student", TRUE)
loadModule("Ged", TRUE)

#' @importFrom utils getFromNamespace
#' @importFrom utils tail
#' @importFrom stats as.ts 
#' @importFrom stats is.ts 
#' @importFrom stats pnorm
#' @importFrom stats pnorm 
#' @importFrom stats AIC 
#' @importFrom stats BIC 
#' @importFrom methods is
.onAttach <- function(libname, pkgname){
  
  MSGARCH_env <- as.environment("package:MSGARCH")
  # a bunch of distribution classes are retrieved from the MSGARCH namespace and assinged to .<<distribution>>_created 
  # they are assingedto envir = "package:MSGARCH"
  assign(".norm_sym_created", new(Class = getFromNamespace("norm_sym", ns = "MSGARCH")), envir = MSGARCH_env)
  assign(".norm_skew_created", new(Class = getFromNamespace("norm_skew", ns = "MSGARCH")), envir = MSGARCH_env)
  assign(".std_sym_created", new(Class = getFromNamespace("std_sym", ns = "MSGARCH")), envir = MSGARCH_env)
  assign(".std_skew_created", new(Class = getFromNamespace("std_skew", ns = "MSGARCH")), envir = MSGARCH_env)
  assign(".ged_sym_created", new(Class = getFromNamespace("ged_sym", ns = "MSGARCH")), envir = MSGARCH_env)
  assign(".ged_skew_created", new(Class = getFromNamespace("ged_skew", ns = "MSGARCH")), envir = MSGARCH_env)
  




}

