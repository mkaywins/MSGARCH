f_nll <- function(vPw, data, spec, do.plm) {
  #print("Call: f_nll")
  # if there were no names for the parameters vPw - assign spec$label to the names
  if (is.null(names(vPw))) {
    vPw <- f_rename_par(vPw, spec)
  }

  # applies transform methods to ensure that parameters during the (optimisation) are always in bounds
  vPn <- f_mapPar(vPw, spec, do.plm) # ERROR: HERE the factors and probs are turned

  # if parameters in spec were specified that should stay fixed, they are added here
  if (isTRUE(spec$fixed.pars.bool)) {
    vPn <- f_add_fixedpar(vPn, spec$fixed.pars)
    vPn <- vPn[spec$label]
  }

  # if regime specific parameters in spec were specified that should stay fixed, they are added here
  if (isTRUE(spec$regime.const.pars.bool)) {
    vPn <- f_add_regimeconstpar(vPn, spec$K, spec$label)
  }
  
  # computes the log likelihood of the model
  dLLK <- Kernel(spec, vPn, data, log = TRUE, do.prior = FALSE)

  # if the log likelihood is not finite - set it to -1e+10
  if (!is.finite(dLLK)) {
    dLLK <- -1e+10
  }
  
  # return negative log likelihood
  return(-dLLK)
}

#' @importFrom stats logLik
#' @export
logLik.MSGARCH_ML_FIT <- function(object, ...){
  # this function provides the 'logLik' method to the output object of FitML() such that AIC and BIC from stats can be computed 
  out = structure(object$loglik, df = dofMSGARCH(object), 
                  nobs = length(object$data))
  return(out)
}