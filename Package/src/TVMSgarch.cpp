#include "TVMSgarch.h"


//===========================================================================//
//===========================================================================//
//================================= TVMSgarch =================================//
//===========================================================================//
//===========================================================================//
RCPP_MODULE(TVMSgarch){
  class_<MSgarch>("MSgarch")
  .constructor<List>()
  .field("name", &MSgarch::name)
  .field("theta0", &MSgarch::theta0)
  .field("Sigma0", &MSgarch::Sigma0)
  .field("label", &MSgarch::label)
  .field("lower", &MSgarch::lower)
  .field("upper", &MSgarch::upper)
  .field("ineq_lb", &MSgarch::ineq_lb)
  .field("ineq_ub", &MSgarch::ineq_ub)
  .field("NbParams", &MSgarch::NbParams)
  .field("NbParamsModel", &MSgarch::NbParamsModel)
  .method("f_sim", &MSgarch::f_sim)
  .method("f_simAhead", &MSgarch::f_simAhead)
  .method("f_get_sd", &MSgarch::get_sd)
  .method("f_set_sd", &MSgarch::set_sd)
  .method("f_get_mean", &MSgarch::get_mean)
  .method("f_set_mean", &MSgarch::set_mean)
  .method("calc_ht", &MSgarch::calc_ht)
  .method("ineq_func", &MSgarch::ineq_func)
  .method("f_pdf", &MSgarch::f_pdf)
  .method("f_pdf_its", &MSgarch::f_pdf_its)
  .method("f_cdf", &MSgarch::f_cdf)
  .method("f_cdf_its", &MSgarch::f_cdf_its)
  .method("f_rnd", &MSgarch::f_rnd)
  .method("f_unc_vol", &MSgarch::f_unc_vol);
  ;
  class_<TVMSgarch>("TVMSgarch")
    .derives<MSgarch>("MSgarch")
    .constructor<List>()
    .field("NbFactors", &TVMSgarch::NbFactors)
    .method("f_get_Pstate", &TVMSgarch::f_get_Pstate)
    .method("eval_model", &TVMSgarch::eval_model)
    .method("get_all_Pt", &TVMSgarch::get_all_Pt)
  ;
}