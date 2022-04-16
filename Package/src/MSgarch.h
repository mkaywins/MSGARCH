#ifndef MSGARCH_H  // include guard
#define MSGARCH_H

#include "SingleRegime.h"
#include <RcppArmadillo.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//======================================== AUXILIARY FUNCTIONS
//============================================//

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//========================================== MS-GARCH class
//===============================================//

// type definition (vector of pointers to Base)
typedef std::vector<Base*> many;
typedef std::vector<volatility> volatilityVector;

// MS-GARCH class
class MSgarch {
  many specs;           // vector of pointers to Base objects
  int K;                // number of models
  
  NumericMatrix P;      // transition-probability matrix
  NumericVector PLast;  // all transition-probability matrix at last step
  NumericVector P0;     // initial distribution of states
  double P_mean;        // mean for the prior on transition-probabilities
  double P_sd;          // sd for the prior on transition-probabilities
  double LND_MIN;       // minimum loglikelihood allowed
  
public:
  std::vector<std::string> name;
  NumericVector theta0;
  NumericVector Sigma0;
  CharacterVector label;
  NumericVector lower;
  NumericVector upper;
  NumericVector ineq_lb;
  NumericVector ineq_ub;
  IntegerVector NbParams;  // number of parameters for each model (excluding the
  // transition probabilities)
  IntegerVector NbParamsModel;
  
  // constructor
  MSgarch(List L) {
    //std::cout << "Called MSgarch constructorasfd..." << std::endl;
    // the input to the constuctor is a list of c++ objects like [[1]] C++ object <0x5653739eb4e0> of class 'sGARCH_norm' <0x56536c192e50>
    // extract pointers to models
    K = L.length();  // number of models
    Environment env;
    for (List::iterator it = L.begin(); it != L.end();
    ++it) {  // loop over models
      env = *it; // returns a pointer when de-referenced
      //std::cout << R_ExternalPtrAddr(env.get(".pointer")) << std::endl;
      specs.push_back(
        static_cast<Base*>(R_ExternalPtrAddr(env.get(".pointer"))));
      //std::cout << "spec " << env << std::endl;
      // casts pointer to sGARCH_norm to other pointer of Base
      // An object pointer can be explicitly converted to an object pointer of a diﬀerent type. When a prvalue v of type “pointer to T1” is converted to the type “pointer to cv T2”
    }
    
    
    
    // loop over models // many is a vector of pointers to Base objects
    for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
      name.push_back((*it)->spec_name());
      
      // creates a numeric vector for variables for each model in the MSGarch
      MyConcatenate(theta0, (*it)->spec_theta0()); // de-reference it to get the value - the value is a pointer to the base object in the "many" vector // the pointer gets (->) the member 'spec_theta0' // spec_theta0 gets the mean of prior distribution and is defined in Base and the mean is declared in sGARCH.h
      MyConcatenate(Sigma0, (*it)->spec_Sigma0()); //...
      MyConcatenate(label, (*it)->spec_label());
      MyConcatenate(lower, (*it)->spec_lower());
      MyConcatenate(upper, (*it)->spec_upper());
      MyConcatenate(ineq_lb, NumericVector::create((*it)->spec_ineq_lb()));
      MyConcatenate(ineq_ub, NumericVector::create((*it)->spec_ineq_ub()));
      NbParams.push_back((*it)->spec_nb_coeffs());
      NbParamsModel.push_back((*it)->spec_nb_coeffs_model());
    }
    P0 = rep(1.0 / K, K);       // initial probabilities are 1/K // returns a numericvector object
    PLast = rep(1.0 / K, K);
    P_mean = 1 / K;
    P_sd = 100;                 // sd on the prior is set to 100
    LND_MIN = log(DBL_MIN) + 1; // initialized the minimum allowed loglikihood - important for HamiltonFilter // DBL_MIN is the smallest double ! 
    
    // add transition-probabilities ( constraints)
    if (K > 1) {
      int NbP = K * (K - 1); // possible state to transition into ??
      NumericVector P_theta0  = rep(1.0 / K, NbP);  // theta0
      NumericVector P_Sigma0  = rep(1.0, NbP);      // Sigma0
      NumericVector P_lower   = rep(0.0, NbP);      // lower
      NumericVector P_upper   = rep(1.0, NbP);      // upper
      NumericVector P_ineq_lb = rep(0.0, K);       // ineq_lb
      NumericVector P_ineq_ub = rep(1.0, K);       // ineq_ub
      CharacterVector P_label(NbP, "P");           // label // creating a charactervector from "NbP" and "P"
      
      // concatenates all stuff above
      MyConcatenate(theta0, P_theta0); // this essentially concatenates theta0 from all the indiviual K models and the transition probabilities
      MyConcatenate(Sigma0, P_Sigma0);
      MyConcatenate(label, P_label);
      MyConcatenate(lower, P_lower);
      MyConcatenate(upper, P_upper);
      MyConcatenate(ineq_lb, P_ineq_lb);
      MyConcatenate(ineq_ub, P_ineq_ub);
      
    }
  }
  
  // set the parameters (including those of the distribution) of all models
  // the last elements of theta should be those of the transition-probability
  // matrix
  // this function should always be called first
  void loadparam(const NumericVector&);
  
  // to be called before 'calc_prior', 'ineq_func' or 'set_vol'
  void prep_ineq_vol() {
    for (many::iterator it = specs.begin(); it != specs.end(); ++it)
      (*it)->spec_prep_ineq_vol(); // for sGarch this can be empty 
  }
  
  // to be called before 'calc_kernel'
  void prep_kernel() {
    for (many::iterator it = specs.begin(); it != specs.end(); ++it)
      (*it)->spec_prep_kernel();
  }
  
  // loglikelihood of a single observation for all models
  NumericVector calc_kernel(const volatilityVector& vol, const double& yi) {
    NumericVector lnd(K);
    int k = 0;
    for (many::iterator it = specs.begin(); it != specs.end(); ++it) { // iterate over models in specs
      lnd[k] = (*it)->spec_calc_kernel(vol[k], yi); // the model specific method for calc_kernel
      k++;
    }
    return lnd;
  }
  
  // initialize all volatilities to their undonditional expected value
  volatilityVector set_vol() {
    volatilityVector vol(K);
    int k = 0;
    for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
      vol[k] = (*it)->spec_set_vol();
      k++;
    }
    return vol;
  }
  
  // increment all volatilities
  void increment_vol(volatilityVector& vol, const double& yim1) {
    int k = 0;
    for (many::iterator it = specs.begin(); it != specs.end(); ++it) { // iterate over all models
      (*it)->spec_increment_vol(vol[k], yim1);                         // increment the volatility - the volatility is passed by reference i.e. the vol struct itself is changed when we increment
      k++;
    }
  }
  
  NumericVector get_mean() {
    NumericVector out;
    for (many::iterator it = specs.begin(); it != specs.end();
    ++it) {  // loop over models
      MyConcatenate(out, (*it)->get_mean());
    }
    return (out);
  }
  
  void set_mean(NumericVector& new_mean) {
    int k = 0;
    NumericVector out(sum(NbParams));
    for (many::iterator it = specs.begin(); it != specs.end();
    ++it) {  // loop over models
      (*it)->set_mean(extract_theta_it(new_mean, k));
      k++;
    }
  }
  
  NumericVector get_sd() {
    NumericVector out;
    for (many::iterator it = specs.begin(); it != specs.end();
    ++it) {  // loop over models
      MyConcatenate(out, (*it)->get_sd());
    }
    return (out);
  }
  
  void set_sd(NumericVector& new_sd) {
    int k = 0;
    NumericVector out(sum(NbParams));
    for (many::iterator it = specs.begin(); it != specs.end();
    ++it) {  // loop over models
      (*it)->set_sd(extract_theta_it(new_sd, k));
      k++;
    }
  }
  
  // extract parameter vector of model 'k', where k is in [0, K-1]
  virtual NumericVector extract_theta_it(const NumericVector& theta, const int& k) {
    int start = MyCumsum(NbParams, k);
    NumericVector theta_it(theta.begin() + start,
                           theta.begin() + start + NbParams[k]);
    return theta_it;
  }
  
  // extract transition-probability from state 'k', where k is in [0, K-1]
  virtual NumericVector extract_P_it(const NumericVector& theta, const int& k) {
    int Tot_NbParams = sum(NbParams);
    NumericVector P_it(theta.begin() + Tot_NbParams + k * (K - 1),
                       theta.begin() + Tot_NbParams + (k + 1) * (K - 1)); // only computing K-1 transition probabilities 
    P_it.push_back(1 - sum(P_it)); // the last one is just 1 - rest
    return P_it;
  }
  
  // simulate a random inovation from model 'k', where k is in [0, K-1]
  double rndgen(const int& k) {
    many::iterator it = specs.begin() + k;
    return (*it)->spec_rndgen(1)[0];
  }
  
  // inequality function
  NumericVector ineq_func(const NumericVector& theta) {
    NumericVector out;
    loadparam(theta);
    prep_ineq_vol();
    for (many::iterator it = specs.begin(); it != specs.end(); ++it)
      out.push_back((*it)->spec_ineq_func());
    if (K > 1) {
      NumericMatrix Psub = P(Range(0, K - 1), Range(0, K - 2));
      for (int i = 0; i < K; i++) out.push_back(sum(Psub(i, _)));
    }
    return out;
  }
  
  // Getter functions for private members
  many get_specs() {return specs; }
  
  int get_K() {return K; }
  
  NumericVector get_PLast() { return PLast; }
  
  NumericMatrix get_P() {return P; }
  
  NumericVector get_P0() {return P0; }
  
  double get_P_mean() {return P_mean; }
  
  double get_P_sd() {return P_sd; }
  
  double get_LND_MIN() {return LND_MIN; }
  
  // Setter functions for private members
  
  void set_specs(many specs_){specs = specs_;}
  
  void set_K(int K_) {K = K_; }
  
  void set_PLast(NumericVector PLast_) {PLast = PLast_;}
  
  void set_P(NumericMatrix P_) {P = P_; }
  
  void set_P0(NumericVector P0_) {P0 = P0_; }
  
  void set_P_mean(double P_mean_) {P_mean = P_mean_; }
  
  void set_P_sd(double P_sd_) {P_sd = P_sd_; }
  
  void set_LIN_MIN(double LND_MIN_) {LND_MIN = LND_MIN_; }
  
  

  
  List f_get_Pstate(const NumericVector&, const NumericVector&);
  
  // check prior
  prior calc_prior(const NumericVector&);
  
  arma::cube calc_ht(NumericMatrix&, const NumericVector&);
  
  NumericVector f_pdf(const NumericVector&, const NumericVector&,
                      const NumericVector&, const bool&);
  
  arma::cube f_pdf_its(const NumericVector&, const NumericVector&,
                       const NumericMatrix&, const bool&);
  
  NumericVector f_cdf(const NumericVector&, const NumericVector&,
                      const NumericVector&, const bool&);
  
  arma::cube f_cdf_its(const NumericVector&, const NumericVector&,
                       const NumericMatrix&, const bool&);
  
  // model simulation
  Rcpp::List f_sim(const int&, const int&, const NumericVector&);
  
  Rcpp::List f_simAhead(const NumericVector&, const int&,  const int&, const NumericVector&,
                        const NumericVector&);
  
  Rcpp::List f_rnd(const int&, const NumericVector&, const NumericVector&);
  // compute loglikelihood matrix
  NumericMatrix calc_lndMat(const NumericVector&);
  
  NumericMatrix f_unc_vol(NumericMatrix&);
  
  // apply Hamilton filter
  double HamiltonFilter(const NumericMatrix&);
  
  // Model evaluation
  NumericVector eval_model(NumericMatrix&, const NumericVector&, const bool&);
};

//---------------------- load parameters of all models  ----------------------//
inline void MSgarch::loadparam(const NumericVector& theta) {
  //std::cout << "Call: loadparam" << std::endl;
  // load the parameters of each model and the transition-probability matrix
  NumericMatrix P_mat(K, K);
  int k = 0;
  for (many::iterator it = specs.begin(); it != specs.end();
  ++it) {  // loop over models

    NumericVector theta_it = extract_theta_it(theta, k);  // parameters of model 'it'
    NumericVector P_it = extract_P_it(theta, k);  // transition probabilities from model 'it'
    (*it)->spec_loadparam(theta_it);
    P_mat(k, _) = P_it; // each row corresponds to a state - the vector of probabilities is assigned to the rows
    k++;
  }
  //std::cout << "transitions probability matrix" << std::endl;
  //Rcout << P_mat << std::endl;
  P = P_mat; // transitions matrix
  
  // -- computing the unconditional probabilities pi -- 
  
  // pi = (1,...,1)(I - P + Q)^-1, where Q is a matrix of ones 
  arma::mat I = arma::eye(K,K); // identitiy matrix
  arma::mat Umat = arma::ones(K,K); // matrix of 1s
  arma::vec Uvec(K); // vector of 0s of size K
  Uvec.fill(1); // make it into a vector of 1s
  arma::mat foo = (I - as<arma::mat>(P_mat) + Umat).t();
  
  arma::vec delta = (foo).i() * Uvec; // .i() returns inverse of foo
  for(int i = 0; i < K; i++){
    P0(i) = delta(i); // set the unconditional probabilities as the initial probabilities P0 = (p, 1-p)
  }
}

//------------------------------ Prior calculation
//------------------------------//
inline prior MSgarch::calc_prior(const NumericVector& theta) {
  // compute prior of individual models
  bool r1_joint = 1;    // joint r1 of the models // we add the true/false to check wheter everything for all models is ok
  double r2_joint = 0;  // joint r2 of the models
  double r3_joint = 0;
  int k = 0;
  prior pr;             // prior is just a struct of r1, r2, r3
  for (many::iterator it = specs.begin(); it != specs.end();
  ++it) {  // loop over models
    NumericVector theta_it =
    extract_theta_it(theta, k);  // parameters of model 'it' as a NumericVector
    NumericVector P_it =
      extract_P_it(theta, k);   // transition probabilities from model 'it'
    pr = (*it)->spec_calc_prior(theta_it);  // prior of model 'it'
    
    r1_joint = r1_joint && pr.r1 && is_true(all((0 < P_it) & (P_it < 1))); // check wheter the parameters for the chosen distributions are right and wheter all elements in the transition probability matrix P_it are in (0,1)
    
    //std::cout << "calc_prior" << std::endl;
    //std::cout << "P_it: " << P_it << std::endl;
    //std::cout << "is_true(all((0 < P_it) & (P_it < 1))): " << is_true(all((0 < P_it) & (P_it < 1))) << std::endl;
    //std::cout << "pr.r1: " << pr.r1 << std::endl;
    
    r2_joint += pr.r2; 
    r3_joint += sum(dunif(P_it, 0.0, 1.0, 1)) + pr.r3;
    k++;
  }
  // return result
  prior out;
  out.r1 = r1_joint;
  out.r2 = ((out.r1) ? r2_joint : -1e10); // return r2_joint if the conditions for r1 are met oww. -1e10
  out.r3 = r3_joint;
  return out;
}

inline NumericMatrix MSgarch::f_unc_vol(NumericMatrix& all_thetas) {
  
  int nb_thetas = all_thetas.nrow();
  volatilityVector vol;
  NumericVector theta_j;
  NumericMatrix ht(nb_thetas, K); // create numeric matrix nb_thetas x K ht 
  
  for (int j = 0; j < nb_thetas; j++) {  // loop over vectors of parameters
    theta_j = all_thetas(j, _);
    loadparam(theta_j);
    prep_ineq_vol();
    vol = set_vol();
    // initialize volatility
    for (int s = 0; s < K; s++) {
      ht(j, s) = vol[s].h;
    }
  }
  return ht;
}

inline arma::cube MSgarch::calc_ht(NumericMatrix& all_thetas,
                                   const NumericVector& y) {
  
  int nb_obs = y.size();
  int nb_thetas = all_thetas.nrow();
  volatilityVector vol;
  NumericVector theta_j;
  arma::cube ht(nb_obs + 1, nb_thetas, K); // create a n+1 x 3 x K tensor/cube
  
  for (int j = 0; j < nb_thetas; j++) {  // loop over vectors of parameters
    theta_j = all_thetas(j, _);
    loadparam(theta_j);
    prep_ineq_vol();
    vol = set_vol();
    // initialize volatility for all states
    for (int s = 0; s < K; s++) {
      ht(0, j, s) = vol[s].h;
    }
    
    for (int i = 1; i <= nb_obs; i++) {  // loop over observations
      increment_vol(vol, y[i - 1]);      // increment all volatilities
      
      for (int s = 0; s < K; s++) {     // loop over all states to update ht for all states s
        ht(i, j, s) = vol[s].h;
      }
    }
  }
  return ht;
}

inline NumericVector MSgarch::f_pdf(const NumericVector& x,
                                    const NumericVector& theta,
                                    const NumericVector& y,
                                    const bool& is_log) {
  // computes volatility
  int s = 0;
  int nx = x.size();
  int ny = y.size();
  double sig;
  NumericVector tmp(nx);
  NumericVector out(nx);
  loadparam(theta);  // load parameters
  prep_ineq_vol();   // prepare functions related to volatility
  volatilityVector vol = set_vol();  // initialize volatility
  
  for (int t = 0; t < ny; t++) increment_vol(vol, y[t]);
  
  HamiltonFilter(calc_lndMat(y));
  
  for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
    sig = sqrt(vol[s].h);
    // computes PDF
    for (int i = 0; i < nx; i++) {
      tmp[i] = (*it)->spec_calc_pdf(x[i] / sig) / sig;  //
      out[i] = out[i] + tmp[i] * PLast[s];
    }
    s++;
  }
  
  if (is_log) {
    for (int i = 0; i < nx; i++) {
      out[i] = log(tmp[i]);
    }
  }
  
  return out;
}

inline arma::cube MSgarch::f_pdf_its(const NumericVector& theta,
                                     const NumericVector& y,
                                     const NumericMatrix& x,
                                     const bool& is_log) {
  // computes volatility
  int s = 0;
  int ny = y.size();
  int nx = x.nrow();
  double sig;
  arma::cube tmp(ny, nx, K);
  loadparam(theta);  // load parameters
  prep_ineq_vol();   // prepare functions related to volatility
  volatilityVector vol = set_vol();  // initialize volatility
  
  for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
    sig = sqrt(vol[s].h);
    for (int ix = 0; ix < nx; ix++) {
      tmp(0,ix , s) = (*it)->spec_calc_pdf(x(ix, 0) / sig) / sig;  //
    }
    s++;
  }
  
  for (int i = 1; i < ny; i++) {
    s = 0;
    increment_vol(vol, y[i - 1]);
    for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
      sig = sqrt(vol[s].h);
      for (int ix = 0; ix < nx; ix++) {
        tmp(i,ix, s) = (*it)->spec_calc_pdf(x(ix, i) / sig) / sig;  //
      }
      s++;
    }
  }
  
  return tmp;
}

inline NumericVector MSgarch::f_cdf(const NumericVector& x,
                                    const NumericVector& theta,
                                    const NumericVector& y,
                                    const bool& is_log) {
  // computes volatility
  int s = 0;
  int nx = x.size();
  int ny = y.size();
  double sig;
  NumericVector tmp(nx);
  NumericVector out(nx);
  loadparam(theta);  // load parameters
  prep_ineq_vol();   // prepare functions related to volatility
  volatilityVector vol = set_vol();  // initialize volatility
  
  for (int t = 0; t < ny; t++) increment_vol(vol, y[t]);
  
  HamiltonFilter(calc_lndMat(y));
  
  for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
    sig = sqrt(vol[s].h);
    // computes CDF
    for (int i = 0; i < nx; i++) {
      tmp[i] = (*it)->spec_calc_cdf(x[i] / sig);
      out[i] = out[i] + tmp[i] * PLast[s];
    }
    s++;
  }
  
  if (is_log) {
    for (int i = 0; i < nx; i++) {
      out[i] = log(tmp[i]);
    }
  }
  
  return out;
}

inline arma::cube MSgarch::f_cdf_its(const NumericVector& theta,
                                     const NumericVector& y,
                                     const NumericMatrix& x,
                                     const bool& is_log) {
  // computes volatility
  int s = 0;
  int ny = y.size();
  double sig;
  int nx = x.nrow();
  arma::cube tmp(ny,nx, K);
  loadparam(theta);  // load parameters
  prep_ineq_vol();   // prepare functions related to volatility
  volatilityVector vol = set_vol();  // initialize volatility
  for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
    sig = sqrt(vol[s].h);
    for (int ix = 0; ix < nx; ix++) {
      tmp(ix, 0, s) = (*it)->spec_calc_cdf(x(ix, 0) / sig);  //
    }
    s++;
  }
  for (int i = 1; i < ny; i++) {
    s = 0;
    increment_vol(vol, y[i - 1]);
    for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
      sig = sqrt(vol[s].h);
      for (int ix = 0; ix < nx; ix++) {
        tmp(i,ix, s) = (*it)->spec_calc_cdf(x(ix, i) / sig);  //
      }
      s++;
    }
  }
  
  return tmp;
}

//------------------------------ Model simulation
//------------------------------//
inline List MSgarch::f_sim(const int& n, const int& m, const NumericVector& theta) {
  // setup
  NumericMatrix y(m,n);  // observations
  NumericMatrix S(m,n);  // states of the Markov chain
  arma::cube CondVol(m,n,K);
  loadparam(theta);       // load parameters
  
  // generate first draw
  double z;            // random innovation from initial state
  prep_ineq_vol();                    // prep for 'set_vol'
  volatilityVector vol;  // initialize all volatilities
  for (int i = 0; i < m; i++) {
    S(i,0) = sampleState(P0);             // sample initial state
    z = rndgen(S(i,0));
    vol = set_vol();
    for (int s = 0; s < K; s++) {
      CondVol(i, 0, s) = sqrt(vol[s].h);
    }
    y(i,0) = z * sqrt(vol[S(i,0)].h);
    for (int t = 1; t < n; t++) {
      S(i,t) = sampleState(P(S(i,t - 1), _));  // sample new state
      z = rndgen(S(i,t));                    // sample new innovation
      increment_vol(vol, y(i,t - 1));        // increment all volatilities
      y(i,t) = z * sqrt(vol[S(i,t)].h);        // new draw
      for (int s = 0; s < K; s++) {
        CondVol(i, t, s) = sqrt(vol[s].h);
      }
    }
  }
  return (List::create(Rcpp::Named("draws") = y, Rcpp::Named("state") = S,  Rcpp::Named("CondVol") = CondVol));
}

inline List MSgarch::f_simAhead(const NumericVector& y, const int& n, const int& m,
                                const NumericVector& theta,
                                const NumericVector& P0_) {
  // setup
  int nb_obs = y.size();  // total number of observations to simulate
  NumericMatrix y_sim(m, n);
  NumericMatrix S(m, n);
  arma::cube CondVol(m,n,K);
  loadparam(theta);  // load parameters
  prep_ineq_vol();   // prep for 'set_vol'
  volatilityVector vol0 = set_vol();
  double z;
  for (int t = 1; t <= nb_obs; t++) {
    increment_vol(vol0, y[t - 1]);  // increment all volatilities
  }
  for (int i = 0; i < m; i++) {
    S(i,0) = sampleState(P0_);           // sample initial state
    z = rndgen(S(i,0));
    y_sim(i,0) = z * sqrt(vol0[S(i,0)].h);  // first draw
  }
  volatilityVector vol = vol0;
  for (int i = 0; i < m; i++) {
    for (int s = 0; s < K; s++) {
      CondVol(i, 0, s) = sqrt(vol[s].h);
    }
    for (int t = 1; t < n; t++) {
      S(i,t) = sampleState(P(S(i,t - 1), _));  // sample new state
      z = rndgen(S(i,t));                    // sample new innovation
      increment_vol(vol, y_sim(i,t - 1));    // increment all volatilities
      y_sim(i,t) = z * sqrt(vol[S(i,t)].h);
      for (int s = 0; s < K; s++) {
        CondVol(i, t, s) = sqrt(vol[s].h);
      }
    }  // new draw
    vol = vol0;
  }
  return (List::create(Rcpp::Named("draws") = y_sim, Rcpp::Named("state") = S,  Rcpp::Named("CondVol") = CondVol));
}

inline List MSgarch::f_rnd(const int& n, const NumericVector& theta,
                           const NumericVector& y) {
  // setup
  int nb_obs = y.size();
  NumericVector draw(n);  // draw
  IntegerVector S(n);     // states of the Markov chain
  loadparam(theta);       // load parameters
  double z;
  prep_ineq_vol();  // prep for 'set_vol'
  volatilityVector vol = set_vol();
  
  for (int i = 1; i <= nb_obs; i++) {  // loop over observations
    increment_vol(vol, y[i - 1]);      // increment all volatilities
  }
  HamiltonFilter(calc_lndMat(y));
  
  // increment over time
  for (int i = 0; i < n; i++) {
    S[i] = sampleState(PLast);        // sample new state
    z = rndgen(S[i]);                 // sample new innovation
    draw[i] = z * sqrt(vol[S[i]].h);  // new draw
  }
  NumericVector yy(draw.begin(), draw.end());
  NumericVector SS(S.begin(), S.end());
  return (
      Rcpp::List::create(Rcpp::Named("draws") = yy, Rcpp::Named("state") = SS));
}

//------------------------------ Compute loglikelihood matrix
//------------------------------//
inline NumericMatrix MSgarch::calc_lndMat(const NumericVector& y) {
  // set up
  int nb_obs = y.size();                        // number os observations
  NumericMatrix lndMat(K, nb_obs - 1);          // initialize a NumericMatrix with dim K x n-1
  
  // initialize
  volatilityVector vol = set_vol();            // initialize volatility to the unconditional expected value // returns a volatilityVector object because for each model there is a volatility to initialize
  prep_kernel();                               // needs to be run before calc_kernel - it is a distribution-member function // e.g. for normal it is empty
  
  // loop over observations
  for (int t = 1; t < nb_obs; t++) {
    increment_vol(vol, y[t - 1]);               // increment all volatilities // the vol object itself will be changed through increment_vol
    lndMat(_, t - 1) = calc_kernel(vol, y[t]);  // calc all kernels // takes in the current volatility and current observation y_t // returns a K size vector of loglikelihoods 
  }
  return lndMat;
}

//-------------------------------------  Hamilton filter
//-------------------------------------//
inline double MSgarch::HamiltonFilter(const NumericMatrix& lndMat) {
  // the input to the hamilton filter is a matrix of loglikelihoods for all observations (col) and states (row)
  //std::cout << "Call: HM Filter" << std::endl;
  int n_step = lndMat.ncol();               // ncol of lndMat = number of observations
  double lnd = 0, min_lnd, delta, sum_tmp;  // 
  NumericVector Pspot, Ppred, lndCol, tmp;  
  
  // first step
  Pspot = clone(P0);             // Prob(St-1 | I(t-1) // P0 are the initial probabilities (pi1, ..., piK)
  Ppred = matrixProd(Pspot, P);  // one-step-ahead Prob(St | I(t-1))
  lndCol = lndMat(_, 0);         // first column of lndMat
  min_lnd = min(lndCol),         // 
    delta =
      ((min_lnd < LND_MIN) ? LND_MIN - min_lnd : 0);  // handle over/under-flows // if min_lnd is too low, then the difference between the bound and the value is taken ow. 0
  
  tmp = Ppred *
    exp(lndCol + delta);  // unormalized one-step-ahead Prob(St | I(t)) // we exponentiate the (by underflow corrected) loglikelihood s.t. we get P(St| I(t)) that was computed before by calc_lndMat and outputted as log
  
  // remaining steps
  for (int t = 1; t < n_step; t++) { // loop over 1,...,n observations
    sum_tmp = sum(tmp);              // the sum of tmp gives the summed probability for all states i.e. "sum prob over all states"
    lnd += -delta + log(sum_tmp);    // !!! increment loglikelihood // we take the log of the summed probabilities (and correct for underflow)
    Pspot = tmp / sum_tmp;           // Prob(St-1 | I(t-1)) / sum 
    Ppred = matrixProd(Pspot, P);    // Prob(St | I(t-1)) one step ahead probability
    lndCol = lndMat(_, t);           // taking the t-th column of the loglik matrix 
    min_lnd = min(lndCol),
      delta = ((min_lnd < LND_MIN) ? LND_MIN - min_lnd
                 : 0);                 // handle over/under-flows // like above ...
    tmp = Ppred * exp(lndCol + delta); // unormalized one-step-ahead Prob(St | I(t))
  }
  sum_tmp = sum(tmp);
  lnd += -delta + log(sum_tmp);  // increment loglikelihood
  Pspot = tmp / sum_tmp;
  PLast = matrixProd(Pspot, P);    // the final transition matrix
  return lnd;
}

inline List MSgarch::f_get_Pstate(const NumericVector& theta,
                                  const NumericVector& y) {
  // init
  loadparam(theta);  // load parameters
  prep_ineq_vol();   // prepare functions related to volatility
  volatilityVector vol = set_vol();   // initialize volatility
  NumericMatrix lndMat = calc_lndMat(y);  // likelihood in each state
  
  int n_step = lndMat.ncol();
  double lnd = 0, min_lnd, delta, sum_tmp;
  NumericVector Pspot, Ppred, lndCol, tmp;
  arma::mat PtmpSpot(n_step + 1, K);
  arma::mat PtmpPred(n_step + 2, K);
  arma::mat PtmpSmooth(n_step + 2, K);
  
  // first step
  Pspot = clone(P0);             // Prob(St | I(t))
  Ppred = matrixProd(Pspot, P);  // one-step-ahead Prob(St | I(t-1))
  for (int i = 0; i < K; i++) {
    PtmpSpot(0, i) = Pspot(i);
  }
  
  for (int i = 0; i < K; i++) {
    PtmpPred(0, i) = Pspot(i);
  }
  for (int i = 0; i < K; i++) {
    PtmpPred(1, i) = Ppred(i);
  }
  lndCol = lndMat(_, 0);
  min_lnd = min(lndCol),
    delta =
      ((min_lnd < LND_MIN) ? LND_MIN - min_lnd : 0);  // handle over/under-flows
  tmp = Ppred *
    exp(lndCol + delta);  // unormalized one-step-ahead Prob(St | I(t))
  
  // remaining steps
  for (int t = 1; t < n_step; t++) {
    sum_tmp = sum(tmp);
    lnd += -delta + log(sum_tmp);  // increment loglikelihood
    Pspot = tmp / sum_tmp;
    Ppred = matrixProd(Pspot, P);
    for (int i = 0; i < K; i++) {
      PtmpSpot(t, i) = Pspot(i);
    }
    for (int i = 0; i < K; i++) {
      PtmpPred(t + 1, i) = Ppred(i);
    }
    lndCol = lndMat(_, t);
    min_lnd = min(lndCol),
      delta = ((min_lnd < LND_MIN) ? LND_MIN - min_lnd
                 : 0);  // handle over/under-flows
    tmp = Ppred * exp(lndCol + delta);
  }
  sum_tmp = sum(tmp);
  lnd += -delta + log(sum_tmp);  // increment loglikelihood
  Pspot = tmp / sum_tmp;
  PLast = matrixProd(Pspot, P);
  
  for (int i = 0; i < K; i++) {
    PtmpSpot(n_step, i) = Pspot(i);
  }
  for (int i = 0; i < K; i++) {
    PtmpPred(n_step + 1, i) = PLast(i);
  }
  for (int i = 0; i < K; i++) {
    PtmpSmooth(n_step + 1, i) = PLast(i);
  }
  arma::mat tmpMat(1, 2);
  arma::mat tmpMat2(2, 1);
  for (int t = n_step; t >= 0; t--) {
    tmpMat = (PtmpSmooth.row(t + 1) / PtmpPred.row(t + 1));
    tmpMat2 = (as<arma::mat>(P) * tmpMat.t());
    PtmpSmooth.row(t) = PtmpSpot.row(t) % tmpMat2.t();
  }
  
  return List::create(
    Rcpp::Named("FiltProb") = PtmpSpot, Rcpp::Named("PredProb") = PtmpPred,
    Rcpp::Named("SmoothProb") = PtmpSmooth, Rcpp::Named("LL") = lndMat);
}

//------------------------------------- Model evaluation
//-------------------------------------//
inline NumericVector MSgarch::eval_model(NumericMatrix& all_thetas,
                                         const NumericVector& y,
                                         const bool& do_prior) {
  //std::cout << "Call: c MSgarch::eval_model..." << std::endl;
  
  //Rcout << all_thetas << std::endl;
  //Rcout << all_thetas << std::endl; 
  // set up
  int nb_thetas = all_thetas.nrow(); // all_thetas = par ie. a0, a1, b, a0, a1, b, p1, p2
  NumericVector lnd(nb_thetas), theta_j(all_thetas.ncol()); // create 2 NumericVector lnd, theta_j
  prior pr;
  double tmp;
  
  // loop over each vector of parameters
  for (int j = 0; j < nb_thetas; j++) { // loop over all models
    theta_j = all_thetas(j, _);  // extract parameters // all_thetas(j,_) = all_thetas[j,] as numeric vector 1xncol
    
    //std::cout << "thetaj: " << std::endl;
    //Rcout << theta_j << std::endl;
    
    loadparam(theta_j);          // load parameters for the object associated with MSgarch::
    prep_ineq_vol();             // needs to be called before 'calc_prior' // for sGarch it does nothing
    pr = calc_prior(theta_j);    // will calculate the model-specific prior - returns a prior struct
    
    if (do_prior == true) {      // usually this is set to false
      lnd[j] = pr.r2 + pr.r3;
    } else {
      lnd[j] = pr.r2;           // add the computed value for the model-specific prior to the vector lnd
    }
    pr = calc_prior(theta_j);   // maybe a duplicate from above ??
    tmp = 0;
    
    // after the prior is computed we need to compute the rest of the loglikelihood w. the rest of the observations
    if (pr.r1) tmp += HamiltonFilter(calc_lndMat(y)); // if the conditions are met for prior.r1, then we add the loglikelihood from the HamiltonFilter to the prior
    lnd[j] += tmp;
  }
  return lnd;
}



#endif  // MSgarch.h
