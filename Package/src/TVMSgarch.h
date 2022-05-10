#ifndef TVMSGARCH_H  // include guard
#define TVMSGARCH_H

#include "MSgarch.h"
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
class TVMSgarch : public MSgarch {
  
public:
  int NbFactors; // number of parameters for the transition probability function 
  
  // constructor
  TVMSgarch(List L) : MSgarch(L) {}
  
  // time-varying transition probability
  NumericMatrix Pt(NumericVector all_factors, NumericVector Z){
    // t: time
    // all_factors: vector of factors 
    // Z: covariate vector of size n+1 for time t
    
    
    int nZ = Z.length(); // length of one covariate vector at time t
    int nfactor = all_factors.length();
    int K = get_K(); // get K from MSgarch class
    NumericVector num, denom, prob;
    double tmp;
    
    NumericMatrix P_mat(K,K); // cerate transition matrix
    std::fill( P_mat.begin(), P_mat.end(), NumericVector::get_na() ); // Na vals
    
    // fill the matrix with the correct probabilities
    //std::cout << "nfactor: " << nfactor << std::endl;
    //std::cout << "nZ: " << nZ << std::endl;
    //std::cout << "K: " << K << std::endl;
    
    for(int i = 0; i < nfactor; i += nZ ){
      // take nfactors each time and consturct the proability by means of multinomial logit
      num.push_back( exp(MyInner(all_factors[Range(i, i + nZ - 1)],  Z) ) );
    }
    
    
    int j = 0;
    //std::cout << "numeratoor:" << std::endl;
    //std::cout << num << std::endl;
    
    for(int i = 0; i < num.length(); i += (K-1), j++){
      // row sums
      NumericVector tmp = num[Range(i, i + K - 2)] / (1 + sum(num[Range(i, i + K - 2)]));
      tmp.push_back(1 / (1 + sum(num[Range(i, i + K - 2)])));
      P_mat(j, _) = tmp; 
    }
    
    
    //std::cout << "-----------P_mat from Pt() --------------" << std::endl;
    //std::cout << P_mat << std::endl;
    
    return(P_mat);
  }
  
  List get_all_Pt(NumericVector& all_thetas, const NumericMatrix& Z){
    
    int n = Z.nrow();
    NbFactors = Z.ncol();
    List all_Pt(n);
    
    NumericVector all_factors = extract_factors(all_thetas);
    
    for (int t = 0; t < n; t++) { 
      NumericMatrix P_t = Pt(all_factors, Z(t,_));
      all_Pt(t) = P_t;
    }
    return (all_Pt);
  }
  
  // set the parameters (including those of the distribution) of all models
  // the last elements of theta should be those of the transition-probability
  // matrix
  // this function should always be called first
  void loadparam(const NumericVector& theta, NumericMatrix Z);
  
  // extract parameter vector of model 'k', where k is in [0, K-1]
  NumericVector extract_theta_it(const NumericVector& theta, const int& k) override {
    //std::cout << "Call: extract_theta_it" << std::endl;
    
    //std::cout << "k: " << k << std::endl;
    int start = MyCumsum(NbParams, k);
    //std::cout << "start: " << start << std::endl;
    
    //std::cout << "theta.begin() + start: " << theta.begin() + start << std::endl;
    
    NumericVector theta_it(theta.begin() + start,
                           theta.begin() + start + NbParams[k]);
    //std::cout << "theta_it: " << theta_it << std::endl;
    return theta_it;
  }
  
  NumericVector extract_factors(const NumericVector& theta){
    //std::cout << "Call: extract_factors_it" << std::endl;
    //Rcout << theta << std::endl;
    
    int K = get_K();
    int Tot_NbParams = sum(NbParams);
    
    
    NumericVector factors_it(theta.begin() + Tot_NbParams,
                             theta.begin() + Tot_NbParams + (K) * (K-1) * NbFactors);
    
    //std::cout << "factors_it: " << std::endl;
    //Rcout << factors_it << std::endl;
    
    return factors_it;
  }
  
  // extract transition-probability from state 'k', where k is in [0, K-1]
  NumericVector extract_P_it(const NumericVector& theta, const int& k) override {
    //std::cout << "Call: extract_P_it" << std::endl;
    //Rcout << theta << std::endl;
    
    int K = get_K();
    int Tot_NbParams = sum(NbParams) + NbFactors * (K-1) * (K);
    //std::cout << "sum(NbParams): " << sum(NbParams) << std::endl;
    //std::cout << "NbFactors: " << NbFactors << std::endl;
    //std::cout << "Tot_NbParams + k * (K - 1): " << Tot_NbParams + k * (K - 1) << std::endl;
    NumericVector P_it(theta.begin() + Tot_NbParams + k * (K - 1),
                       theta.begin() + Tot_NbParams + (k + 1) * (K - 1)); // only computing K-1 transition probabilities 
    P_it.push_back(1 - sum(P_it)); // the last one is just 1 - rest
    
    //std::cout << "----------------- P_it: -------------" << std::endl;
    //std::cout << P_it << std::endl;
    
    return P_it;
  }

  // compute loglikelihood matrix
  NumericMatrix calc_lndMat(const NumericVector&);
  
  // apply Hamilton filter
  double HamiltonFilter(const NumericMatrix& lndMat,
                        const NumericVector& all_factors,
                        const NumericMatrix& Z);
  
  // get state probabilities
  List f_get_Pstate(const NumericVector&, 
                    const NumericVector&,
                    const NumericMatrix& Z);
  
  // model simulation
  Rcpp::List f_sim(const int&, const int&, const NumericVector&, const NumericMatrix& Z);
  
  // model simulation
  Rcpp::List f_simAhead(const NumericVector&, const int&,  const int&, const NumericVector&,
                        const NumericVector&, const NumericMatrix& Z);
  
  // Model evaluation
  NumericVector eval_model(NumericMatrix& all_thetas,
                           const NumericVector& y,
                           const NumericMatrix& Z,
                           const bool& do_prior);
};

//---------------------- load parameters of all models  ----------------------//
inline void TVMSgarch::loadparam(const NumericVector& theta, NumericMatrix Z) {
  //std::cout << "Call: ----------TVMSgarch loadparam --------------" << std::endl;
  // theta example: a0, a1, b, a0, a1, b g0 g1 g0 g1 g0 g1 g0 g1 p1 p2
  
  
  int K = get_K();             // get K from MSgarch class
  many specs = get_specs();    // get specs from MSgarch class
  NumericVector P0 = get_P0(); // get P0 from MSgarch class
  
  // get the number of factors
  NbFactors = Z.ncol(); 
  
  // get the number of model parameters (in total)
  int Total_NbParams = sum(NbParams);
  
  // initializing the transition probability matrix
  NumericMatrix P_mat(K, K);
  int k = 0;
  for (many::iterator it = specs.begin(); it != specs.end(); ++it) {  // loop over models
    //std::cout << "loop loadparam" << std::endl;
    NumericVector theta_it = extract_theta_it(theta, k);  // parameters of model 'it'
    NumericVector P_it = extract_P_it(theta, k);  // transition probabilities from model 'it'
    (*it)->spec_loadparam(theta_it); // set the model params in each model
    P_mat(k, _) = P_it; // each row corresponds to a state - the vector of probabilities is assigned to the rows
    k++;
  }
  
  // (initial) transitions matrix
  //P = P_mat;
  set_P(P_mat);
  
  // -- computing the unconditional probabilities pi -- 
  
  // pi = (1,...,1)(I - P + Q)^-1, where Q is a matrix of ones 
  arma::mat I = arma::eye(K,K);     // identitiy matrix
  arma::mat Umat = arma::ones(K,K); // matrix of 1s
  arma::vec Uvec(K);                // vector of 0s of size K
  Uvec.fill(1);                     // make it into a vector of 1s
  arma::mat foo = (I - as<arma::mat>(P_mat) + Umat).t();
  
  arma::vec delta = (foo).i() * Uvec; // .i() returns inverse of foo
  for(int i = 0; i < K; i++){
    P0(i) = delta(i); // set the unconditional probabilities as the initial probabilities P0 = (p, 1-p)
  }
  
  set_P0(P0);
}





//------------------------------ Model simulation
//------------------------------//
inline List TVMSgarch::f_sim(const int& n, const int& m, const NumericVector& theta,
                             const NumericMatrix& Z) {
  // setup
  NumericMatrix y(m,n);  // observations
  NumericMatrix S(m,n);  // states of the Markov chain
  NumericMatrix P_t;
  NumericVector P0 = get_P0();       // get P0      from MSgarch class
  NumericMatrix P = get_P();         // get P       from MSgarch class
  int K = get_K();
  NumericVector all_factors = extract_factors(theta);   // to get all the factors from the vector all_theta
  arma::cube CondVol(m,n,K);
  loadparam(theta, Z);       // load parameters
  
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
      P_t = Pt(all_factors, Z(t,_));     // transition prob at time t
      S(i,t) = sampleState(P_t(S(i,t - 1), _));  // sample new state
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

inline List TVMSgarch::f_simAhead(const NumericVector& y, const int& n, const int& m,
                                  const NumericVector& theta, const NumericVector& P0_,
                                  const NumericMatrix& Z) {
  // setup
  int nb_obs = y.size();  // total number of observations to simulate
  NumericMatrix y_sim(m, n);
  NumericMatrix S(m, n);
  NumericMatrix P_t;

  int K = get_K();
  NumericVector all_factors = extract_factors(theta);   // to get all the factors from the vector all_theta
  arma::cube CondVol(m,n,K);
  loadparam(theta, Z);               // load parameters
  prep_ineq_vol();                   // prep for 'set_vol'
  NumericMatrix P = get_P();         // get P       from MSgarch class
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


//------------------------------ Compute loglikelihood matrix
//------------------------------//
inline NumericMatrix TVMSgarch::calc_lndMat(const NumericVector& y) {
  // set up
  int nb_obs = y.size();                        // number os observations
  int K = get_K();                              // get K from MSgarch class
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
inline double TVMSgarch::HamiltonFilter(const NumericMatrix& lndMat,
                                        const NumericVector& all_theta,
                                        const NumericMatrix& Z) {

  //std::cout << "Call: HM Filter" << std::endl;
  int n_step = lndMat.ncol();               // ncol of lndMat = number of observations
  
  double lnd = 0, min_lnd, delta, sum_tmp;  // 
  NumericVector Pspot, Ppred, lndCol, tmp;  
  NumericMatrix P_t;
  NumericVector P0 = get_P0();       // get P0      from MSgarch class
  NumericMatrix P = get_P();         // get P       from MSgarch class
  double LND_MIN = get_LND_MIN();    // get LND_MIN from MSgarch class
  NumericVector PLast = get_PLast(); // get PLast   from MSgarch class
  NumericVector all_factors = extract_factors(all_theta);   // to get all the factors from the vector all_theta
  
  // first step
  Pspot = clone(P0);
  
  
  // first step
  Pspot = clone(P0);             // Prob(St-1 | I(t-1) // P0 are the initial probabilities (pi1, ..., piK)
  Ppred = matrixProd(Pspot, P);  // one-step-ahead Prob(St | I(t-1))
  lndCol = lndMat(_, 0);         // first column of lndMat
  min_lnd = min(lndCol),         // 
    delta =
      ((min_lnd < LND_MIN) ? LND_MIN - min_lnd : 0);  // handle over/under-flows // if min_lnd is too low, then the difference between the bound and the value is taken ow. 0
  
  tmp = Ppred *            // this gives the numerator of (11.2) in the book
    exp(lndCol + delta);  // unormalized one-step-ahead Prob(St | I(t)) // we exponentiate the (by underflow corrected) loglikelihood s.t. we get P(St| I(t)) that was computed before by calc_lndMat and outputted as log
  
  // remaining steps
  for (int t = 1; t < n_step; t++) { // loop over 1,...,n observations
    sum_tmp = sum(tmp);              // the sum of tmp gives the summed probability for all states i.e. "sum prob over all states"
    lnd += -delta + log(sum_tmp);    // !!! increment loglikelihood // we take the log of the summed probabilities (and correct for underflow)
    Pspot = tmp / sum_tmp;           // Prob(St-1 | I(t-1)) / sum 
    P_t = Pt(all_factors, Z(t,_));     // transition prob at time t
    Ppred = matrixProd(Pspot, P_t);   // Prob(St | I(t-1)) one step ahead probability
    lndCol = lndMat(_, t);           // taking the t-th column of the loglik matrix 
    min_lnd = min(lndCol),
      delta = ((min_lnd < LND_MIN) ? LND_MIN - min_lnd
                 : 0);                 // handle over/under-flows // like above ...
    tmp = Ppred * exp(lndCol + delta); // unormalized one-step-ahead Prob(St | I(t))
  }
  sum_tmp = sum(tmp);
  lnd += -delta + log(sum_tmp);  // increment loglikelihood
  Pspot = tmp / sum_tmp;
  PLast = matrixProd(Pspot, P_t);    // the final transition matrix
  
  set_PLast(PLast);
  
  return lnd;
}


inline List TVMSgarch::f_get_Pstate(const NumericVector& theta,
                                  const NumericVector& y,
                                  const NumericMatrix& Z) {
  // init
  loadparam(theta, Z);  // load parameters
  prep_ineq_vol();   // prepare functions related to volatility
  volatilityVector vol = set_vol();   // initialize volatility
  NumericMatrix lndMat = calc_lndMat(y);  // likelihood in each state
  
  
  NumericVector P0 = get_P0();       // get P0      from MSgarch class
  NumericMatrix P = get_P();         // get P       from MSgarch class
  NumericMatrix P_t = P;
  double LND_MIN = get_LND_MIN();    // get LND_MIN from MSgarch class
  NumericVector PLast = get_PLast(); // get PLast   from MSgarch class
  NumericVector all_factors = extract_factors(theta);   // to get all the factors from the vector all_theta
  int K = get_K();
  
  
    
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
  for (int t = 1; t < n_step; t++) { // loop over 1,...,n observations
    sum_tmp = sum(tmp);              // the sum of tmp gives the summed probability for all states i.e. "sum prob over all states"
    lnd += -delta + log(sum_tmp);    // !!! increment loglikelihood // we take the log of the summed probabilities (and correct for underflow)
    Pspot = tmp / sum_tmp;           // Prob(St-1 | I(t-1)) / sum 
    P_t = Pt(all_factors, Z(t,_));   // transition prob at time t
    Ppred = matrixProd(Pspot, P_t);  // Prob(St | I(t-1)) one step ahead probability
    
    for (int i = 0; i < K; i++) {
      PtmpSpot(t, i) = Pspot(i);
    }
    
    for (int i = 0; i < K; i++) {
      PtmpPred(t + 1, i) = Ppred(i);
    }
    
    lndCol = lndMat(_, t);           // taking the t-th column of the loglik matrix 
    min_lnd = min(lndCol),
      delta = ((min_lnd < LND_MIN) ? LND_MIN - min_lnd
                 : 0);                 // handle over/under-flows // like above ...
    tmp = Ppred * exp(lndCol + delta); // unormalized one-step-ahead Prob(St | I(t))
  }
  
  sum_tmp = sum(tmp);
  lnd += -delta + log(sum_tmp);  // increment loglikelihood
  Pspot = tmp / sum_tmp;
  PLast = matrixProd(Pspot, P_t); 
  
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
    P_t = Pt(all_factors, Z(t,_));     // transition prob at time t
    tmpMat2 = (as<arma::mat>(P_t) * tmpMat.t());
    PtmpSmooth.row(t) = PtmpSpot.row(t) % tmpMat2.t();
  }
  
  return List::create(
    Rcpp::Named("FiltProb") = PtmpSpot, Rcpp::Named("PredProb") = PtmpPred,
    Rcpp::Named("SmoothProb") = PtmpSmooth, Rcpp::Named("LL") = lndMat);
}

//------------------------------------- Model evaluation
//-------------------------------------//
inline NumericVector TVMSgarch::eval_model(NumericMatrix& all_thetas,
                                         const NumericVector& y,
                                         const NumericMatrix& Z,
                                         const bool& do_prior) {

  // all_theta = a0, a1, b, a0, a1, b, g00, g01, g10, g11, p1, p2
  //Rcout << all_thetas << std::endl;
  //Rcout << y << std::endl;
  //Rcout << Z << std::endl;
  //Rcout << do_prior << std::endl;
  // set up
  int nb_thetas = all_thetas.nrow(); // all_thetas = par ie. a0, a1, b, a0, a1, b, g00, g01, g10, g11
  NumericVector lnd(nb_thetas), theta_j(all_thetas.ncol()); // create 2 NumericVector lnd, theta_j
  prior pr;
  double tmp;
  
  // loop over each vector of parameters
  for (int j = 0; j < nb_thetas; j++) { // loop over all models
    theta_j = all_thetas(j, _);  // extract parameters // all_thetas(j,_) = all_thetas[j,] as numeric vector 1xncol
    loadparam(theta_j, Z);       // load parameters for the object associated with MSgarch::
    prep_ineq_vol();             // needs to be called before 'calc_prior' // for sGarch it does nothing
    pr = calc_prior(theta_j);    // will calculate the model-specific prior - returns a prior struct
    if (do_prior == true) {      // usually this is set to false
      lnd[j] = pr.r2 + pr.r3;
    } else {
      lnd[j] = pr.r2;           // add the computed value for the model-specific prior to the vector lnd
    }
    pr = calc_prior(theta_j);   // maybe a duplicate from above ??
    
    //std::cout << "r1: " << pr.r1 << std::endl;
    //std::cout << "r2: " << pr.r2 << std::endl;
    //std::cout << "r3: " << pr.r3 << std::endl;
    
    tmp = 0;
    // after the prior is computed we need to compute the rest of the loglikelihood w. the rest of the observations
    if (pr.r1) tmp += HamiltonFilter(calc_lndMat(y), theta_j, Z); // if the conditions are met for prior.r1, then we add the loglikelihood from the HamiltonFilter to the prior
    lnd[j] += tmp;
  }
  return lnd;
}



#endif  // TVMSgarch.h
