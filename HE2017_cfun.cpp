// A conditioned Wright-Fisher diffusion and its application to inferring natural selection and allele age from time series data of allele frequencies
// Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

// this version is unable to handle missing values in DNA data

// C functions

#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]
#include <math.h>

using namespace Rcpp;
using namespace std;
// using namespace arma;

/********** WFM **********/
// Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::drowvec simulateOLWFMS_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const double& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare the mutant allele frequency trajectory
  arma::drowvec ale_frq_pth(arma::uword(lst_gen - int_gen) + 1);
  
  // initialise the mutant allele frequency in generation 0
  ale_frq_pth(0) = int_frq;
  
  // declare the fitness
  arma::dcolvec fts = arma::zeros<arma::dcolvec>(3);
  fts(0) = 1.0;
  fts(1) = 1.0 - sel_cof * dom_par;
  fts(2) = 1.0 - sel_cof;
  
  // declare and initialise the genotype frequencies during a single generation of the life cycle
  arma::dcolvec gen_frq = arma::zeros<arma::dcolvec>(3);
  
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    // random union of gametes
    gen_frq(0) = ale_frq_pth(t - 1) * ale_frq_pth(t - 1);
    gen_frq(1) = 2 * ale_frq_pth(t - 1) * (1 - ale_frq_pth(t - 1));
    gen_frq(2) = (1 - ale_frq_pth(t - 1)) * (1 - ale_frq_pth(t - 1));
    
    // viability selection
    gen_frq = (fts % gen_frq) / arma::accu(fts % gen_frq);
    
    // meiosis (calculate the sampling probability)
    double prob = gen_frq(0) + gen_frq(1) / 2;
    
    // reproduction (the Wright-Fisher sampling)
    ale_frq_pth(t) = R::rbinom(2 * pop_siz, prob) / 2 / pop_siz;
  }
  
  return ale_frq_pth;
}
/*************************/


/********** WFD **********/
// Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
// [[Rcpp::export]]
arma::drowvec simulateOLWFDS_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const double& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // rescale the selection coefficient
  double scl_sel_cof = 2 * pop_siz * sel_cof;
  
  // declare delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  // declare delta W
  arma::drowvec dW = pow(dt, 0.5) * arma::randn<arma::drowvec>(arma::uword(lst_gen - int_gen) * ptn_num);
  
  // declare the mutant allele frequency trajectory
  arma::drowvec ale_frq_pth(arma::uword(lst_gen - int_gen) * ptn_num + 1);
  
  // initialise the mutant allele frequency in generation 0
  ale_frq_pth(0) = int_frq;
  
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient
    double mu = scl_sel_cof * ale_frq_pth(t - 1) * (1 - ale_frq_pth(t - 1)) * ((1 - dom_par) - (1 - 2 * dom_par) * ale_frq_pth(t - 1));
    
    // calculate the diffusion coefficient
    double sigma = pow(ale_frq_pth(t - 1) * (1 - ale_frq_pth(t - 1)), 0.5);
    
    // proceed the Euler-Maruyama scheme
    ale_frq_pth(t) = ale_frq_pth(t - 1) + mu * dt + sigma * dW(t - 1);
    
    // remove the noise from the numerical techniques
    if (ale_frq_pth(t) < 0) {
      ale_frq_pth(t) = 0;
    }
    if (ale_frq_pth(t) > 1) {
      ale_frq_pth(t) = 1;
    }
  }
  
  return ale_frq_pth;
}
/*************************/



/********** KBE **********/
// Solve the Kolmogorov backward equation associated with the Wright-Fisher diffusion using the Crank–Nicolson method
// [[Rcpp::export]]
List solveKBE_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const double& int_frq, const arma::uword& ptn_num_gen, const arma::uword& ptn_num_frq, const arma::dcolvec& frq_grd, const arma::uword& thr_num_gen, const double& thr_tot_mas) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // rescale the selection coefficient
  double scl_sel_cof = 2 * pop_siz * sel_cof;
  
  arma::dcolvec mu = scl_sel_cof * frq_grd % (1 - frq_grd) % ((1 - dom_par) - (1 - 2 * dom_par) * frq_grd);
  arma::dcolvec Sigma = frq_grd % (1 - frq_grd);
  
  double delta_gen = 1.0 / (2 * pop_siz) / ptn_num_gen;
  arma::dcolvec delta_frq = diff(frq_grd, 1);
  
  arma::dcolvec frq_bin = frq_grd.tail(ptn_num_frq - 1) - frq_grd.head(ptn_num_frq - 1);
  frq_bin = frq_bin / 2;
  frq_bin.insert_rows(0, 1);
  frq_bin.head(1) = delta_frq.head(1) / 2;
  frq_bin.insert_rows(frq_bin.n_elem, 1);
  frq_bin.tail(1) = delta_frq.tail(1) / 2;
  
  arma::dmat tridiag_fwd = arma::zeros<arma::dmat>(ptn_num_frq + 1, ptn_num_frq + 1);
  for (arma::uword i = 1; i < ptn_num_frq; i++) {
    tridiag_fwd(i, i - 1) = delta_gen * ((Sigma(i) / 2) - mu(i) * delta_frq(i - 1)) / (delta_frq(i - 1) * (delta_frq(i - 1) + delta_frq(i)));
    tridiag_fwd(i, i) = 1 - delta_gen * (Sigma(i) / 2) / (delta_frq(i - 1) * delta_frq(i));
    tridiag_fwd(i, i + 1) = delta_gen * ((Sigma(i) / 2) + mu(i) * delta_frq(i)) / (delta_frq(i) * (delta_frq(i - 1) + delta_frq(i)));
  }
  
  arma::dmat tridiag_bkd = arma::zeros<arma::dmat>(ptn_num_frq + 1, ptn_num_frq + 1);
  for (arma::uword i = 1; i < ptn_num_frq; i++) {
    tridiag_bkd(i, i - 1) = -delta_gen * ((Sigma(i) / 2) - mu(i) * delta_frq(i - 1)) / (delta_frq(i - 1) * (delta_frq(i - 1) + delta_frq(i)));
    tridiag_bkd(i, i) = 1 + delta_gen * (Sigma(i) / 2) / (delta_frq(i - 1) * delta_frq(i));
    tridiag_bkd(i, i + 1) = -delta_gen * ((Sigma(i) / 2) + mu(i) * delta_frq(i)) / (delta_frq(i) * (delta_frq(i - 1) + delta_frq(i)));
  }
  tridiag_bkd(0, 0) = 1;
  tridiag_bkd(0, 1) = -1;
  tridiag_bkd(ptn_num_frq, ptn_num_frq - 1) = -1;
  tridiag_bkd(ptn_num_frq, ptn_num_frq) = 1;
  
  arma::dmat tridiag = tridiag_bkd.i() * tridiag_fwd;
  
  arma::drowvec sol_age = arma::zeros<arma::drowvec>(1);
  arma::dmat sol_frq = arma::zeros<arma::dmat>(ptn_num_frq + 1, 1);
  
  sol_frq(arma::find(frq_grd > int_frq).min()) = 1.0 / frq_bin(arma::find(frq_grd > int_frq).min());
  
  double tot_mas = arma::sum(sol_frq % frq_bin);
  
  for(arma::uword k = 1; k < thr_num_gen * ptn_num_gen + 1; k++) {
    arma::dcolvec sol_tmp  = sol_frq.col(0);
    
    tot_mas = arma::sum(sol_tmp % frq_bin);
    if (tot_mas > thr_tot_mas) {
      sol_tmp = tridiag * sol_tmp;
      sol_tmp.elem(arma::find(sol_tmp < 0)).zeros();
      tot_mas = arma::sum(sol_tmp % frq_bin);
    } else {
      break;
    }
    
    sol_age.insert_cols(0, 1);
    sol_age.head(1) = sol_tmp.head(1) * frq_bin.head(1);
    
    sol_tmp.head(1) = 0;
    double tot_mas_nonzero = arma::sum(sol_tmp % frq_bin);
    
    if (tot_mas < tot_mas_nonzero) {
      sol_age.shed_col(0);
      break;
    } else {
      sol_tmp.tail(1) = 0;
      sol_tmp = sol_tmp * tot_mas_nonzero / arma::sum(sol_tmp % frq_bin);
      
      sol_frq.insert_cols(0, 1);
      sol_frq.col(0) = sol_tmp;
    }
  }
  
  return List::create(Named("sol_frq", sol_frq), 
                      Named("sol_age", sol_age));
}

// Compute the transition probability density function of the Wright-Fisher diffusion
// [[Rcpp::export]]



// Solve the Kolmogorov backward equation associated with the conditioned Wright-Fisher diffusion using the Crank–Nicolson method
// [[Rcpp::export]]
List solveCondKBE_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const double& int_frq, const arma::uword& ptn_num_gen, const arma::uword& ptn_num_frq, const arma::dcolvec& frq_grd, const arma::uword& thr_num_gen, const double& thr_tot_mas) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // rescale the selection coefficient
  double scl_sel_cof = 2 * pop_siz * sel_cof;
  
  arma::dcolvec mu = scl_sel_cof * frq_grd % (1 - frq_grd) % ((1 - dom_par) - (1 - 2 * dom_par) * frq_grd);
  arma::dcolvec Sigma = frq_grd % (1 - frq_grd);
  
  double delta_gen = 1.0 / (2 * pop_siz) / ptn_num_gen;
  arma::dcolvec delta_frq = diff(frq_grd, 1);
  
  arma::dcolvec frq_bin = frq_grd.tail(ptn_num_frq - 1) - frq_grd.head(ptn_num_frq - 1);
  frq_bin = frq_bin / 2;
  frq_bin.insert_rows(0, 1);
  frq_bin.head(1) = delta_frq.head(1) / 2;
  frq_bin.insert_rows(frq_bin.n_elem, 1);
  frq_bin.tail(1) = delta_frq.tail(1) / 2;
  
  arma::dmat tridiag_fwd = arma::zeros<arma::dmat>(ptn_num_frq + 1, ptn_num_frq + 1);
  for (arma::uword i = 1; i < ptn_num_frq; i++) {
    tridiag_fwd(i, i - 1) = delta_gen * ((Sigma(i) / 2) - mu(i) * delta_frq(i - 1)) / (delta_frq(i - 1) * (delta_frq(i - 1) + delta_frq(i)));
    tridiag_fwd(i, i) = 1 - delta_gen * (Sigma(i) / 2) / (delta_frq(i - 1) * delta_frq(i));
    tridiag_fwd(i, i + 1) = delta_gen * ((Sigma(i) / 2) + mu(i) * delta_frq(i)) / (delta_frq(i) * (delta_frq(i - 1) + delta_frq(i)));
  }
  
  arma::dmat tridiag_bkd = arma::zeros<arma::dmat>(ptn_num_frq + 1, ptn_num_frq + 1);
  for (arma::uword i = 1; i < ptn_num_frq; i++) {
    tridiag_bkd(i, i - 1) = -delta_gen * ((Sigma(i) / 2) - mu(i) * delta_frq(i - 1)) / (delta_frq(i - 1) * (delta_frq(i - 1) + delta_frq(i)));
    tridiag_bkd(i, i) = 1 + delta_gen * (Sigma(i) / 2) / (delta_frq(i - 1) * delta_frq(i));
    tridiag_bkd(i, i + 1) = -delta_gen * ((Sigma(i) / 2) + mu(i) * delta_frq(i)) / (delta_frq(i) * (delta_frq(i - 1) + delta_frq(i)));
  }
  tridiag_bkd(0, 0) = 1;
  tridiag_bkd(0, 1) = -1;
  tridiag_bkd(ptn_num_frq, ptn_num_frq - 1) = -1;
  tridiag_bkd(ptn_num_frq, ptn_num_frq) = 1;
  
  arma::dmat tridiag = tridiag_bkd.i() * tridiag_fwd;
  
  arma::drowvec sol_age = arma::zeros<arma::drowvec>(1);
  arma::dmat sol_frq = arma::zeros<arma::dmat>(ptn_num_frq + 1, 1);
  
  sol_frq(arma::find(frq_grd > int_frq).min()) = 1.0 / frq_bin(arma::find(frq_grd > int_frq).min());
  
  double tot_mas = arma::sum(sol_frq % frq_bin);
  
  for(arma::uword k = 1; k < thr_num_gen * ptn_num_gen + 1; k++) {
    arma::dcolvec sol_tmp  = sol_frq.col(0);
    
    tot_mas = arma::sum(sol_tmp % frq_bin);
    if (tot_mas > thr_tot_mas) {
      sol_tmp = tridiag * sol_tmp;
      sol_tmp.elem(arma::find(sol_tmp < 0)).zeros();
      tot_mas = arma::sum(sol_tmp % frq_bin);
    } else {
      break;
    }
    
    sol_age.insert_cols(0, 1);
    sol_age.head(1) = sol_tmp.head(1) * frq_bin.head(1);
    
    sol_tmp.head(1) = 0;
    double tot_mas_nonzero = arma::sum(sol_tmp % frq_bin);
    
    if (tot_mas < tot_mas_nonzero) {
      sol_age.shed_col(0);
      break;
    } else {
      sol_tmp.tail(1) = 0;
      sol_tmp = sol_tmp * tot_mas_nonzero / arma::sum(sol_tmp % frq_bin);
      
      sol_frq.insert_cols(0, 1);
      sol_frq.col(0) = sol_tmp;
    }
  }
  
  return List::create(Named("sol_frq", sol_frq), 
                      Named("sol_age", sol_age));
}

// Compute the transition probability density function of the conditioned Wright-Fisher diffusion
// [[Rcpp::export]]



/*************************/


/******* Posterior *******/
// Compute the posterior probability distribution for the selection coefficient and the allele age under the Wright-Fisher diffusion
// [[Rcpp::export]]



// Compute MAP estimates for the selection coefficient and the allele age under the Wright-Fisher diffusion
// [[Rcpp::export]]



// Compute the posterior probability distribution for the selection coefficient and the allele age under the conditioned Wright-Fisher diffusion
// [[Rcpp::export]]



// Compute MAP estimates for the selection coefficient and the allele age under the conditioned Wright-Fisher diffusion
// [[Rcpp::export]]



/*************************/
