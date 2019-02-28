// A conditioned Wright-Fisher diffusion and its application to inferring natural selection and allele age from time series data of allele frequencies
// Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

// this version is able to handle missing values in DNA data

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
