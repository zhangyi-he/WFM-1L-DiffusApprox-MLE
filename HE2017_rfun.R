#' @title A conditioned Wright-Fisher diffusion and its application to inferring natural selection and allele age from time series data of allele frequencies
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' this version is unable to handle missing values in DNA data

#' R functions

#install.packages("inline")
library("inline")
#install.packages("Rcpp")
library("Rcpp")
#install.packages("RcppArmadillo")
library("RcppArmadillo")

#install.packages("compiler")
library("compiler")
#enableJIT(1)

# call C++ functions
sourceCpp("./Code/Code v1.1/HE2017_cfun.cpp")

################################################################################

#' Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_ale_frq the initial mutant allele frequency of the population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory

#' Standard version
simulateOLWFMS <- function(sel_cof, dom_par, pop_siz, int_ale_frq, int_gen, lst_gen) {
  ale_frq_pth <- simulateOLWFMS_arma(sel_cof, dom_par, pop_siz, int_ale_frq, int_gen, lst_gen)
  ale_frq_pth <- as.vector(ale_frq_pth)
  
  return(ale_frq_pth)
}
#' Compiled version
cmpsimulateOLWFMS <- cmpfun(simulateOLWFMS)

########################################

#' Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_ale_frq the initial mutant allele frequency of the population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method
#' @param data_augmentation = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

#' Standard version
simulateOLWFDS <- function(sel_cof, dom_par, pop_siz, int_ale_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE) {
  ale_frq_pth <- simulateOLWFDS_arma(sel_cof, dom_par, pop_siz, int_ale_frq, int_gen, lst_gen, ptn_num)
  ale_frq_pth <- as.vector(ale_frq_pth)
  
  # return the simulated sample trajectory without data augmentation
  if (data_augmentation == FALSE) {
    ale_frq_pth <- ale_frq_pth[(0:(lst_gen - int_gen)) * ptn_num + 1]
  }
  
  return(ale_frq_pth)
}
#' Compiled version
cmpsimulateOLWFDS <- cmpfun(simulateOLWFDS)

########################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_gen the initial generation
#' @param int_ale_frq the initial mutant allele frequency of the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_chr_cnt the count of the chromosomes drawn from the population at all sampling time points
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Standard version
simulateHMM <- function(model, sel_cof, dom_par, pop_siz, int_gen, int_ale_frq, smp_gen, smp_chr_cnt, ...) {
  if (is.na(int_gen)) {
    fst_gen <- min(smp_gen)
    lst_gen <- max(smp_gen)
  } else {
    fst_gen <- min(smp_gen, int_gen)
    lst_gen <- max(smp_gen, int_gen)
  }
  
  # generate the population allele frequency trajectory
  if (model == "WFM") {
    repeat {
      pop_ale_frq <- cmpsimulateOLWFMS(sel_cof, dom_par, pop_siz, int_ale_frq, fst_gen, lst_gen)
      if (tail(pop_ale_frq, 1) > 0 & tail(pop_ale_frq, 1) < 1) {
        break
      }
    }
  }
  if (model == "WFD") {
    repeat {
      pop_ale_frq <- cmpsimulateOLWFDS(sel_cof, dom_par, pop_siz, int_ale_frq, fst_gen, lst_gen, ptn_num, data_augmentation = FALSE)
      if (tail(pop_ale_frq, 1) > 0 & tail(pop_ale_frq, 1) < 1) {
        break
      }
    }
  }
  
  # generate the sample allele counts at all sampling time points
  smp_ale_cnt <- numeric(length(smp_gen))
  for (k in 1:length(smp_gen)) {
    smp_ale_cnt[k] <- rbinom(1, size = smp_chr_cnt[k], prob = pop_ale_frq[smp_gen[k] - fst_gen + 1])
  }
  
  return(list(int_gen = int_gen, 
              smp_gen = smp_gen, 
              smp_chr_cnt = smp_chr_cnt, 
              smp_ale_cnt = smp_ale_cnt, 
              pop_ale_frq = pop_ale_frq))
}
#' Compiled version
cmpsimulateHMM <- cmpfun(simulateHMM)

########################################



#' Solve the Kolmogorov backward equation (KBE) associated with the one-locus Wright-Fisher diffusion with selection using the Crank–Nicolson method (CRM)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_gen the initial generation
#' @param int_ale_frq the initial mutant allele frequency of the population
#' @param ptn_num_gen the number of the subintervals divided per generation in the Crank–Nicolson method
#' @param ptn_num_frq the number of the subintervals divided for the allele frequency in the Crank–Nicolson method
#' @param thr_gen_num the threshold of the maximum of the allele age
#' @param thr_tot_mas the threshold of the remaining probability mass 

#' Standard version
runCRM <- function(sel_cof, dom_par, pop_siz, int_gen, int_ale_frq, ptn_num_gen, ptn_num_frq, thr_gen_num, thr_tot_mas) {
  frq_grd_dtn <- head(rev(1 - log(1:(ptn_num_frq + 1)) / log(ptn_num_frq + 1)), round(ptn_num_frq / 3) + 1)
  frq_grd_fxn <- tail(log(1:(ptn_num_frq + 1)) / log(ptn_num_frq + 1), round(ptn_num_frq / 3) + 1)
  frq_grd <- c(frq_grd_dtn, seq(tail(frq_grd_dtn, 1), head(frq_grd_fxn, 1), length.out = ptn_num_frq - 2 * round(ptn_num_frq / 3) + 1), frq_grd_fxn)
  frq_grd <- frq_grd[!duplicated(frq_grd)]
  
  # run the CRM
  CRM <- runCRM_arma(sel_cof, dom_par, pop_siz, int_ale_frq, ptn_num_gen, ptn_num_frq, frq_grd, thr_gen_num, thr_tot_mas)
  sol_ale_frq <- CRM$sol_frq
  sol_ale_age <- as.vector(CRM$sol_age)
  
  gen_grd <- int_gen - ((ncol(sol_ale_frq) - 1):0) / ptn_num_gen
  
  return(list(gen_grd = gen_grd, 
              frq_grd = frq_grd, 
              sol_ale_frq = sol_ale_frq, 
              sol_ale_age = sol_ale_age))
}
#' Compiled version
cmprunCRM <- cmpfun(runCRM)

########################################

#' Run the Kolmogorov backward equation (KBE) associated with the one-locus Wright-Fisher diffusion with selection using the Crank–Nicolson method (CRM)
#' Parameter settings
#' @param sim_sel_cof the selection coefficient generated from the particle marginal Metropolis-Hastings
#' @param sim_dom_par the dominance parameter generated from the particle marginal Metropolis-Hastings
#' @param sim_pop_siz the number of the diploid individuals in the population generated from the particle marginal Metropolis-Hastings
#' @param sim_obs_pop_ale_frq the mutant allele frequency of the population at the sampling time point when the mutant allele was first observed in the sample generated from the particle marginal Metropolis-Hastings
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the number of the diploid individuals drawn from the population at all sampling time points
#' @param smp_ale_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num_gen the number of the subintervals divided per generation in the Crank–Nicolson method
#' @param ptn_num_frq the number of the subintervals divided for the allele frequency in the Crank–Nicolson method
#' @param thr_gen_num the threshold of the maximum of the allele age
#' @param thr_tot_mas the threshold of the remaining probability mass
#' @param sim_num the number of the potential values of the allele age generated

#' Standard version
runKBE <- function(sim_sel_cof, sim_dom_par, sim_pop_siz, sim_obs_pop_ale_frq, smp_gen, smp_siz, smp_ale_cnt, ptn_num_gen, ptn_num_frq, thr_gen_num, thr_tot_mas, sim_num) {
  smp_gen <- floor(smp_gen)
  
  smp_gen <- smp_gen[1:min(which(smp_ale_cnt > 0))]
  smp_siz <- smp_siz[1:min(which(smp_ale_cnt > 0))]
  smp_ale_cnt <- smp_ale_cnt[1:min(which(smp_ale_cnt > 0))]
  
  frq_grd_dtn <- head(rev(1 - log(1:(ptn_num_frq + 1)) / log(ptn_num_frq + 1)), round(ptn_num_frq / 3) + 1)
  frq_grd_fxn <- tail(log(1:(ptn_num_frq + 1)) / log(ptn_num_frq + 1), round(ptn_num_frq / 3) + 1)
  frq_grd <- c(frq_grd_dtn, seq(tail(frq_grd_dtn, 1), head(frq_grd_fxn, 1), length.out = ptn_num_frq - 2 * round(ptn_num_frq / 3) + 1), frq_grd_fxn)
  frq_grd <- frq_grd[!duplicated(frq_grd)]
  
  KBE <- runKBE_arma(sim_sel_cof, sim_dom_par, sim_pop_siz, sim_obs_pop_ale_frq, smp_gen, smp_siz, smp_ale_cnt, ptn_num_gen, ptn_num_frq, frq_grd, thr_gen_num, thr_tot_mas, sim_num)
  sim_ale_age <- KBE$ale_age
  
  return(list(sim_ale_age = sim_ale_age))
}
#' Compiled version
cmprunKBE <- cmpfun(runKBE)

########################################

#' Run the backward procedure (for the inference of natural selection and allele age only)
#' Parameter settings
#' @param sim_sel_cof the selection coefficient generated from the particle marginal Metropolis-Hastings
#' @param sim_dom_par the dominance parameter generated from the particle marginal Metropolis-Hastings
#' @param sim_pop_siz the number of the diploid individuals in the population generated from the particle marginal Metropolis-Hastings
#' @param sim_obs_pop_ale_frq the mutant allele frequency of the population at the sampling time point when the mutant allele was first observed in the sample generated from the particle marginal Metropolis-Hastings
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the number of the diploid individuals drawn from the population at all sampling time points
#' @param smp_ale_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num_gen the number of the subintervals divided per generation in the Crank–Nicolson method
#' @param ptn_num_frq the number of the subintervals divided for the allele frequency in the Crank–Nicolson method
#' @param thr_gen_num the threshold of the maximum of the allele age
#' @param thr_tot_mas the threshold of the remaining probability mass
#' @param grd_num the number of the grids in the kernel density estimation
#' @param sim_num the number of the potential values of the population genetic parameters generated

#' Standard version
runBkdProcedure <- function(sim_sel_cof, sim_dom_par, sim_pop_siz, sim_obs_pop_ale_frq, smp_gen, smp_siz, smp_ale_cnt, ptn_num_gen, ptn_num_frq, thr_gen_num, thr_tot_mas, grd_num, sim_num) {
  smp_gen <- floor(smp_gen)
  
  smp_gen <- smp_gen[1:min(which(smp_ale_cnt > 0))]
  smp_siz <- smp_siz[1:min(which(smp_ale_cnt > 0))]
  smp_ale_cnt <- smp_ale_cnt[1:min(which(smp_ale_cnt > 0))]
  
  frq_grd_dtn <- head(rev(1 - log(1:(ptn_num_frq + 1)) / log(ptn_num_frq + 1)), round(ptn_num_frq / 3) + 1)
  frq_grd_fxn <- tail(log(1:(ptn_num_frq + 1)) / log(ptn_num_frq + 1), round(ptn_num_frq / 3) + 1)
  frq_grd <- c(frq_grd_dtn, seq(tail(frq_grd_dtn, 1), head(frq_grd_fxn, 1), length.out = ptn_num_frq - 2 * round(ptn_num_frq / 3) + 1), frq_grd_fxn)
  frq_grd <- frq_grd[!duplicated(frq_grd)]
  
  # run the KBE
  KBE <- runKBE_arma(sim_sel_cof, sim_dom_par, sim_pop_siz, sim_obs_pop_ale_frq, smp_gen, smp_siz, smp_ale_cnt, ptn_num_gen, ptn_num_frq, frq_grd, thr_gen_num, thr_tot_mas, sim_num)
  sim_ale_age <- KBE$ale_age
  sim_ale_age <- as.vector(t(apply(sim_ale_age, 2, 
                                   function(ale_age) {
                                     pdf <- density(ale_age)
                                     return(pdf$x[which(pdf$y == max(pdf$y))])
                                   })))
  
  #sim_ale_age <- as.vector(t(sim_ale_age))
  
  #sim_sel_cof <- rep(sim_sel_cof, sim_num)
  #sim_dom_par <- rep(sim_dom_par, sim_num)
  #sim_pop_siz <- rep(sim_pop_siz, sim_num)
  #sim_obs_pop_ale_frq <- rep(sim_obs_pop_ale_frq, sim_num)
  
  # estimate the joint posterior probability density function for the selection coefficient and the allele age
  if (length(sim_sel_cof) < 1e+05) {
    jnt_pdf <- kde2d(sim_sel_cof, sim_ale_age, n = grd_num)
    sel_cof_grd <- jnt_pdf$x
    ale_age_grd <- jnt_pdf$y
  } else {
    jnt_pdf <- kde2d(tail(sim_sel_cof, 1e+05), tail(sim_ale_age, 1e+05), n = grd_num)
    sel_cof_grd <- jnt_pdf$x
    ale_age_grd <- jnt_pdf$y
  }
  
  return(list(jnt_pdf = jnt_pdf, 
              sim_sel_cof = sim_sel_cof, 
              sim_dom_par = sim_dom_par, 
              sim_pop_siz = sim_pop_siz, 
              sim_ale_age = sim_ale_age, 
              sim_obs_pop_ale_frq = sim_obs_pop_ale_frq))
}
#' Compiled version
cmprunBkdProcedure <- cmpfun(runBkdProcedure)

################################################################################
