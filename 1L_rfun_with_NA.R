#' @title A conditioned Wright-Fisher diffusion and its application to inferring natural selection and allele age from time series data of allele frequencies
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' this version is able to handle missing values in DNA data

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
sourceCpp("1L_cfun_with_NA.cpp")

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
#' @param missing = TRUE/FALSE (return the observations with missing values or not)
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_gen the initial generation
#' @param int_ale_frq the initial mutant allele frequency of the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_chr_cnt the count of the chromosomes drawn from the population at all sampling time points
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Standard version
simulateHMM <- function(model, missing, sel_cof, dom_par, pop_siz, int_gen, int_ale_frq, smp_gen, smp_chr_cnt, ...) {
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
  if (missing == TRUE) {
    smp_ale_cnt <- numeric(length(smp_gen))
    mis_ale_cnt <- numeric(length(smp_gen))
    for (k in 1:length(smp_gen)) {
      smp_ale_cnt[k] <- rbinom(1, size = smp_chr_cnt[k], prob = pop_ale_frq[smp_gen[k] - fst_gen + 1])
      mis_ale_cnt[k] <- sample(0:round(smp_chr_cnt[k] / 4 / k), size = 1)
      smp_ale_cnt[k] <- smp_ale_cnt[k] - sample(0:min(smp_ale_cnt[k], mis_ale_cnt[k]), size = 1)
    }
    
    return(list(int_gen = int_gen, 
                smp_gen = smp_gen, 
                smp_chr_cnt = smp_chr_cnt, 
                smp_ale_cnt = smp_ale_cnt, 
                mis_ale_cnt = mis_ale_cnt, 
                pop_ale_frq = pop_ale_frq))
  } else {
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
}
#' Compiled version
cmpsimulateHMM <- cmpfun(simulateHMM)

########################################



################################################################################
