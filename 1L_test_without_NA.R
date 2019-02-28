#' @title A conditioned Wright-Fisher diffusion and its application to inferring natural selection and allele age from time series data of allele frequencies
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' this version is unable to handle missing values in DNA data

source("1L_rfun_without_NA.R")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("ggplot2")  
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

################################################################################

#' Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_ale_frq the initial mutant allele frequency of the population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory

sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_ale_frq <- 2e-01
int_gen <- 0
lst_gen <- 500

ale_frq_pth <- cmpsimulateOLWFMS(sel_cof, dom_par, pop_siz, int_ale_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, ale_frq_pth, type = 'l', lwd = 1.5, 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A mutant allele frequency generated with the Wright-Fisher model")

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

sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_ale_frq <- 2e-01
int_gen <- 0
lst_gen <- 500
ptn_num <- 1e+01

ale_frq_pth <- cmpsimulateOLWFDS(sel_cof, dom_par, pop_siz, int_ale_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE)

t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / pop_siz
plot(t, ale_frq_pth, type = 'l', lwd = 1.5, 
     xlab = "Time", ylab = "Allele frequency", 
     main = "A mutant allele frequency generated with the Wright-Fisher diffusion")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_ale_frq <- 2e-01
int_gen <- 0
lst_gen <- 500
ptn_num <- 1e+01
sim_num <- 1e+06

sim_ale_frq_WFM <- numeric(sim_num)
sim_ale_frq_WFD <- numeric(sim_num)
for (i in 1:sim_num) {
  print(i)
  sim_ale_frq_WFM[i] <- cmpsimulateOLWFMS(sel_cof, dom_par, pop_siz, int_ale_frq, int_gen, lst_gen)[(lst_gen - int_gen) + 1]
  sim_ale_frq_WFD[i] <- cmpsimulateOLWFDS(sel_cof, dom_par, pop_siz, int_ale_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)[(lst_gen - int_gen) + 1]
}

save(sel_cof, sel_cof, pop_siz, int_ale_frq, int_gen, lst_gen, ptn_num, sim_num, sim_ale_frq_WFM, sim_ale_frq_WFD, 
     file = "TEST_1L_comparison_WFM_and_WFD.rda")

load("TEST_1L_comparison_WFM_and_WFD.rda")

pdf(file = "TEST_1L_comparison_WFM_and_WFD.pdf", width = 20, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(sim_ale_frq_WFM, breaks = seq(min(sim_ale_frq_WFM, sim_ale_frq_WFD), max(sim_ale_frq_WFM, sim_ale_frq_WFD), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_ale_frq_WFM, sim_ale_frq_WFD), max(sim_ale_frq_WFM, sim_ale_frq_WFD)), 
     xlab = "Allele frequency", main = paste("Histograms of the allele frequency in generation", lst_gen, "under the W-F model and the W-F diffusion"))
hist(sim_ale_frq_WFD, breaks = seq(min(sim_ale_frq_WFM, sim_ale_frq_WFD), max(sim_ale_frq_WFM, sim_ale_frq_WFD), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
dev.off()

################################################################################

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

#' Simulate the dataset under the Wright-Fisher model 
model <- "WFM"
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_gen <- -200
int_ale_frq <- 1 / (2 * pop_siz)
smp_gen <- (0:5) * 100
smp_chr_cnt <- rep(100, 6)

sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_gen, int_ale_frq, smp_gen, smp_chr_cnt)
int_gen <- sim_HMM_WFM$int_gen
smp_gen <- sim_HMM_WFM$smp_gen
smp_chr_cnt <- sim_HMM_WFM$smp_chr_cnt
smp_ale_cnt <- sim_HMM_WFM$smp_ale_cnt
smp_ale_frq <- sim_HMM_WFM$smp_ale_cnt / sim_HMM_WFM$smp_chr_cnt
pop_ale_frq <- sim_HMM_WFM$pop_ale_frq

k <- min(int_gen, smp_gen):max(int_gen, smp_gen)
plot(k, pop_ale_frq, type = 'l', lwd = 1.5, 
     xlim = c(min(int_gen, smp_gen), max(int_gen, smp_gen)), ylim = c(min(smp_ale_frq, pop_ale_frq), max(smp_ale_frq, pop_ale_frq)), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the mutant allele generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq, col = 'red', pch = 17, cex = 1)

####################

#' Simulate the dataset under the Wright-Fisher diffusion
model <- "WFD"
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_gen <- -200
int_ale_frq <- 1 / (2 * pop_siz)
smp_gen <- (0:5) * 100
smp_chr_cnt <- rep(100, 6)
ptn_num <- 1e+01

sim_HMM_WFD <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_gen, int_ale_frq, smp_gen, smp_chr_cnt, ptn_num)
int_gen <- sim_HMM_WFD$int_gen
smp_gen <- sim_HMM_WFD$smp_gen
smp_chr_cnt <- sim_HMM_WFD$smp_chr_cnt
smp_ale_cnt <- sim_HMM_WFD$smp_ale_cnt
smp_ale_frq <- sim_HMM_WFD$smp_ale_cnt / sim_HMM_WFD$smp_chr_cnt
pop_ale_frq <- sim_HMM_WFD$pop_ale_frq

k <- min(int_gen, smp_gen):max(int_gen, smp_gen)
plot(k, pop_ale_frq, type = 'l', lwd = 1.5, 
     xlim = c(min(int_gen, smp_gen), max(int_gen, smp_gen)), ylim = c(min(smp_ale_frq, pop_ale_frq), max(smp_ale_frq, pop_ale_frq)), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the mutant allele generated with the Wright-Fisher diffusion")
points(smp_gen, smp_ale_frq, col = 'red', pch = 17, cex = 1)

################################################################################

#' Generate a simulated dataset under the Wright-Fisher model
test_seed <- 28
set.seed(test_seed)

model <- "WFM"
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_gen <- -200
int_ale_frq <- 1 / (2 * pop_siz)
smp_gen <- (0:5) * 100
smp_chr_cnt <- rep(100, 6)

sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_gen, int_ale_frq, smp_gen, smp_chr_cnt)
int_gen <- sim_HMM_WFM$int_gen
smp_gen <- sim_HMM_WFM$smp_gen
smp_chr_cnt <- sim_HMM_WFM$smp_chr_cnt
smp_ale_cnt <- sim_HMM_WFM$smp_ale_cnt
smp_ale_frq <- sim_HMM_WFM$smp_ale_cnt / sim_HMM_WFM$smp_chr_cnt
pop_ale_frq <- sim_HMM_WFM$pop_ale_frq

save(sel_cof, dom_par, pop_siz, int_gen, smp_gen, smp_chr_cnt, smp_ale_cnt, smp_ale_frq, pop_ale_frq, 
     file = "TEST_1L_simulated_dataset.rda")

load("TEST_1L_simulated_dataset.rda")

pdf(file = "TEST_1L_simulated_dataset.pdf", width = 10, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
k <- min(int_gen, smp_gen):max(int_gen, smp_gen)
plot(k, pop_ale_frq, type = 'l', lwd = 1.5, 
     xlim = c(min(int_gen, smp_gen), max(int_gen, smp_gen)), ylim = c(min(smp_ale_frq, pop_ale_frq), max(smp_ale_frq, pop_ale_frq)), 
     xlab = "Generation", ylab = "Allele frequency", main = "A simulated dataset generated with the W-F model")
points(smp_gen, smp_ale_frq, col = 'red', pch = 17, cex = 1)
dev.off()

########################################



################################################################################
