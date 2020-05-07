#!/usr/bin/env Rscript
source('lkhd_funcs.R')

if (!interactive()) {
  args = commandArgs(trailingOnly=TRUE)
  i_params = as.numeric(args[1]); i_sample = as.numeric(args[2]); n_chunks = as.numeric(args[3])
  i_chunk = as.numeric(args[4])
  cat('i_params i_sample n_chunks i_chunk =', i_params, i_sample, n_chunks, i_chunk, '\n')
} else {
  i_params = 223; i_sample = 19; n_chunks = 100; i_chunk = 85
}

dt=0.0001
base_dir = 'demography/samples8'
load(paste(base_dir,'/sample',i_params,'.RData',sep=''))
if (i_params==0) {
  y2 = 138000; ss = seq(0, 0.004, length.out=201); T2 = round(y2/8/2/N/dt)
} else if (i_params==1) {
  y2 = 56000; ss = seq(0.006, 0.019, length.out=201); T2 = round(y2/8/2/N/dt)
} else {
  T2 = 2000; ss = seq(0, 0.02, length.out=601)
}
#i_params = 1; i_sample = 1; alphas = seq(100,300)
#i_params = 1; i_sample = 2; alphas = seq(300,500)
#i_params = 1; i_sample = 3; alphas = seq(600,1000,2)

data4 = process_data_demog(counts, params_data, i_sample)
cat("data: ", data4[1,], "\n")
cat(data4[3,], "\n")
grid_option = 2; xs=gen_grid(grid_option)
#alphas = seq(0,150,1)
Nt = params_data$Nt
n1 = length(Nt$t)
alphas = ss * Nt$N[1] *2;
ts = (1:T2)*dt
lkhd_surface = matrix(0, nrow=length(alphas), ncol=T2)
for (i_alpha in seq(i_chunk, length(alphas), n_chunks)) {
  cat('alpha=', alphas[i_alpha], '\n')
  alpha_demog = list('alpha'=alphas[i_alpha],
                     't'=c(Nt$t[n1]-rev(Nt$t[2:n1])+1, Inf)/Nt$N[1]/2,
                     'rho'=c(rev(Nt$N[1:(n1-1)])/Nt$N[1], 1))
  #out1 = calc_lkhd_old(xs, 0.5, alphas[i_alpha], data4, T2, 'boundary', TRUE, 500000, 64, FALSE)
  out1 = calc_lkhd_demog(xs, 0.5, alpha_demog, data4, T2, 'boundary',
                         TRUE, 500000, 64, FALSE)
  lkhd_surface[i_alpha,] = out1$lkhd1_t
  cat(max(lkhd_surface[i_alpha,], "\n"))
}
#if (!interactive()) {
  file_name = file=paste(base_dir,'/lkhd','_',
    i_params, '_', i_sample, '_', i_chunk,'.RData', sep='')
  cat("saving to", file_name, "\n")
  save(lkhd_surface, alphas, ts, i_chunk, file=file_name)
#}
