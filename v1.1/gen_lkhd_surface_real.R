#!/usr/bin/env Rscript
source('lkhd_funcs.R')

if (!interactive()) {
  args = commandArgs(trailingOnly=TRUE)
  i_params = as.numeric(args[1]); i_sample = as.numeric(args[2]); n_chunks = as.numeric(args[3])
  i_chunk = as.numeric(args[4]); mode_N = as.numeric(args[5])
  cat('i_params i_sample n_chunks i_chunk =', i_params, i_sample, n_chunks, i_chunk, '\n')
} else {
  i_params = 1; i_sample = 2; n_chunks = 400; i_chunk = 202; mode_N = 0
}

dt=0.0001
base_dir = 'samples0'
load(paste(base_dir,'/sample',i_params,'.RData',sep=''))
#Nt = list('t'=c(-34557,-27756,-22028,-17205,-13142,-9721,-1144,0),
#           'N'=c(171778,124078,82404,107519,195019,110000,16000,16000))
Nt2 = list('t'=c(-34557,-27756,-22028,-17205,-13142,-9721,-1144,0),
           'N'=c(171778,124078,82404,107519,195019,110000,16000,16000))
p11 = seq(-6644, -1145, 10)
N11 = 110/16*exp(-20.956/32000*(p11 + 6644)) * 16000
Nt = list('t'=c(Nt2$t[1:6], p11, Nt2$t[7:8]),
           'N'=c(Nt2$N[1:6], N11, Nt2$N[7:8]))
params_data = list('Nt'=Nt, 'N0'=16000, 's'=NULL, 'h'=params_data[1,3],
                   'k_last'=34557, 'last_smp_gbp'=params_data[1,7])
if (mode_N>0) {
  params_data$Nt = list('t'=c(1,params_data$k_last), 'N'=c(mode_N, mode_N))
  params_data$N0 = mode_N
  Nt = params_data$Nt
}

data4 = process_data_demog(counts, params_data, i_sample)
cat("data: ", data4[1,], "\n")
cat("data: ", data4[3,], "\n")
grid_option = 2; xs=gen_grid(grid_option)
#alphas = seq(0,150,1)
n1 = length(Nt$t)
if ("N0" %in% names(params_data)) {N0 = params_data$N0
} else {N0 = Nt$N[1]}
if (i_params==0) {
  ss = seq(0, 0.005, length.out=401); T2 = round(20000/2/N0/dt)
} else if (i_params==1) {
  ss = seq(0.005, 0.02, length.out=401); T2 = round(8000/2/N0/dt)
}
ts = (1:T2)*dt * N0 * 2
lkhd_surface = matrix(0, nrow=length(ss), ncol=T2)
for (i_alpha in seq(i_chunk, length(ss), n_chunks)) {
  cat('s=', ss[i_alpha], '\n')
  alpha_demog = list('alpha'=ss[i_alpha]*N0*2,
                     't'=c(Nt$t[n1]-rev(Nt$t[2:n1])+1, Inf)/N0/2,
                     'rho'=c(rev(Nt$N[1:(n1-1)])/N0, NA))
  #out1 = calc_lkhd_old(xs, 0.5, alphas[i_alpha], data4, T2, 'boundary', TRUE, 500000, 64, FALSE)
  out1 = calc_lkhd_demog(xs, params_data$h, alpha_demog, data4, T2, 'boundary',
                         TRUE, 500000, 100, FALSE)
  lkhd_surface[i_alpha,] = out1$lkhd1_t
  cat(max(lkhd_surface[i_alpha,]), "\n")
}
#if (!interactive()) {
  file_name = file=paste(base_dir,'/lkhd','_',
    i_params, '_', i_sample, '_', i_chunk,'.RData', sep='')
  cat("saving to", file_name, "\n")
  save(lkhd_surface, ss, ts, i_chunk, file=file_name)
#}
