#!/usr/bin/env Rscript
source('lkhd_funcs.R')

if (!interactive()) {
  args = commandArgs(trailingOnly=TRUE)
  i_params = as.numeric(args[1]); n_chunks = as.numeric(args[2]); i_chunk = as.numeric(args[3])
  mode_h = as.numeric(args[4]); base_dir = args[5]; grid_option = as.numeric(args[6])
  mode_age = as.numeric(args[7])
  x1_option = 1001
  cat(i_params, n_chunks, i_chunk, mode_h, base_dir, '\n')
} else {
  i_params = 1; n_chunks = 1; i_chunk = 1; mode_h = 0; base_dir = 'samples0'
  grid_option = 2; mode_age = 1
}
alphas0 = seq(-200, 1000, length.out=7)
load(paste(base_dir,'/sample',i_params,'.RData',sep=''))
#xs=gen_grid(grid_option); cat('grid points: ', 1/xs[1:5], '\n')
#data4 = process_data(counts, params_data, i_chunk)
#out1 = calc_lkhd_outer(xs, 1, -300, data4, 'age', 'boundary', TRUE, 5, 4)

section_search = function(xs, nt_factor_max, h, alphas, mode_age, y_in, N) {
  alpha_ceil = 0.04*N*2; alpha_floor = -0.02*N*2
  l = length(alphas); d_alpha = alphas[2]-alphas[1]
  y1 = matrix(NA, nrow=l, ncol=2)
  y1[1,] = y_in[1, c(4,2)]; y1[(l+1)/2,] = y_in[2, c(4,2)]; y1[l,] = y_in[3, c(4,2)]
  for (i_alpha in 1:l) {
    if (is.na(y1[i_alpha,1])) {
      out1 = calc_lkhd_outer(xs, h, alphas[i_alpha], data4,
                             mode_age, 'boundary', TRUE, 5, 4, nt_factor_max)
      y1[i_alpha, ] = c(out1$lkhd1, out1$t0_max1)
    }
    cat('alpha=', alphas[i_alpha], ', lkhd=', y1[i_alpha,1], ', t0_max=', y1[i_alpha,2], '\n')
  }
  i_max = which.max(y1[,1]); alpha_max = alphas[i_max]
  if (i_max>1 && i_max<l) {d_alpha = d_alpha/3}
  alphas_out = seq(alpha_max-3*d_alpha, alpha_max+3*d_alpha, length.out=l)
  if (i_max==1) {
    if (min(alphas_out)<=alpha_floor) {
      alphas_out = seq(alpha_floor, alphas_out[l], length.out=l)
      y_out = rbind(c(NA,NA,NA,NA), c(NA,NA,NA,NA),
                    c(alphas[i_max+3], y1[i_max+3,2], h, y1[i_max+3,1]))
    } else {
      y_out = rbind(c(NA,NA,NA,NA), c(alphas[i_max], y1[i_max,2], h, y1[i_max,1]),
                    c(alphas[i_max+3], y1[i_max+3,2], h, y1[i_max+3,1]))
    }
  } else if (i_max==l) {
    if (max(alphas_out)>=alpha_ceil) {
      alphas_out = seq(alphas_out[1], alpha_ceil, length.out=l)
      y_out = rbind(c(alphas[i_max-3], y1[i_max-3,2], h, y1[i_max-3,1]),
                    c(NA,NA,NA,NA), c(NA,NA,NA,NA))
    } else {
      y_out = rbind(c(alphas[i_max-3], y1[i_max-3,2], h, y1[i_max-3,1]),
                  c(alphas[i_max], y1[i_max,2], h, y1[i_max,1]), c(NA,NA,NA,NA))
    }
  } else {
    y_out = rbind(c(alphas[i_max-1], y1[i_max-1,2], h, y1[i_max-1,1]),
                  c(alphas[i_max], y1[i_max,2], h, y1[i_max,1]),
                  c(alphas[i_max+1], y1[i_max+1,2], h, y1[i_max+1,1]))
  }
  return(list(alphas=alphas_out, y=y_out, y1=y1, i_max=i_max, alphas_old=alphas))
}

n_reals = min(200,dim(params_data)[1])
out = matrix(NA, nrow=n_reals, ncol=8)
for (i_sample in seq(i_chunk, n_reals, n_chunks)) {
  cat(i_sample, i_chunk, n_reals, n_chunks, '\n')
  data4 = process_data(counts, params_data, i_sample)
  cat('data: ', data4[3,], '\n')
  if (mode_h==0) {
    hs = params_data[i_sample,3]
  } else {
    hs = c(0, 0.5, 1)
  }
  for (grid_option1 in grid_option:1) {
    xs=gen_grid(grid_option1)
    cat('##### grid option ', grid_option1, ', points: ', 1/xs[1:5], '\n')
    if (grid_option1>1) nt_factor_max = 50
    else nt_factor_max = 200
    out2 = matrix(NA, nrow=length(hs), ncol=4)
    next_grid = FALSE
    for (i_h in 1:length(hs)) {
      h = hs[i_h]
      alphas1 = alphas0
      y5 = matrix(NA, nrow=3, ncol=4)
      for (i1 in 1:5) {
        cat('alphas1: ', alphas1, '\n')
        y3 = section_search(xs, nt_factor_max, h, alphas1, 0, y5, params_data[i_sample,1])
        if (max(y3$y[,4], na.rm=T)==0) {next_grid=TRUE; break}
        y5 = y3$y
        alphas1 = y3$alphas
        if (alphas1[2]-alphas1[1]<100) break
      }
      if (next_grid) break
      alphas2 = alphas1
      y5 = matrix(NA, nrow=3, ncol=4)
      for (i1 in 1:6) {
        cat('alphas2: ', alphas2, '\n')
        y4 = section_search(xs, nt_factor_max, h, alphas2, mode_age, y5, params_data[i_sample,1])
        if (max(y4$y[,4], na.rm=T)==0) {next_grid=TRUE; break}
        y5 = y4$y
        alphas2 = y4$alphas
        if (alphas2[length(alphas2)]-alphas2[1]<3) break
      }
      if (next_grid) break
      out2[i_h,] = c(y4$alphas_old[y4$i_max], y4$y1[y4$i_max,2], h, y4$y1[y4$i_max,1])
    }
    if (!next_grid) break
  }
  out[i_sample,] = c(out2[which.max(out2[,4]),], out2[out2[,3]==params_data[i_sample,3],])
}

if (interactive()) {
  cat(out[i_chunk,1:4])
  save(out, file=paste(base_dir,'/out',grid_option,'_',mode_age,'_',
                       i_params, '_', mode_h, '.RData', sep=''))
} else {
  if (n_chunks>1) {
    cat('output:', out[i_chunk,1:4], '\n')
    save(out, file=paste(base_dir,'/out',grid_option,'_',mode_age,'_',
                         i_params, '_', mode_h, '_', i_chunk,'.RData', sep=''))
  } else {
    save(out, file=paste(base_dir,'/out',grid_option,'_',mode_age,'_',
                         i_params, '_', mode_h, '.RData', sep=''))
  }
}
