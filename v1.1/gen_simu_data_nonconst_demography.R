sim_wf = function(Nt, s, h) {
  fitness = c(1, 1-s*h, 1-s)
  len_Nt = length(Nt$t)
  max_len = Nt$t[len_Nt]
  repeat {
    x = rep(NA, max_len)
    x[1] = 1/2/Nt$N[1]
    N = Nt$N[1]
    i_Nt = 1
    for (k in 1:(max_len-1)) {
      if (k==Nt$t[i_Nt+1]) {
        i_Nt = i_Nt+1; N = Nt$N[i_Nt]; #cat("N change to", N, " at ", k, "\n")
      }
      ff = x[k]
      geno_freq = fitness * c(ff*ff, 2*ff*(1-ff), (1-ff)*(1-ff))
      geno_freq = geno_freq / sum(geno_freq)
      ale_freq = geno_freq[1] + geno_freq[2]/2
      x[k+1] = rbinom(1, size=2*N, prob=ale_freq)/2/N
      if (x[k+1]==0) {x[(k+1):max_len] = 0; break}
      if (x[k+1]==1) {x[(k+1):max_len] = 1; break}
    }
    if (tail(x, 1)>0) break
  }
  return(x)
}

gen_wf_data = function(Nt, s, h, smp_size1, smp_t1, smp_size2, smp_t2) {
  n_ts1 = length(smp_size1)
  n_ts2 = length(smp_size2)
  max_len = tail(Nt$t, 1)
  traj = matrix(NA, nrow=max_len, ncol=m)
  counts1 = array(NA, dim=c(3,n_ts1,m))
  counts2 = array(NA, dim=c(3,n_ts2,m))
  for (i in 1:m) {
    repeat {
      x = sim_wf(Nt, s, h)
      smp_cnt1 = numeric(n_ts1)
      for (k in 1:n_ts1) {
        smp_cnt1[k] = rbinom(1, size=smp_size1[k], prob=x[smp_t1[k]])
      }
      smp_cnt2 = numeric(n_ts2)
      for (k in 1:n_ts2) {
        smp_cnt2[k] = rbinom(1, size=smp_size2[k], prob=x[smp_t2[k]])
      }
      if ((max(smp_cnt1)>=1) && (max(smp_size1 - smp_cnt1)>0)
          && (max(smp_cnt2)>=1) && (max(smp_size2 - smp_cnt2)>0)) {break}
    }
    traj[,i] = x
    counts1[1, ,i] = smp_t1
    counts1[2, ,i] = smp_size1
    counts1[3, ,i] = smp_cnt1
    counts2[1, ,i] = smp_t2
    counts2[2, ,i] = smp_size2
    counts2[3, ,i] = smp_cnt2
    if (i %% 20 == 0) {cat(i, "samples generated\n")}
  }
  params = list('Nt'=Nt, 's'=s, 'h'=h, 'k_last'=max_len, 'last_smp_gbp'=0)
  return(list('counts1'=counts1, 'counts2'=counts2, 'params'=params, 'traj'=traj))
}

par(mfrow=c(2,3))
s_list = c(0.0025, 0.005, 0.01, 0.015, 0.02, 0)
h = 0.5
m = 120 # number of replicates desired
i_base_dir = 9
for (i in 1:length(s_list)) {
  if (i<=5) {
    old_data_file = sprintf('../samples5/sample22%d.RData', i)
    load(old_data_file)
    k_first = params_data[1,6]; k_last=params_data[1,5]
  } else {
    k_first = 2000; k_last=7000
  }
  cat("k_first, last=", k_first, " ", k_last, "\n")
  i_params = 220 + i
  cat("i_params=", i_params, "\n")
  if (i_base_dir==1) {
    Nt = list('t'=round(c(1, k_last/2, k_last/4*3, k_last)),
              'N'=c(16000, 8000, 32000, 32000))
  } else if (i_base_dir==3) {
    Nt = list('t'=round(c(1, k_last/2, k_last/4*3, k_last)),
              'N'=c(16000, 8000, 16000, 16000))
  } else if (i_base_dir==5) {
    Nt = list('t'=round(c(1, k_last/2, k_last/4*3, k_last)),
                'N'=c(8000, 16000, 16000, 16000))
  } else if (i_base_dir==7) {
    Nt = list('t'=round(c(1, k_last/2, k_last/4*3, k_last)),
              'N'=c(16000, 32000, 32000, 32000))
  } else if (i_base_dir==9) {
    Nt = list('t'=round(c(1, k_last/2, k_last/4*3, k_last)),
              'N'=c(100000, 16000, 16000, 16000))
  } else if (i_base_dir==11) {
    Nt = list('t'=round(c(1, k_last)),
              'N'=c(16000, 16000))
  }
  s = s_list[i]
  n_ts1=10; smp_size1=rep(120,n_ts1); smp_t1=round(seq(k_first,k_last,length.out=n_ts1))
  n_ts2=60; smp_size2=rep(20,n_ts2); smp_t2=round(seq(k_first,k_last,length.out=n_ts2))
  out = gen_wf_data(Nt, s, h, smp_size1, smp_t1, smp_size2, smp_t2)
  counts = out$counts1
  traj = out$traj
  params_data = out$params; if (i_base_dir==9) params_data$N0=16000
  matplot(counts[3,,seq(1, m,5)]/counts[2,,seq(1, m,5)], type='l', ylab='',
          main=paste(i_params, ', ', 's=', s, '\n', sep=''))
  points(rowMeans(counts[3,,]/counts[2,,]), pch=21, col='red', bg='red', cex=1)
  file_name = sprintf('samples%d/sample%d.RData', i_base_dir, i_params)
  save(counts, params_data, traj, file=file_name)
  counts = out$counts2
  traj = out$traj
  params_data = out$params; if (i_base_dir==9) params_data$N0=16000
  matplot(counts[3,,seq(1, m,5)]/counts[2,,seq(1, m,5)], type='l', ylab='',
          main=paste(i_params, ', ', 's=', s, '\n', sep=''))
  points(rowMeans(counts[3,,]/counts[2,,]), pch=21, col='red', bg='red', cex=1)
  file_name = sprintf('samples%d/sample%d.RData', i_base_dir+1, i_params)
  save(counts, params_data, traj, file=file_name)
}

par(mfrow=c(2,3))
s_list = c(0.0025, 0.005, 0.01, 0.015, 0.02, 0)
h = 0.5
m = 120 # number of replicates desired
for (i_base_dir in 1:8) {
  for (i_params in c(221,222,223,224,225,226)) {
    load(sprintf('samples%d/sample%d.RData', i_base_dir, i_params))
    N0=params_data$Nt$N[1]; s=params_data$s; k_last=params_data$k_last
    t0_true = k_last/2/N0
    load(sprintf('samples%d/res2_1_%d_0_0.RData', i_base_dir, i_params))
    s1=median(res[,1])/2/N0; t1=median(res[,2])
    load(sprintf('samples%d/res2_1_%d_0_%d.RData', i_base_dir, i_params, N0))
    s2=median(res[,1])/2/N0; t2=median(res[,2])
    cat(sprintf('base_dir=%d %d, s=%0.3f %0.3f %0.3f (%0.1f %0.1f), t0=%0.3f %0.3f %0.3f (%0.1f %0.1f)\n',
          i_base_dir, i_params, s, s1, s2, (s1-s)/s*100, (s2-s)/s*100,
          t0_true, t1, t2, (t1-t0_true)/t0_true*100, (t2-t0_true)/t0_true*100))
    if (i_base_dir==3 || i_base_dir==4) {
      N0 = 8000
      load(sprintf('samples%d/res2_1_%d_0_%d.RData', i_base_dir, i_params, N0))
      s3=median(res[,1])/2/N0; t3=median(res[,2])
      cat(sprintf('s=%0.3f (%0.1f), t0=%0.3f (%0.1f)\n',
                  s3, (s3-s)/s*100, t3, (t3-t0_true)/t0_true*100))
    }
    if (i_base_dir==7 || i_base_dir==8 ) {
      N0 = 32000
      load(sprintf('samples%d/res2_1_%d_0_%d.RData', i_base_dir, i_params, N0))
      s3=median(res[,1])/2/N0; t3=median(res[,2])
      cat(sprintf('s=%0.3f (%0.1f), t0=%0.3f (%0.1f)\n',
                  s3, (s3-s)/s*100, t3, (t3-t0_true)/t0_true*100))
    }
    if (i_base_dir==9 || i_base_dir==10 ) {
      N0 = 16000
      load(sprintf('samples%d/res2_1_%d_0_%d.RData', i_base_dir, i_params, N0))
      s3=median(res[,1]); t3=median(res[,2])/N0/2
      cat(sprintf('s=%0.3f (%0.1f), t0=%0.3f (%0.1f)\n',
                  s3, (s3-s)/s*100, t3, (t3-t0_true)/t0_true*100))
    }
  }
}
cat('base_dir, sample file, s=true, estimates using true demography, estimates constant pop size (percentage bias)')
for (i_base_dir in 9:10) {
  for (i_params in c(221,222,223,224,225,226)) {
    load(sprintf('samples%d/sample%d.RData', i_base_dir, i_params))
    N0=params_data$Nt$N[1]; s=params_data$s; k_last=params_data$k_last
    t0_true = k_last
    load(sprintf('samples%d/res2_1_%d_0_0.RData', i_base_dir, i_params))
    s1=median(res[,1]); t1=median(res[,2])
    load(sprintf('samples%d/res2_1_%d_0_%d.RData', i_base_dir, i_params, N0))
    s2=median(res[,1]); t2=median(res[,2])
    cat(sprintf('base_dir=%d %d, s=%0.3f %0.3f %0.3f (%0.1f %0.1f), t0=%0.3f %0.3f %0.3f (%0.1f %0.1f)\n',
                i_base_dir, i_params, s, s1, s2, (s1-s)/s*100, (s2-s)/s*100,
                t0_true, t1, t2, (t1-t0_true)/t0_true*100, (t2-t0_true)/t0_true*100))
    N0 = 16000
    load(sprintf('samples%d/res2_1_%d_0_%d.RData', i_base_dir, i_params, N0))
    s3=median(res[,1], na.rm=TRUE); t3=median(res[,2], na.rm=TRUE)
    cat(sprintf('s=%0.3f (%0.1f), t0=%0.3f (%0.1f)\n',
                s3, (s3-s)/s*100, t3, (t3-t0_true)/t0_true*100))

  }
}
