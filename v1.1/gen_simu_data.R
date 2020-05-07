source('simu_funcs.R')

# params_data: matrix, each row giving:
# (population size N, number of sample points, h, s, allele age in diff time)
N=16000; dt=0.0001; m=100; grid_option=2; mode_age=1

task = 'show_estimates' # compare_data, gen_data, show_data, show_estimates
task_type='wf' # wf or diffusion
# mode2=1: fix sampling times, =2: find appropriate sampling time according to s and h
#samples_id = 1; mode2=2; f2=0.9; f1=0.1; smsz=60; n_ts=10; f0=0
#samples_id = 2; mode2=2; f2=0.9; f1=0.1; smsz=20; n_ts=60; f0=0
#samples_id = 3; mode2=2; f2=0.9; f1=0.1; smsz=10; n_ts=60; f0=0
#samples_id = 4; mode2=2; f2=0.9; f1=0.1; smsz=15; n_ts=60; f0=0
#samples_id = 5; mode2=2; f2=0.9; f1=0.1; smsz=120; n_ts=10; f0=0
#samples_id = 6; mode2=2; f2=0.9; f1=0.1; smsz=60; n_ts=10; f0=0
#samples_id = 7; mode2=2; f2=0.9; f1=0.1; smsz=90; n_ts=10; f0=0
#samples_id = 9; m=400
#samples_id = 11; m=200
#mode2=1; f0=0
#samples_id = 21; mode2=1; T1=2000; T2=5000; smsz=20; n_ts=60; f0=0
#samples_id = 22; mode2=1; T1=2000; T2=5000; smsz=120; n_ts=10; f0=0
samples_id = 23; mode2=1; T1=2000; T2=5000; smsz=60; n_ts=10; f0=0

base_dir = paste('samples', samples_id, '/', sep='')
par(mfrow=c(3,4))

if (samples_id==21 || samples_id==22) {
  s_list = 0; h_list = 0.5; last_smp_gbp=0; smp_siz = list(random=0, val=smsz)
} else if (samples_id==23) {
  s_list = 0; h_list = 0.5; last_smp_gbp=0; smp_siz = list(random=3, val=smsz)
} else if (samples_id==8 || samples_id==10 || samples_id==12 || samples_id==14) {
  load('samples0/sample0.RData')
  counts0 = counts[,,2]; params_data0 = params_data[2,]
  n_ts = params_data0[2]; s_list = 56.79012/N/2; h_list = params_data0[3]; last_smp_gbp=params_data0[7]
  if (samples_id==8) {
    T2 = -counts0[1,1]; T1=round(0.2293203*N*2-T2)
    smp_siz = list(random=2, smp_gen=counts0[1,], smsz=counts0[2,])
  } else if (samples_id==10) {
    T2 = 11311; T1 = 2814
    smp_siz = list(random=2, smp_gen=round(-counts0[1,]/counts0[1,1]*T2), smsz=counts0[2,])
  } else if (samples_id==12){
    smp_siz = list(random=2, smp_gen=counts0[1,], smsz=rep(10, n_ts))
  } else {
    smp_siz = list(random=2, smp_gen=counts0[1,], smsz=rep(20, n_ts))
  }
} else if (samples_id==9 || samples_id==11 || samples_id==13 || samples_id==15 ||
           samples_id==17) {
  load('samples0/sample1.RData')
  counts0 = counts[,,2]; params_data0 = params_data[2,]
  n_ts = params_data0[2]; s_list = 413.1687/N/2; h_list = params_data0[3]; last_smp_gbp=params_data0[7]
  if (samples_id==9) {
    T2 = -counts0[1,1]; T1=round(0.08122031*N*2-T2)
    smp_siz = list(random=2, smp_gen=counts0[1,], smsz=counts0[2,])
  } else if (samples_id==11) {
    T2 = 2905; T1 = 1739
    smp_siz = list(random=2, smp_gen=round(-counts0[1,]/counts0[1,1]*T2), smsz=rep(5, n_ts))
  } else if (samples_id==13) {
    smp_siz = list(random=2, smp_gen=counts[1,,2], smsz=rep(10, n_ts))
  } else if (samples_id==15) {
    smp_siz = list(random=2, smp_gen=counts[1,,2], smsz=rep(20, n_ts))
  } else {
    smp_siz = list(random=2, smp_gen=counts[1,,2], smsz=rep(50, n_ts))
  }
} else {
  last_smp_gbp=0; s_list = c(0.0025, 0.005, 0.01, 0.015, 0.02); h_list = c(0, 0.5, 1)
  smp_siz = list(random=0, val=smsz);
}

for (i_h in 1:length(h_list)) {
  for (i_s in 1:length(s_list)) {
    s = s_list[i_s]; h = h_list[i_h]
    #if (s==0 && h!=0.5) next
    if (task_type=='wf') {
      i_params = 200 + i_h*10 + i_s
    } else {
      i_params = 100 + i_h*10 + i_s
    }
    file_name = paste(base_dir, 'sample',i_params,'.RData',sep='')
    if (task=='gen_data') {
      if (task_type=='wf') {
        if (mode2==2) {
          if (samples_id==2 || samples_id==3 || samples_id==4) {
            out_ts = calc_avg_path(f1, f2, h, s)
            out = gen_samples_wf(out_ts$T1, out_ts$T2, n_ts, smp_siz, 
                                 N, s, h, m, f0, last_smp_gbp)
          } else if (samples_id==5 || samples_id==6 || samples_id==7) {
            load(paste('samples', samples_id-3, '/sample',i_params,'.RData',sep=''))
            out = gen_samples_wf(params_data[1,6], params_data[1,5]-params_data[1,6],
                                 n_ts, smp_siz, N, s, h, m, f0, last_smp_gbp)
          } else if (samples_id==1) {
            load(paste('samples5/sample',i_params,'.RData',sep=''))
            smp_siz = list(random=3, val=smsz);
            out = gen_samples_wf(params_data[1,6], params_data[1,5]-params_data[1,6],
                                 n_ts, smp_siz, N, s, h, m, f0, last_smp_gbp)
          }
        } else {
          out = gen_samples_wf(T1, T2, n_ts, smp_siz, N, s, h, m, f0, last_smp_gbp)
        }
      } else {
        s_info=c(s,s)
        out = gen_samples_diffusion(t0, t1, n_ts, smp_siz, s_info, h, m=m)
      }
      counts = out$counts; params_data = out$params_data; traj_list = out$traj_list
      n_discarded = out$n_discarded
      save(counts, params_data, traj_list, n_discarded, file=file_name)
    } else if (task=='show_data') {
      if (file.exists(file_name)) load(file_name)
      else cat(file_name, 'does not exist\n')
    }
    if (task=='gen_data' || task=='show_data') {
      if (TRUE || s!=0) {
        matplot(counts[3,,seq(1, dim(params_data)[1],5)]/counts[2,,seq(1, dim(params_data)[1],5)],
                type='l', ylab='',
                main=paste('file_id=', i_params, ', h=', h, ', s=', s, '\n', sep=''))
        points(rowMeans(counts[3,,]/counts[2,,]), pch=21, col='red', bg='red', cex=1)
      }
      cat('file_id=', i_params, ', h=', h, ', s=', s, ', num_all_one_paths = ',
          n_discarded, ', sampling times 1st last =', params_data[1,6], params_data[1,5], '\n')
      next
    }
    if (task=='show_estimates') {
      res_file_name = paste(base_dir, 'res', grid_option, '_', mode_age, '_', i_params,'_0.RData',sep='')
      if (file.exists(res_file_name)) {
        load(res_file_name)
        load(file_name)
        s_true = params_data[1,4]; t1_true = params_data[1,5]
        ii = which(apply(counts[3,,]/counts[2,,], 2, max)>0.1)
        cat('h=', h, ', s=', s, ', bias sel_cof =',
            median(res[ii,1], na.rm=T)/2/N-s_true, '(',
            round((median(res[ii,1], na.rm=T)/2/N-s_true)/s_true*100), '%), bias age=',
            median(res[ii,2]*2*N), '(',
            round((median(res[ii,2])*2*N-t1_true)/t1_true*100), '%), #NA=', sum(is.na(res[ii,1])), '\n')
        if (s!=0 || h==0.5) {
          boxplot(res[ii,1]/2/N, main = paste('h=', h, ', s=', s, '\n', sep=''))
          ogn_pt = 1 - 0.4; end_pt = 1 + 0.4
          segments(x0 = ogn_pt, x1 = end_pt, y0 = s_true, col = "red", lty = 2, lwd = 2)
          if (!is.na(res[1,2])) {
            boxplot(res[ii,2]*2*N, main = paste('h=', h, ', s=', s, '\n', sep=''))
            ogn_pt = 1 - 0.4; end_pt = 1 + 0.4
            segments(x0 = ogn_pt, x1 = end_pt, y0 = t1_true, col = "red", lty = 2, lwd = 2)
          }
        }
      }
      next
    }
    if (task=='compare_data') {
      compare_sample_data(i_h*10 + i_s, base_dir)
    }
  }
}

#load('res224_0.RData')
#load('samples/sample224.RData')
#summary(res[,1])
#xmin=apply(counts[3,,], 2, min); xmax=apply(counts[3,,], 2, max)
#ii = which(xmax/xmin>1.4)
#summary(res[ii,1])

if (FALSE) {
  x = matrix(0, nrow=1000, ncol=9)
  y1 = rep(NA, 1000)
  y2 = rep(NA, 1000)
  for (i in 1:1000) {
    x[i,] = counts[3,,i]/counts[2,,i]
    y1[i] = sum(counts[2,,i])
    y2[i] = sum(counts[3,,i])
  }
}


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
