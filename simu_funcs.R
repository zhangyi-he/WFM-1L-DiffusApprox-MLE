gen_ts_sample = function(t0, t1, n_ts, ts_sample_in=NA) {
  if (!is.na(ts_sample_in)) {return(ts_sample_in)}
  ts_sample = round(seq(t0,t1,length.out=n_ts)/dt)*dt
  return(ts_sample)
}

# m=(how many paths are needed?), n=(how many to simulate every run)
gen_samples_diffusion = function(t0, t1, n_ts, smp_siz, s_info, h_info, m=100, n=20000) {
  cat('s_info:', s_info, ', h:', h_info, '\n')
  T1 = round(t1/dt);
  ts_sample = gen_ts_sample(t0,t1,n_ts)
  traj = matrix(0, nrow=m, ncol=T1)
  params_data = matrix(0, nrow=m, ncol=5)
  counts = array(NA, dim=c(3,n_ts,m))
  k = 0 # keep track of kept paths
  n_all_ones = 0 # keep track of all 1 paths discarded
  flag_stop = FALSE
  for (j in 1:30) {
    ss = runif(n, min=s_info[1], max=s_info[2]); alphas = 2*N*ss
    hs = sample(h_info, n, replace=T)
    # flag indicating whether change in selection has occurred, 1 means change yet to occur
    flag_switch = rep(1, n)
    x = matrix(0, nrow=n, ncol=T1)
    x[,1] = 1/N/2
    for (t in 2:T1) {
      if (t%%100==0) {
        i_alive = which(x[,t-1]>0) # & x[,t-1]<1)
        if (length(i_alive)>=2) {
          x = x[i_alive,]; hs = hs[i_alive]; ss = ss[i_alive]; alphas = alphas[i_alive]
        }
      }
      y1 = x[,t-1] * (1-x[,t-1]) * dt
      x[,t] = pmin(1,pmax(0, x[,t-1] + alphas*y1*(1-hs-(1-2*hs)*x[,t-1])
                          + sqrt(y1)*rnorm(dim(x)[1])))
    }
    cat('j=', j, dim(x)[1], 'paths alive\n')
    if (dim(x)[1]<=1) {next}
    for (i in 1:dim(x)[1]) {
      if (x[i,round(ts_sample[n_ts]/dt)]==0) {next}
      if (smp_siz$random==0) {smsz = rep(smp_siz$val, n_ts)}
      else if (smp_siz$random==1) {smsz = rpois(n_ts,smp_siz$val)+1}
      else {smsz = smp_siz$val}
      mutant_counts_true = mapply(function(smsz1, p) 
        rbinom(1, smsz1, p), smsz, x[i,round(ts_sample/dt)])
      if (min(smsz-mutant_counts_true)<0) {cat('something bad')}
      if (max(mutant_counts_true/smsz)>0) {
        if (max(smsz-mutant_counts_true)==0) {
          n_all_ones = n_all_ones + 1
        } else {
          k = k + 1
          traj[k,] = x[i,]
          params_data[k,] = c(N,n_ts,hs[i],ss[i],t1)
          counts[1,,k] = ts_sample
          counts[2,,k] = smsz
          counts[3,,k] = mutant_counts_true
          if (k>=m) {flag_stop=TRUE; break}
        }
      }
    }
    if (flag_stop) {break}
  }
  cat('num all 1 paths = ', n_all_ones, '\n')
  return(list(counts=counts, params_data=params_data, traj=traj, n_discarded=n_all_ones))
}

simulate_wf = function(T1, T2, N, s, h, f0=0) {
  if (mode2==3) {max_len = T2*10; T1s=NA}
  else {max_len = T1+T2; T1s=T1}
  frq_pth = rep(NA, max_len); frq_pth[1] = f0
  fitness = c(1, 1-s*h, 1-s)
  for (k in 1:(max_len-1)) {
    if (is.na(T1s) && T1<1 && frq_pth[k]>=T1) {
      if (k+T2>max_len) {break}
      T1s=k
    }
    ff = frq_pth[k]
    gen_frq = c(ff*ff, 2*ff*(1-ff), (1-ff)*(1-ff))
    gen_frq = fitness * gen_frq
    gen_frq = gen_frq / sum(gen_frq)
    ale_frq = gen_frq[1] + gen_frq[2]/2
    frq_pth[k+1] = rbinom(1, size=2*N, prob=ale_frq)/2/N
    if (frq_pth[k+1]==0) {
      frq_pth[(k+1):max_len] = 0
      break}
    if (frq_pth[k+1]==1) {
      frq_pth[(k+1):max_len] = 1
      break
    }
    if (T1<1 && !is.na(T1s) && k>T1s+T2) break
  }
  return(list(frq_pth=frq_pth, T1s=T1s, T2s=T2))
}

simulate_hmm = function(T1, T2, n_ts, smsz, N, s, h, f0) {  
  n_all_ones = 0 # keep track of all 1 paths discarded
  repeat {
    out = simulate_wf(T1, T2, N, s, h, f0)
    if (!is.na(out$T1s) && (is.na(tail(out$frq_pth,1)) || tail(out$frq_pth,1)>0)) {
      smp_cnt = numeric(n_ts)
      smp_gen = round(seq(out$T1s, out$T1s+out$T2s, length.out=n_ts))
      for (k in 1:n_ts) {
        smp_cnt[k] = rbinom(1, size=smsz[k], prob=out$frq_pth[smp_gen[k]])
      }
      if (max(smp_cnt)>=1) {
        if (max(smsz-smp_cnt)==0) {
          n_all_ones = n_all_ones + 1
        } else {
          break
        }
      }
    }
  }
  return(list(smp_gen = smp_gen, smp_siz = smsz, smp_cnt = smp_cnt,
              pop_frq = out$frq_pth, n_discarded = n_all_ones))
}

gen_samples_wf = function(T1, T2, n_ts, smp_siz, N, s, h, m=100, f0=0) {
  traj_list = list()
  #traj = matrix(0, nrow=m, ncol=T1)
  params_data = matrix(0, nrow=m, ncol=6)
  counts = array(NA, dim=c(3,n_ts,m))
  n_all_ones = 0
  for (i in 1:m) {
    if (smp_siz$random==0) {smsz = rep(smp_siz$val, n_ts)}
    else if (smp_siz$random==1) {smsz = rpois(n_ts,smp_siz$val)+1}
    else {smsz = smp_siz$val}
    if (f0==0) {
      out = simulate_hmm(T1, T2, n_ts, smsz, N, s, h, 1/2/N)
    } else {
      out = simulate_hmm(T1, T2, n_ts, smsz, N, s, h, f0)
    }
    n_all_ones = n_all_ones + out$n_discarded
    counts[1,,i] = out$smp_gen
    counts[2,,i] = out$smp_siz
    counts[3,,i] = out$smp_cnt
    params_data[i,] = c(N,n_ts,h,s,out$smp_gen[n_ts],out$smp_gen[1])
    traj_list[[i]] = out$pop_frq[1:out$smp_gen[n_ts]]
    if (i%%10==0) cat(i, out$n_discarded, n_all_ones, ' || ')
  }
  cat('\n', 'num all 1 paths discarded = ', n_all_ones, '\n')
  return(list(counts=counts, params_data=params_data, traj_list=traj_list, n_discarded=n_all_ones))
}

compare_sample_data = function(jp, base_dir) {
  load(paste(base_dir, 'sample2',jp,'.RData',sep='')); traj1=traj
  load(paste(base_dir, 'sample1',jp,'.RData',sep=''));
  #n=2000; n1=round(n/3.2); cat(mean(traj1[,n1]), mean(traj[,n]), '\n')
  plot((1:ncol(traj1))/ncol(traj1), colMeans(traj1), type='l',
       main=paste(jp, 'red is diffusion'))
  lines((1:ncol(traj))/ncol(traj), colMeans(traj), col='red')
}
