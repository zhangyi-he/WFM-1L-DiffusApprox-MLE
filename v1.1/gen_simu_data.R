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
