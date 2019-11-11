source('simu_funcs.R')

# params_data: matrix, each row giving:
# (population size N, number of sample points, h, s, allele age in diff time)
N=16000; dt=0.0001; m=100; grid_option=2; mode_age=1

task = 'show_estimates' # compare_data, gen_data, show_data, show_estimates
task_type='wf' # wf or diffusion
if (task=='show_estimates') {par(mfrow=c(3,4))} else {par(mfrow=c(3,4))}
#base_dir = 'samples1/'; mode2=3; T2=1500; T1=0.05; smsz=20; n_ts=60; f0=0 # mode2=3 unprogrammed
#base_dir = 'samples2/'; mode2=2; f2=0.9; f1=0.1; smsz=20; n_ts=60; f0=0
#base_dir = 'samples3/'; mode2=2; f2=0.9; f1=0.1; smsz=10; n_ts=60; f0=0
#base_dir = 'samples4/'; mode2=2; f2=0.9; f1=0.1; smsz=15; n_ts=60; f0=0
base_dir = 'samples5/'; mode2=2; f2=0.9; f1=0.1; smsz=120; n_ts=10; f0=0
#base_dir = 'samples6/'; mode2=2; f2=0.9; f1=0.1; smsz=60; n_ts=10; f0=0
#base_dir = 'samples7/'; mode2=2; f2=0.9; f1=0.1; smsz=90; n_ts=10; f0=0
#base_dir = 'samples8/'; mode2=1; T1=2000; T2=4000; smsz=90; n_ts=10; f0=0
smp_siz = list(random=0, val=smsz);

calc_det_path = function(f1, f2, h, s) {
  ys = rep(0, 100000)
  if (h!=1) ys[1]=1/2/N
  else ys[1]=0.01
  out=list(T1=NA, T2=NA)
  for (t in 2:length(ys)) {
    ys[t] = ys[t-1] + 2*N*s*ys[t-1]*(1-ys[t-1])*(1-h-(1-2*h)*ys[t-1]) * dt
    if (is.na(out$T1) && ys[t]>f1) out$T1=round(t*dt*2*N)
    if (ys[t]>f2) {out$T2=round(t*dt*2*N-out$T1); break}
  }
  out$y=ys[1:t]
  return(out)
}

calc_avg_path = function(f1, f2, h, s) {
  out_ts = calc_det_path(f1, f2, h, s)
  T2 = (out_ts$T1+out_ts$T2)*2
  out = gen_samples_wf(1, T2, 5, smp_siz, N, s, h, 100)
  paths = matrix(unlist(out$traj_list), ncol = length(out$traj_list[[1]]), byrow = TRUE)
  mean_path = colMeans(paths)
  T1 = which(mean_path>f1)[1]
  T12 = which(mean_path>f2)[1]
  return(list(T1=T1, T2=T12-T1, mean_path=mean_path))
}
#z = simulate_hmm(T1, T2, 10, rep(20, n_ts), N, 0.01, 0.5, 1/2/N)

if (as.numeric(substr(base_dir, 8, 8))>=8) {
  s_list = 0; h_list = 0.5
} else {
  s_list = c(0.0025, 0.005, 0.01, 0.015, 0.02); h_list = c(0, 0.5, 1)
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
          samples_id = as.numeric(substr(base_dir, 8, 8))
          if (samples_id==2 || samples_id==3 || samples_id==4) {
            out_ts = calc_avg_path(f1, f2, h, s)
            out = gen_samples_wf(out_ts$T1, out_ts$T2, n_ts, smp_siz, N, s, h, m, f0)
          } else {
            load(paste(substr(base_dir,1,7), samples_id-3, '/sample',i_params,'.RData',sep=''))
            out = gen_samples_wf(params_data[1,6], params_data[1,5]-params_data[1,6],
                                 n_ts, smp_siz, N, s, h, m, f0)
          }
        } else {
          out = gen_samples_wf(T1, T2, n_ts, smp_siz, N, s, h, m, f0)
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
        matplot(counts[3,,seq(1, dim(params_data)[1],5)], type='l', ylab='',
                main=paste('file_id=', i_params, ', h=', h, ', s=', s, '\n', sep=''))
        points(rowMeans(counts[3,,]), pch=21, col='red', bg='red', cex=1)
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
        cat('h=', h, ', s=', s, ', bias sel_cof =',
            median(res[,1], na.rm=T)/2/N-s_true, '(',
            round((median(res[,1], na.rm=T)/2/N-s_true)/s_true*100), '%), bias age=',
            median(res[,2]*2*N), '(',
            round((median(res[,2])*2*N-t1_true)/t1_true*100), '%), #NA=', sum(is.na(res[,1])), '\n')
        if (s!=0 || h==0.5) {
          boxplot(res[,1]/2/N, main = paste('h=', h, ', s=', s, '\n', sep=''))
          ogn_pt = 1 - 0.4; end_pt = 1 + 0.4
          segments(x0 = ogn_pt, x1 = end_pt, y0 = s_true, col = "red", lty = 2, lwd = 2)
          if (!is.na(res[1,2])) {
            boxplot(res[,2]*2*N, main = paste('h=', h, ', s=', s, '\n', sep=''))
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
