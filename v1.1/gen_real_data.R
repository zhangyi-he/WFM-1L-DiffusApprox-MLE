load('1locus_without_missingtimepoint.RData')
data1 = rbind(-asip_unmiss[,'times']/8, asip_unmiss[,'tothap'], asip_unmiss[,'mutant_cnt'])
data1 = data1[,3:dim(data1)[2]]
last_smp_ybp1 = -data1[1,ncol(data1)]
data1[1,] = round(data1[1,] + last_smp_ybp1)
data2 = rbind(-mc1r_unmiss[,'times']/8, mc1r_unmiss[,'tothap'], mc1r_unmiss[,'mutant_cnt'])
data2 = data2[,3:dim(data2)[2]]
last_smp_ybp2 = -data2[1,ncol(data2)]
data2[1,] = round(data2[1,] + last_smp_ybp2)
data1_all = array(NA, dim=c(3,dim(data1)[2],3))
data2_all = array(NA, dim=c(3,dim(data2)[2],3))
for (i in 1:3) {
  data1_all[,,i] = data1; data2_all[,,i] = data2
}

pop_sizes = c(8000, 16000, 32000)
#counts = list(ts=rbind(data1[1,], data1[1,]), smp_siz=rbind(data1[2,], data1[2,]), counts=rbind(data1[3,], data1[3,]))
#params_data = rbind(c(pop_siz,dim(data1)[2],1,NA,NA), 
#                    c(pop_siz,dim(data1)[2],1,NA,NA))
counts = data1_all; params_data = matrix(NA, nrow=3, ncol=7)
for (i in 1:3) {params_data[i,] = c(pop_sizes[i],dim(data1)[2],1,NA,NA,NA,last_smp_ybp1)}
save(counts, params_data, file='samples0/sample0.RData')
#counts = list(ts=rbind(data2[1,], data2[1,]), smp_siz=rbind(data2[2,], data2[2,]), counts=rbind(data2[3,], data2[3,]))
#params_data = rbind(c(pop_siz,dim(data2)[2],1,NA,NA), 
#                    c(pop_siz,dim(data2)[2],1,NA,NA))
counts = data2_all; params_data = matrix(NA, nrow=3, ncol=7)
for (i in 1:3) {params_data[i,] = c(pop_sizes[i],dim(data2)[2],1,NA,NA,NA,last_smp_ybp2)}
save(counts, params_data, file='samples0/sample1.RData')

################# generate bootstrap samples ################# 
m=600

# i_params<100: sample times specified
i_params = 31
if (i_params==30) {
  load('samples0/sample0.RData')
} else if (i_params==31) {
  load('samples0/sample1.RData')
}
n_ts=ncol(counts[,,2]);
params_data0 = params_data[2,]
counts0 = counts[,,2]
counts = array(NA, dim=c(3,n_ts,m))
params_data = matrix(NA, nrow=m, ncol=7)
for (i in 1:m) {
  counts[1,,i] = counts0[1,]
  counts[2,,i] = counts0[2,]
  counts[3,,i] = mapply(function(smsz1, p) 
    min(smsz1, max(0, round(rbinom(1, smsz1, p) + rnorm(1, mean=0, sd=sqrt(smsz1)/5)))),
    counts0[2,], counts0[3,]/counts0[2,]) 
  params_data[i,] = params_data0
}
save(counts, params_data, file=paste('samples0/sample',i_params,'.RData',sep=''))

#sample0: [,1]   [,2] [,3]         [,4]
#[1,]  28.80658 0.2920406    1 2.310579e-40
#[2,]  56.79012 0.2293203    1 7.556035e-41
#[3,] 113.58025 0.1678602    1 3.217199e-41
#sample1:> [,1]   [,2] [,3]         [,4]
#[1,] 200.0000 0.11694063    1 6.138891e-37
#[2,] 413.1687 0.08122031    1 2.084360e-37
#[3,] 867.4897 0.05976016    1 1.642861e-37


################# convert time units of randomly grouped samples ################# 
i_params = 1
m = 200
load(paste('samples0/sample', i_params, '.RData', sep=''))
counts_old = counts[,,2]
params_data_old = params_data[2,]
n_ts_old = ncol(counts_old)
n_ts_new = 9
counts = array(NA, dim=c(3, n_ts_new, m))
params_data = matrix(NA, nrow=m, ncol=7)
i = 1
while (i<=m) {
  t_cuts = c(-Inf, sort(runif(n_ts_new - 1, counts_old[1,1], counts_old[1,n_ts_old])), Inf)
  for (k in 1:n_ts_new) {
    js = which((counts_old[1,]>=t_cuts[k]) & (counts_old[1,]<t_cuts[k+1]))
    if (length(js)==0) break
    counts[1,k,i] = mean(counts_old[1,js])
    counts[2,k,i] = sum(counts_old[2,js])
    counts[3,k,i] = sum(counts_old[3,js])
  }
  if (!is.na(counts[1,n_ts_new,i])) {
    params_data[i,c(1,3)] = params_data_old[c(1,3)]
    params_data[i,2] = n_ts_new
    params_data[i,7] = params_data_old[7] - counts[1,n_ts_new,i]
    i=i+1
  }
}
save(counts, params_data, file=paste('samples0/sample', i_params+10, '.RData', sep=''))

######## randomly grouped samples, same as above, but with 1100 replicates ######## 
i_params = 0
m = 1100
load(paste('samples0/sample', i_params, '.RData', sep=''))
counts_old = counts[,,2]
params_data_old = params_data[2,]
n_ts_old = ncol(counts_old)
n_ts_new = 9
counts = array(NA, dim=c(3, n_ts_new, m))
params_data = matrix(NA, nrow=m, ncol=7)
i = 1
while (i<=m) {
  t_cuts = c(-Inf, sort(runif(n_ts_new - 1, counts_old[1,1], counts_old[1,n_ts_old])), Inf)
  for (k in 1:n_ts_new) {
    js = which((counts_old[1,]>=t_cuts[k]) & (counts_old[1,]<t_cuts[k+1]))
    if (length(js)==0) break
    counts[1,k,i] = mean(counts_old[1,js])
    counts[2,k,i] = sum(counts_old[2,js])
    counts[3,k,i] = sum(counts_old[3,js])
  }
  if (!is.na(counts[1,n_ts_new,i])) {
    params_data[i,c(1,3)] = params_data_old[c(1,3)]
    params_data[i,2] = n_ts_new
    params_data[i,7] = params_data_old[7] - counts[1,n_ts_new,i]
    i=i+1
  }
}
save(counts, params_data, file=paste('samples0/sample', i_params+40, '.RData', sep=''))

