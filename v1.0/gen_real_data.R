load('1locus_without_missingtimepoint.RData')
data1 = rbind(-asip_unmiss[,'times']/8, asip_unmiss[,'tothap'], asip_unmiss[,'mutant_cnt'])
data1 = data1[,3:dim(data1)[2]]
data1[1,] = round(data1[1,] - data1[1,ncol(data1)])
data2 = rbind(-mc1r_unmiss[,'times']/8, mc1r_unmiss[,'tothap'], mc1r_unmiss[,'mutant_cnt'])
data2 = data2[,3:dim(data2)[2]]
data2[1,] = round(data2[1,] - data2[1,ncol(data2)])
data1_all = array(NA, dim=c(3,dim(data1)[2],3))
data2_all = array(NA, dim=c(3,dim(data2)[2],3))
for (i in 1:3) {
  data1_all[,,i] = data1; data2_all[,,i] = data2
}

pop_sizes = c(8000, 16000, 32000)
#counts = list(ts=rbind(data1[1,], data1[1,]), smp_siz=rbind(data1[2,], data1[2,]), counts=rbind(data1[3,], data1[3,]))
#params_data = rbind(c(pop_siz,dim(data1)[2],1,NA,NA), 
#                    c(pop_siz,dim(data1)[2],1,NA,NA))
counts = data1_all; params_data = matrix(NA, nrow=3, ncol=6)
for (i in 1:3) {params_data[i,1:3] = c(pop_sizes[i],dim(data1)[2],1)}
save(counts, params_data, file='sample0.RData')
#counts = list(ts=rbind(data2[1,], data2[1,]), smp_siz=rbind(data2[2,], data2[2,]), counts=rbind(data2[3,], data2[3,]))
#params_data = rbind(c(pop_siz,dim(data2)[2],1,NA,NA), 
#                    c(pop_siz,dim(data2)[2],1,NA,NA))
counts = data2_all; params_data = matrix(NA, nrow=3, ncol=6)
for (i in 1:3) {params_data[i,1:3] = c(pop_sizes[i],dim(data2)[2],1)}
save(counts, params_data, file='sample1.RData')

################# generate bootstrap samples ################# 
m=200; N=32000

# i_params<100: sample times specified
i_params = 21
if (i_params==20) {
  load('sample0.RData')
} else if (i_params==21) {
  load('sample1.RData')
}
n_ts=ncol(counts[,,2]);
params_data0 = params_data[2,]
counts0 = counts[,,2]
counts = array(NA, dim=c(3,n_ts,m))
params_data = matrix(NA, nrow=m, ncol=6)
for (i in 1:m) {
  counts[1,,i] = counts0[1,]
  counts[2,,i] = counts0[2,]
  counts[3,,i] = mapply(function(smsz1, p) rbinom(1, smsz1, p), counts0[2,], counts0[3,]/counts0[2,])
  params_data[i,] = params_data0
}
save(counts, params_data, file=paste('sample',i_params,'.RData',sep=''))

#sample0: out> [1,]   28 0.225    1 9.324602e-36   28 0.225    1 9.324602e-36
#sample1:> out [1,]  216 0.0789    1 3.33806e-32  216 0.0789    1 3.33806e-32

################# convert time units of randomly grouped samples ################# 
load('samples0/sample10_old.RData')
params_data = cbind(params_data, rep(NA, nrow(params_data)))
counts[1,,] = round(counts[1,,] * 32000)
save(counts, params_data, file=paste('sample10.RData',sep=''))

load('samples0/sample11_old.RData')
params_data = cbind(params_data, rep(NA, nrow(params_data)))
counts[1,,] = round(counts[1,,] * 32000)
save(counts, params_data, file=paste('sample11.RData',sep=''))
