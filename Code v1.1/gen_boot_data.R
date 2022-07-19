i_base_dir = 30 # bootstrap samples for sample0

if (i_base_dir == 30) {
  load('samples0/sample0.RData')
} else if (i_base_dir == 31) {
  load('samples0/sample1.RData')
}

counts0 = counts
params_data0 = params_data
params_data = params_data0[2,,drop=FALSE]
counts1 = counts[,,2]
n1 = sum(counts1[2,])
ber_raw = matrix(0, 2, n1)
j0 = 0
for (i in 1:dim(counts1)[2]) {
  j9 = j0 + counts1[2,i]
  ber_raw[1, (j0+1):j9] = counts1[1,i]
  if (counts1[3,i]>=1) {
    ber_raw[2, (j0+1):(j0+counts1[3,i])] = 1
  }
  #if (j0+1+counts1[3,i]<=n1) ber_raw[2, (j0+1+counts1[3,i]):j9] = 0
  j0 = j9
}

for (i_params in 1:1100) {
  if (i_params %% 100 == 0) cat(i_params, "\n")
  ber_resampled = ber_raw[, sort(sample.int(n1, replace=TRUE))]
  ts = unique(ber_resampled[1,])
  counts = array(0, c(3, length(ts), 1))
  counts[1, ,1] = ts
  for (i in 1:length(ts)) {
    js = which(ber_resampled[1,]==ts[i])
    counts[2,i,1] = length(js)
    counts[3,i,1] = sum(ber_resampled[2,js])
  }
  params_data[1,2] = length(ts)
  file_name = sprintf('samples%d/sample%d.RData', i_base_dir, i_params)
  save(counts, params_data, file=file_name)
}

# generate bootstrap data, resample from the 65 time points
i_base_dir = 40 # bootstrap samples for sample0

if (i_base_dir == 40) {
  load('samples0/sample0.RData')
} else if (i_base_dir == 41) {
  load('samples0/sample1.RData')
}

counts0 = counts
params_data0 = params_data
params_data = params_data0[2,,drop=FALSE]
counts1 = counts[,,2]
n_ts = dim(counts1)[2]

for (i_params in 1:1100) {
  if (i_params %% 100 == 0) cat(i_params, "\n")
  bin_resampled = counts1[, sort(sample.int(n_ts, replace=TRUE))]
  ts = unique(bin_resampled[1,])
  counts = array(0, c(3, length(ts), 1))
  counts[1, ,1] = ts
  for (i in 1:length(ts)) {
    js = which(bin_resampled[1,]==ts[i])
    counts[2,i,1] = sum(bin_resampled[2,js])
    counts[3,i,1] = sum(bin_resampled[3,js])
  }
  params_data[1,2] = length(ts)
  file_name = sprintf('samples%d/sample%d.RData', i_base_dir, i_params)
  save(counts, params_data, file=file_name)
}
