process_data = function(counts, params_data, i_sample) {
  n_ts = params_data[i_sample,2]
  if (length(dim(counts))==2) {
    data4 = counts[, n_ts:1]
  } else {
    data4 = counts[, n_ts:1, i_sample]
  }
  data4[1,] = -data4[1,] + data4[1,1]
  if (data4[1,n_ts]<3) {
    data4[1,] = round(data4[1,], digits=4)
  } else {
    data4[1,] = round(data4[1,]/2/params_data[i_sample,1], digits=4)
  }
  return(data4)
}

gen_grid = function(grid_option) {
  if (grid_option==1) {return(seq(0, 1, length.out=1001))}
  if (grid_option==2) {
    #x2 = 0.05; n1=70; n2 = 15; xs_piece1 = exp(seq(-6, log(x2), length.out = n2))
    x2 = 0.03; n1=120; n2 = 30; xs_piece1 = exp(seq(-6, 0, length.out = n2) + log(x2))
    xs = c(0, xs_piece1[1:(n2-1)], seq(x2,1-x2,length.out=n1), 1-rev(xs_piece1[1:(n2-1)]), 1)
    return(xs)
  }
  if (grid_option==3) {
    x2 = 0.03; n1=120; n2 = 30; xs_piece1 = exp(seq(-6.5, 0, length.out = n2) + log(x2))
    xs = c(0, xs_piece1[1:(n2-1)], seq(x2,1-x2,length.out=n1), 1-rev(xs_piece1[1:(n2-1)]), 1)
    return(xs)
  }
  x2 = 0.5; n1=80; #xs_piece1 = exp(seq(-8, log(x2), length.out = n1))
  xs_piece1 = exp(-seq(-3, 0, length.out=n1+1)^2+log(0.5)) # stored in resb..
  return( c(0, xs_piece1[1:n1], 0.5, 1-rev(xs_piece1[1:n1]), 1) )
  x3 = 0.02; n1=30
  #xs_piece1 = exp(-seq(-3, -0.5, length.out=n1+1)^2)
  xs_piece1 = exp(seq(-8, -0, length.out=n1+1))
  xs_piece2 = xs_piece1[1:n1]/xs_piece1[n1+1]*x3
  xs = c(0, xs_piece2, seq(x3, 1-x3, length.out=38), 1-rev(xs_piece2), 1)
  return(xs)
}

gen_diff_coeff = function(xs) {
  nx = length(xs)
  dxs = diff(xs)
  #dxs2=(xs[3:nx]-xs[1:(nx-2)])/2
  b1 = (dxs[1:(nx-2)] / dxs[2:(nx-1)])^2; a1 = -b1 + 1
  abc1 = dxs[1:(nx-2)] * (dxs[1:(nx-2)] / dxs[2:(nx-1)]+1)
  b0 = ((dxs[2]+dxs[1]) / dxs[1])^2; a0 = -b0 + 1
  abc0 = (dxs[2]+dxs[1]) * ((dxs[2]+dxs[1])/dxs[1]-1)
  b2 = dxs[1:(nx-2)] / dxs[2:(nx-1)]; a2 = -b2 - 1
  abc2 = dxs[1:(nx-2)] * (dxs[2:(nx-1)]+dxs[1:(nx-2)]) /2
  b3 = - (dxs[2]+dxs[1]) / dxs[1]; a3 = -b3 - 1
  abc3 = (dxs[2]+dxs[1]) * dxs[2] / 2
  return(list('a0'=a0, 'b0'=b0, 'abc0'=abc0, 'a1'=a1, 'b1'=b1, 'abc1'=abc1, 
              'a2'=a2, 'b2'=b2, 'abc2'=abc2, 'a3'=a3, 'b3'=b3, 'abc3'=abc3))
}

diff1 = function(temp, dxsc, nx) {
  d1 = rep(0, nx)
  if (length(dxsc)==1 && is.na(dxsc)) {
    d1[2:(nx-1)] = (temp[3:nx]-temp[1:(nx-2)]) * nx / 2
    d1[1] = (-3*temp[1] + 4*temp[2] - temp[3]) * nx / 2
    d1[nx] = (-3*temp[nx] + 4*temp[nx-1] - temp[nx-2]) * nx / 2
  } else {
    d1[2:(nx-1)] = (dxsc$a1*temp[2:(nx-1)]+dxsc$b1*temp[3:nx]-temp[1:(nx-2)]) / dxsc$abc1
    d1[1] = (dxsc$a0*temp[1]+dxsc$b0*temp[2]-temp[3]) / dxsc$abc0
    d1[nx] = -(dxsc$a0*temp[nx]+dxsc$b0*temp[nx-1]-temp[nx-2]) / dxsc$abc0
  }
  return(d1)
}

diff2 = function(temp, dxsc, nx) {
  d2 = rep(0, nx)
  if (length(dxsc)==1 && is.na(dxsc)) {
    nx2 = nx * nx
    d2[2:(nx-1)] = (-2*temp[2:(nx-1)]+temp[3:nx]+temp[1:(nx-2)]) * nx2
    d2[1] = (temp[1]-2*temp[2]+temp[3]) * nx2
    d2[nx] = (temp[nx]-2*temp[nx-1]+temp[nx-2]) * nx2
  } else {
    d2[2:(nx-1)] = (dxsc$a2*temp[2:(nx-1)]+dxsc$b2*temp[3:nx]+temp[1:(nx-2)]) / dxsc$abc2
    d2[1] = (dxsc$a3*temp[1]+dxsc$b3*temp[2]+temp[3]) / dxsc$abc3
    d2[nx] = (dxsc$a3*temp[nx]+dxsc$b3*temp[nx-1]+temp[nx-2]) / dxsc$abc3
  }
  return(d2)
}

calc_lkhd_outer = function(xs, h, alpha, data4, mode_age=0, mode_absorb='boundary', cond=TRUE,
                           cutoff_factor=5, nt_factor=1, nt_factor_max=100) {
  nt_factor1 = nt_factor
  repeat {
    cat('nt_factor=', nt_factor1, '\n')
    out = calc_lkhd(xs, h, alpha, data4, mode_age, mode_absorb, cond, cutoff_factor, nt_factor1)
    if (nt_factor1>nt_factor_max || !is.na(out$lkhd1)) {break}
    nt_factor1 = nt_factor1 * 2
  }
  if (is.na(out$lkhd1)) return(list(lkhd1=0, t0_max1=0))
  return(out)
}

#load('samples2/sample233.RData')
#xs=gen_grid(1); h=1; alpha=-640; data4=process_data(counts, params_data, 83)
#mode_absorb='boundary'; mode_age=1; cond=TRUE; cutoff_factor=5; nt_factor=64
# mode_age: 1 'age' or 0 'no_age'
calc_lkhd = function(xs, h, alpha, data4, mode_age=0,
                     mode_absorb='boundary', cond=TRUE, cutoff_factor=5, nt_factor=4) {
  nx = length(xs)
  if (max(diff(xs))-min(diff(xs))<1e-7) {
    dxsc = NA; dxsc_int = NA
  } else {
    dxsc = gen_diff_coeff(xs); dxsc_int = gen_diff_coeff(xs[2:nx])
  }
  dxs2=c((xs[3:nx]-xs[1:(nx-2)])/2, xs[nx]-xs[nx-1])
  dt = 0.0001 ###################### check this is the same as in sim_sel
  nt = max(1, round(abs(alpha)/100))*nt_factor # *4 for fine grid
  t8 = ceiling(data4[1,max(which(data4[3,]>0))]/dt)
  if (mode_age!=0) {T2 = t8 * 10} else {T2 = ceiling(data4[1,ncol(data4)]/dt)}
  mass_absorbed = rep(0, T2); max_absorbed = 0
  bx = alpha*xs*(1-xs)*(1-h-(1-2*h)*xs); sigmax2 = xs*(1-xs)

  C = matrix(0, nrow=nx, ncol=T2)
  dlogC = matrix(0, nrow=nx-1, ncol=T2)
  C[,1] = 1; C[1,1]=0; #C[nx,1]=0
  j_ts = 2; flag_smp_use = rep(0, ncol(data4)); flag_smp_use[1] = 1
  
  drift1 = rep(0, nx)
  p2 = matrix(0, nrow=nx, ncol=T2)
  p2[,1] = dbinom(data4[3,1], data4[2,1], xs)
  for (t in 2:T2) {
    if (cond==TRUE) {
      temp = C[,t-1]
      for (i in 1:nt) {
        dC = (bx*diff1(temp, dxsc, nx) + sigmax2/2*diff2(temp, dxsc, nx))*dt/nt
        temp = pmin(1, pmax(0, temp + dC))
      }
      C[,t] = temp
      logC = log(C[2:nx,t-1]); dlogC[,t] = diff1(logC, dxsc_int, nx-1)
      dlogC[!is.finite(dlogC[,t]),t] = 0
    }
    temp = p2[,t-1]
    drift1[2:nx] = dlogC[,t]*sigmax2[2:nx]
    for (i in 1:nt) {
      d1 = diff1(temp, dxsc, nx); d2 = diff2(temp, dxsc, nx)
      temp = temp + ((bx+drift1)*d1 + sigmax2/2*d2)*dt/nt
    }
    p2[,t] = pmax(0, temp)
    dp2 = diff(p2[,t])
    i_flips = which(dp2[2:(nx-1)]*dp2[1:(nx-2)]<0)
    # check the shape of C and p2 are still good, otherwise return NA result for finer t grid
    if (min(diff(C[,t]))< -max(abs(C[,t]))/10 || max(p2[,t])>2 ||
        sum(c(abs(dp2[i_flips]), abs(dp2[i_flips+1]))>max(abs(dp2))/100)>20) {
      return(list('lkhd1'=NA, 't0_max1'=0, 't_last'=t))
    }
    if (t>t8) {
      sigmax2sqrt = sqrt(sigmax2*dt)
      if (mode_absorb=='normal') {
        mass_absorbed[t] = sum(p2[2:nx,t]*sapply(2:nx,
          function(ix) pnorm(0, mean=xs[ix]-(bx[ix]+drift1[ix])*dt, sd=sigmax2sqrt[ix])))
      } else {
        mass_absorbed[t] = max(0, d1[1]) * dt
      }
      max_absorbed = max(max_absorbed, mass_absorbed[t])
      if (is.na(max_absorbed)) {
        return(list('lkhd1'=NA))
      }
      if (t>t8+10 && (mass_absorbed[t]<max_absorbed/cutoff_factor)) break #is.na(max_absorbed)
    }
    if (j_ts<=ncol(data4) && flag_smp_use[j_ts]==0 && t>=round(data4[1,j_ts]/dt)) {# && t<=t8) {
      p2[,t] = p2[,t] * dbinom(data4[3,j_ts], data4[2,j_ts], xs)
      flag_smp_use[j_ts] = 1; j_ts = j_ts + 1
    }
    p2[1,t] = 0
    #if (t%%100==0) cat('t=', t, max(p2[,t]), '\n')
  }
  if (mode_age==0) {
    return(list('t0_max1'=NA, 'lkhd1'=sum(dxs2*p2[2:nx,t8]),
                'lkhd_t1'=NA, 'C'=C, 'p2'=p2, 't8'=t8))
  }

  if (FALSE) {
    cdf = c(0,cumsum(dxs2*p2[2:nx,t8])); max_cdf = max(cdf)
    n_interp = 100
    points_p = seq(1/n_interp/2, 1, 1/n_interp)
    points_x = matrix(NA, n_interp, T2)
    weight_per_point = rep(NA, T2)
    new_points = matrix(NA, nrow=n_interp, ncol=n_interp)
    points_x[,t8] = approx(cdf, xs, seq(max_cdf/n_interp/2, max_cdf, max_cdf/n_interp))$y
    weight_per_point[t8] = max_cdf
    for (t in (t8+1):T2) {
      sigmax3 = points_x[,t-1]*(1-points_x[,t-1])
      drift3 = alpha*sigmax3*(1-h-(1-2*h)*points_x[,t-1])*2
      if (cond==TRUE) {
        drift4 = approx(xs[2:nx], dlogC[,t]*sigmax2[2:nx], points_x[,t-1], rule=2)$y
      } else {
        drift4 = rep(0, length(drift3))
      }
      sigmax4 = sqrt(sigmax3*dt)
      for (ix in 1:n_interp) {
        new_points[ix,] = qnorm(points_p, mean=points_x[ix,t-1]-(drift4[ix]+drift3[ix])*dt,
                                sd=sigmax4[ix])
      }
      mass_absorbed[t] = weight_per_point[t-1]*sum(new_points<0)/n_interp/n_interp
      max_absorbed = max(max_absorbed, mass_absorbed[t])
      if (t>t8+10 && mass_absorbed[t]<max_absorbed/5) break
      if (FALSE && j_ts<=ncol(data4) && flag_smp_use[j_ts]==0 && t>=round(data4[1,j_ts]/dt)) {
        new_points1 = new_points[new_points>=0 & new_points<=1]
        weights1 = sapply(new_points1, function(x)
          sum(dbinom(seq(data4[3,j_ts],data4[3,j_ts]+data4[4,j_ts]), data4[2,j_ts], x)))
        new_points2 = wtd.quantile(new_points1, weights=weights1, probs=points_p)
        attr(new_points2, 'names') = NULL
        points_x[,t] = new_points2
        weight_per_point[t] = weight_per_point[t-1]*sum(weights1)/n_interp/n_interp
        flag_smp_use[j_ts] = 1; j_ts = j_ts + 1
      } else {
        new_points1 = new_points[new_points>=0 & new_points<=1]
        points_x[,t] = quantile(new_points1, probs=points_p)
        weight_per_point[t] = weight_per_point[t-1]*length(new_points1)/n_interp/n_interp
      }
    }
  }
  if (is.na(max(mass_absorbed)) || max(mass_absorbed)==0) {
    return(list('t0_max1'=NA, 'lkhd1'=0, 'lkhd1_t'=mass_absorbed, 'C'=C, 'p2'=p2, 't8'=t8,
                'dlogC'=dlogC))
  }
  T0_max = which.max(mass_absorbed)
  #if (interactive()) cat(T0_max, mass_absorbed[T0_max], '\n')
  return(list('t0_max1'=T0_max*dt, 'lkhd1'=mass_absorbed[T0_max], 'lkhd1_t'=mass_absorbed,
              'C'=C, 'p2'=p2, 't8'=t8, 'dlogC'=dlogC, 'sigmax2'=sigmax2))
}

