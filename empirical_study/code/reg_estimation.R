library(locpol)
library(mlrv)
library(PLRModels)
library(expm)
library(MASS)
#### functions ####
# tvReg on designed point 
tv_beta = function(y_dat,x_dat,hr,t){
  n = nrow(x_dat)
  p = ncol(x_dat)
  Ts = (1:n)/n
  S_l.tv = function(t,l){
    out=matrix(0,nrow= p,ncol = p)
    for (i in 1:n) {
      out = out+x_dat[i,] %*% t(x_dat[i,]) * ((Ts[i]-t)/hr)^l * Kernel_ep( (Ts[i]-t)/hr )
    }
    out/(n*hr)
  }
  R_l.tv = function(t,l){
    out=matrix(0,nrow=p,ncol=1)
    for (i in 1:n) {
      out = out+x_dat[i,] * y_dat[i] * ((Ts[i]-t)/hr)^l * Kernel_ep( (Ts[i]-t)/hr )
    }
    out/(n*hr)
  }
  
  SS = rbind( cbind(S_l.tv(t,0),S_l.tv(t,1)), cbind(S_l.tv(t,1),S_l.tv(t,2)))
  RR = rbind(R_l.tv(t,0),R_l.tv(t,1))
  
  beta = solve(SS) %*% RR
  beta[1:p]
}


hat_M = function(x_dat,hr,t){
  n = nrow(x_dat)
  p = ncol(x_dat)
  t_star = max(c(hr,min( c(t,1-hr) ) ) )
  Ts = (1:n)/n
  out=matrix(0,nrow= p,ncol = p)
  for (i in 1:n) {
    out = out+x_dat[i,] %*% t(x_dat[i,]) * Kernel_ep( (Ts[i]-t_star)/hr )
  }
  as.matrix(out/(n*hr))
}

GCV.compute = function(y_dat,x_dat,hr){ 
  # Y_hat=X(n,np) * (I_p O_p)(np,2np) * SS^-1(2np,2np) * W(2np,n) * Y(n,1)
  n = nrow(x_dat)
  p = ncol(x_dat)
  Ts = (1:n)/n
  S_l.tv = function(t,l){
    out=matrix(0,nrow= p,ncol = p)
    for (i in 1:n) {
      out = out+x_dat[i,] %*% t(x_dat[i,]) * ((Ts[i]-t)/hr)^l * Kernel_ep( (Ts[i]-t)/hr )
    }
    out/(n*hr)
  }
  X_n_np = matrix(0, nrow = n, ncol = n*p)
  for (i in 1:n) X_n_np[i, ((i-1)*p+1):(i*p) ] = x_dat[i,]
  I_pO_p = matrix(0, nrow = n*p, ncol = 2*n*p)
  for (i in 1:n) I_pO_p[((i-1)*p+1):(i*p), ( 2*(i-1)*p+1 ):( 2*i*p-p ) ] = diag(1,p,p) 
  SS_2np = matrix(0,nrow = 2*n*p, ncol = 2*n*p)
  for (i in 1:n) {
    t = Ts[i]
    SS = rbind( cbind(S_l.tv(t,0),S_l.tv(t,1)), cbind(S_l.tv(t,1),S_l.tv(t,2)))
    SS_2np[ (2*(i-1)*p+1):(2*i*p) , (2*(i-1)*p+1):(2*i*p)] = solve(SS)
  }
  W_2np_n = matrix(0,nrow = 2*n*p,ncol = n)
  for (i in 1:n) {
    W_2np_n[(  2*(i-1)*p+1 ):( (2*i-1)*p ),] = t(x_dat) * Kernel_ep( (Ts-Ts[i])/hr ) / (n*hr)
    W_2np_n[( (2*i-1)*p+1 ):( 2*i*p ),] = t(x_dat) * ((Ts-Ts[i])/hr)  * Kernel_ep( (Ts-Ts[i])/hr ) / (n*hr)
  }
  Q_h = X_n_np %*% I_pO_p %*% SS_2np %*% W_2np_n
  y_hat = Q_h %*% y_dat 
  GCV = mean((y_hat-y_dat)^2) / (1 - sum(diag(Q_h))/n )^2
  return(GCV)
}

GCV.selector = function(y_dat,x_dat){
  candidate = rule_of_thumb(y_dat, x_dat )
  hrs = seq(candidate[1],candidate[2],length.out=10)
  ll.gcv = c()
  for (h in hrs) {
    ll.gcv = c(ll.gcv, GCV.compute(y_dat = y_dat,x_dat = x_dat, hr= h))
  }
  hr.select = hrs[which.min(ll.gcv)]
  c(candidate[1],hr.select,candidate[2])
}


process = function(ts){
  # missing check
  ts = ts[!is.na(ts$yyyy),]
  years = seq(min(ts$yyyy), min( c(max(ts$yyyy),2024) ),by=1 )
  df = matrix(NA,nrow = 12*length(years),ncol = 3)
  for (i in 1:length(years)){
    df[((i-1)*12+1):(i*12) , 1] = years[i]
    df[ ((i-1)*12+1):(i*12) , 2] = seq(1,12,1)
  }
  for (i in 1:nrow(df)){
    if(length(ts[ (ts$yyyy==df[i,1] & ts$mm==df[i,2]) ,3])==1){
      df[i,3] = ts[ (ts$yyyy==df[i,1] & ts$mm==df[i,2]) ,3]
    }
  } 
  
  # interpolation
  df_inter = data.frame(year=df[,1],month=df[,2],tmax=df[,3])
  for (m in which(is.na(df[,3])) ) {
    month_m = df[m,2]
    year_m = df[m,1]
    ts_m = ts(df[df[,2]==month_m,3],start=years[1],end=years[length(years)])
    df_inter[m,3] = approx(x=time(ts_m),y=ts_m,xout = year_m)$y
  }
  
  # eliminate seasonal effect
  # 若第一年不满，去除
  while ( sum(!is.na(df_inter[df_inter$year==years[1],3]))<12 ) {
    df_inter = df_inter[df_inter$year>years[1],]
    years = years[-1]
  }
  # 若最后一年不满，去除
  while ( sum(!is.na(df_inter[df_inter$year==years[length(years)],3]))<12 ) {
    df_inter = df_inter[df_inter$year<years[length(years)],]
    years = years[-length(years)]
  }
  
  ts_data = ts(df_inter$tmax,start=c(years[1],1),end=c(years[length(years)],12),frequency=12)
  
  stl_data = stl(ts_data,s.window = "periodic")
  seasonal_adj_data = ts_data - stl_data$time.series[,1]
  return( list(ts_data=ts_data,stl_data=stl_data$time.series,seasonal_adj=seasonal_adj_data) )
}

# 临时函数，求反函数放进去 m_I^-1(t) 基于特殊的Kernel, 可以解析地算出
sum_jack <- function(lower){
  if(lower >= 1){
    return(0)
  }else if(lower < -1){
    return(1)
  }else if(lower >= -1 & lower < 1){
    out = 0.5 - 0.75 * lower + 0.25 * lower^3
    return(out)
  }
}

# local linear estimator with jack-knife
locpoly_jack <- function(x, y, xeval, hr){
  XY_data <- as.data.frame(cbind(x,y))
  m_hat <- locpol(y~x,data=XY_data,deg = 1,bw = hr,kernel = EpaK,xeval = xeval )
  m_hat <- m_hat$lpFit
  m_hat_2 <- locpol(y~x,data=XY_data,deg = 1,bw = hr/sqrt(2),kernel = EpaK, xeval = xeval )
  m_hat_2 <- m_hat_2$lpFit
  m_tilde <- 2*m_hat_2[,2]-m_hat[,2]
  dm_est <-  2*m_hat_2[,3]-m_hat[,3]
  out <- list(m=m_tilde,dm=dm_est)
  out
}

# inverse monotone estimator
## Ts is vector, lambda=0 without penalization
inverse_jack <- function(Ts, ll, N, hd, hr){ # Ts support set, ll local linear estimation
  out = c()
  for(t in Ts){
    temp <- sapply((ll-t)/hd, sum_jack)
    out = c(out,mean(temp))
  }
  out
}

# kernel
Kernel_ep <- function(u){
  support <- as.numeric(abs(u)<=1)
  0.75*(1-u^2)*support
}
## kernel with jack-knife
Kernel_tilde <- function(u){
  2*sqrt(2)*Kernel_ep(sqrt(2)*u)-Kernel_ep(u)
}
## kernel with jack-knife and boundary
## output K(k/n,i/N)


S_l <- function(Ts,l,hr){ # 
  u1 = (1-Ts)/hr; l1 = -Ts/hr
  u1[u1 >1] = 1
  l1[l1 < -1] = -1
  out = 0.75*(u1^(l+1)/(l+1)-u1^(l+3)/(l+3)) - 0.75*(l1^(l+1)/(l+1)-l1^(l+3)/(l+3))
}

## kernel start involved in SCB
K_star <- function(n,N,hr){ # output a n*N matrix indicating K_hr(i/N-j/n)
  # k from 1 to n
  tk = (1:n)/n
  # i from 1 to N
  ti = (1:N)/N
  S_0 = S_l(ti,l=0,hr=hr)
  S_1 = S_l(ti,l=1,hr=hr)
  S_2 = S_l(ti,l=2,hr=hr)
  SS = S_0*S_2-S_1^2
  out = matrix(NA,nrow = n,ncol = N)
  for (k in 1:n) {
    us = (tk[k]-ti)/hr
    K = Kernel_ep(us)
    out[k,] = (S_2*K-S_1*K*us)/SS
  }
  out/(n*hr)
}

K_star_tilde <- function(n,N,hr){
  out = 2*K_star(n,N,hr=hr/sqrt(2))-K_star(n,N,hr=hr)
  out
}

hat_M = function(x_dat,hr,t){
  n = nrow(x_dat)
  p = ncol(x_dat)
  t_star = max(c(hr,min( c(t,1-hr) ) ) )
  Ts = (1:n)/n
  out=matrix(0,nrow= p,ncol = p)
  for (i in 1:n) {
    out = out + x_dat[i,] %*% t(x_dat[i,]) * Kernel_ep( (Ts[i]-t_star)/hr )
  }
  out/(n*hr)# symmetric
}

# estimate long-run covariance matrix
long_run_var_est_high = function(residuals,t.out,m,tau){ # residuals n*p matrix
  residuals = as.matrix(residuals)
  n = nrow(residuals)
  X = (1:n)/n
  
  delta = list()
  for (i in 1:n) {
    Qs = colSums(residuals[max(1,i-m):min(n,i+m),]) 
    delta[[i]] = Qs %*% t(Qs)/(2*m+1)
  } 
  
  # NW estimation
  W = matrix(NA,nrow = n,ncol = length(t.out))
  Y.out = list()
  for (t in 1:length(t.out)) W[,t] = Kernel_ep( (t.out[t]-X)/tau ) 
  for (t in 1:length(t.out)) {
    Y.out[[t]] = Reduce('+', Map('*',W[,t],delta))/sum(W[,t])
    Y.out[[t]] = sqrtm(Y.out[[t]])
  }
  return(Y.out)
}


Weights_star = function(ll,support,hd){
  W = matrix(NA,nrow = N,ncol = length(support))
  W_star = W
  for (t in 1:length(support)) W[,t] = Kernel_ep( (ll[2:(N+1)]-support[t])/hd )
  for (t in 1:length(support)) W_star[,t] = W[,t]/sum(W[,t])
  return(W_star)
}

#### preprocess sunshine and temperature data ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
files = dir('../raw_data')
# not include locations.csv
files = files[files!='locations.csv']
ts_all.tmax = list()
ts_all.sun = list()
# prerocess data in heathrow and lerwick (i = 12 and 14)
names = c()
for (i in c(12,14)) {
  name = substr(files[i],start = 1,stop = nchar(files[i])-8)
  names = c(names,name)
  raw_data = read.table(paste0('../raw_data/',files[i]),skip = 7,col.names=c('yyyy', 'mm', 'tmax','tmin','af','rain','sun')
                        ,fill = T)
  raw_data = apply(raw_data, 2, as.numeric)
  
  ts.tmax = raw_data[,1:3]
  ts.tmax = as.data.frame(ts.tmax)
  
  ts.sun = raw_data[,c(1,2,7)]
  ts.sun = as.data.frame(ts.sun)
  
  stl_ts.tmax = process(ts.tmax)
  stl_ts.sun = process(ts.sun)

  ts_all.tmax[[name]] = stl_ts.tmax$seasonal_adj
  ts_all.sun[[name]] = stl_ts.sun$seasonal_adj
}


# choose time series ranged from 1979 to 2023 (45 years)
ts_all = list(tmax=ts_all.tmax,sun=ts_all.sun)
ts_all.trunc = ts_all

names_all = ts_all

for (i in 1:2) {
  starts = c()
  ends = c()
  for (j in names) {
    starts = c(starts,start(ts_all[[i]][[j]])[1] )
    ends = c(ends,end(ts_all[[i]][[j]])[1] )
  }
  ts_all.trunc[[i]] = ts_all[[i]][[1]][time(ts_all[[i]][[1]]) < 2024 & time(ts_all[[i]][[1]]) >= 1979]
  for (j in which(ends>=2023 & starts<=1979) ) {
    ts_all.trunc[[i]] = cbind(ts_all.trunc[[i]],  ts_all[[i]][[j]][time(ts_all[[i]][[j]]) < 2024 & time(ts_all[[i]][[j]]) >= 1979])
  }
  ts_all.trunc[[i]] = data.frame(ts_all.trunc[[i]][,-1])
  colnames(ts_all.trunc[[i]]) = names
}



#### estimation ####

# about 14 min
sun_reg_estimations = list()
for (i in 1:2) {
  n = length(ts_all.trunc$tmax[,i])
  X = (1:n)/n
  N = 4000
  xevals = seq(0,1,by=1/N)
  tmax = ts_all.trunc$tmax[,i]
  sun = ts_all.trunc$sun[,i]
  XX = cbind(rep(1,n),sun)
  hr = GCV.selector(tmax,XX)[2]
  hd = hr^3
  fit1 = sapply(xevals, FUN = tv_beta,y_dat = tmax,x_dat = XX, hr=hr/sqrt(2))
  fit2 = sapply(xevals, FUN = tv_beta,y_dat = tmax,x_dat = XX, hr=hr)
  ll = t(2*fit1)-t(fit2)
  y.hat = ll[(1:n)*floor(N/n),1]+ll[(1:n)*floor(N/n),2] * sun
  
  residuals = y.hat-tmax
  residuals_m = XX*residuals
  hat_Sigma = long_run_var_est_high(residuals = residuals_m,t.out = X,m = floor(n^(2/7)), tau = n^(-1/7))
  hat_M_station = list()
  for (j in 1:n) {
    hat_M_station[[j]] = hat_M(x_dat = XX,hr=hr,t=X[j])
  }
  Sigma = rep(NA , n)
  for (j in 1:n) {
    Sigma_temp = solve(hat_M_station[[j]])  %*% ( (hat_Sigma[[j]])^2 ) %*% solve(hat_M_station[[j]])
    Sigma[j] = t(c(0,1)) %*% Sigma_temp %*% c(0,1)
    Sigma[j] = sqrt(Sigma[j])
  }
  # monotone estimation
  support = seq( min(ll[-1,2]), max(ll[-1,2]), length.out=n)
  monotone = inverse_jack(Ts=support,ll=ll[,2],N=N,hd=hd,hr=hr)
  weights = Weights_star(ll=ll[,2],support = support,hd=hd)
  sun_reg_estimations[[i]] = list(monotone=monotone,support=support,weight=weights,ll=ll,hr.gcv=hr,data=cbind(tmax,XX),Sigma=Sigma)
}

save(sun_reg_estimations,file= '../estimation/regression_sun.RData')



