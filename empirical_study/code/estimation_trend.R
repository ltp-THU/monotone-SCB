library(locpol)
library(PLRModels)
library(rstudioapi)
library(mlrv)
#### functions ####
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
Weights_star = function(ll,support,hd){
  W = matrix(NA,nrow = N,ncol = length(support))
  W_star = W
  for (t in 1:length(support)) W[,t] = Kernel_ep( (ll[2:(N+1)]-support[t])/hd )
  for (t in 1:length(support)) W_star[,t] = W[,t]/sum(W[,t])
  return(W_star)
}

GCV.selector = function(Ys){
  p = ncol(Ys)
  n = nrow(Ys)
  X = (1:n)/n
  hrs.opt = rep(NA,p)
  candidate = matrix(NA,p,2)
  for (j in 1:p) {
    candidate[j,] = rule_of_thumb(Ys[,j],as.matrix(X) )
    hrs.opt[j] = np.gcv(cbind(Ys[,j],X),
                        h.seq = seq(candidate[j,1]*n^(1/4-1/6),candidate[j,2]*n^(1/6-1/7)  ,length.out=10), estimator = 'LLP')$h.opt
  }
  return(list(lower = candidate[,1],opt=hrs.opt,upper=candidate[,2]))
}


#### high-dim estimation ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ts_all.trunc = read.csv('../estimation/data_processed.csv')
ts_all.trunc = ts_all.trunc[,-1]
# estimation
monotones = ts_all.trunc
supports = ts_all.trunc
residuals = ts_all.trunc
n = nrow(ts_all.trunc)
X = (1:n)/n
lls = ts_all.trunc
Weights_ls = list()
hrs.opt = GCV.selector(ts_all.trunc)$opt
hds = c()
for (j in 1:ncol(ts_all.trunc)) {
  Y = ts_all.trunc[,j]
  N = 4000
  hr = hrs.opt[j]
  hd = hr^2
  hds = c(hds,hd)
  xevals = seq(0,1,by=1/N)
  XY_data <- as.data.frame(cbind(X,Y))
  colnames(XY_data) = c('X','Y')
  fit1 = locpol(Y~X,data=XY_data,deg = 1,bw = hr/sqrt(2),kernel = EpaK, xeval = xevals )
  fit2 = locpol(Y~X,data=XY_data,deg = 1,bw = hr,kernel = EpaK,xeval = xevals )
  ll = 2*fit1$lpFit[,2]-fit2$lpFit[,2]
  supports[,j] = seq(min(ll[-1]),max(ll[-1]),length.out=n)
  monotones[,j] = inverse_jack(Ts=supports[,j],ll=ll,N=N,hd=hd,hr=hr)
  lls[,j] = ll[ (1:n)*floor(N/n) ]
  residuals[,j] = ll[ (1:n)*floor(N/n) ]-Y 
  Weights_ls[[j]] = Weights_star(ll=ll,support = supports[,j],hd=hd)
}

# export data
export = list(monotone=monotones,residual=residuals,support=supports,
              ll=lls,weight=Weights_ls,hr.opt=hrs.opt,hds=hds)
save(export,file = '../estimation/high-dim_2023.RData')

