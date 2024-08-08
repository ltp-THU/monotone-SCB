library(foreach)
library(doParallel)
library(Matrix)
library(PLRModels)
library(expm)
library(MASS)
#### functions ####
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



#### functions for bootstrap ####
Gaussian_gen = function(Sigma,n,Kds){ 
  # Sigma is the long run variance over n
  # Kds is a n*N matrix obtained from K_star_tilde(n,N,hr)
  V = rnorm(n)
  out = t(Kds) %*% (Sigma * V)
  out
}

# bootstrap for once: bootstrap a n*1 vector
bootstrap = function(Kds, n, N, hr, Sigma, ...){
  Gk = Gaussian_gen(Sigma=Sigma,n=n,Kds=Kds)
  Wi = Weights_star(...)
  abs(t(Wi) %*% Gk)
}

bootstrap.zhou = function(n,hr,xevals){
  V = rnorm(n)
  hj = hr/sqrt(2)
  G_p = rep(NA,length(xevals))
  X = (1:n)/n
  for (i in 1:length(xevals) ) {
    G_p[i] = sum(( 2*Kernel_ep( (X-xevals[i])/hj )/(n*hj) - Kernel_ep( (X-xevals[i])/hr )/(n*hr) )*V )
  }
  max(abs(G_p[xevals > hr & xevals < 1-hr]))
}




#### bootstrap ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load('../estimation/regression_sun.RData')

n = length(sun_reg_estimations[[1]]$monotone)
X = (1:n)/n
N = 4000
xevals = seq(0,1,by=1/N)
boot_sample = list()
B = 2000

Sigma_ls = list()

# this may take 10 min
for (i in c(1,2)) {
  data_ls = sun_reg_estimations[[i]]
  tmax = data_ls$data[,1]
  sun = data_ls$data[,3]
  hr = data_ls$hr.gcv
  hd = hr^3
  Kds_stat = K_star_tilde(n=n,N=N,hr=hr)
  bound.m = data_ls$monotone > hr & data_ls$monotone < 1-hr
  time1 = Sys.time()
  boots.all = matrix(NA,nrow = B,ncol = 3)
  for (b in 1:B) {
    set.seed( 5000*(i-1)+b )
    Gk = Gaussian_gen(Sigma=data_ls$Sigma,n=n,Kds=Kds_stat)
    boots = abs(t(data_ls$weight) %*% Gk)
    max_boots = max(boots,na.rm = TRUE)
    max_boots.b = max(boots[bound.m],na.rm = TRUE)
    max_boots.zhou = bootstrap.zhou(n=n,hr=hr,xevals = seq(0,1,by=1/N))
    boots.all[b,] = c(max_boots,max_boots.b,max_boots.zhou)
  }
  time2 = Sys.time()
  print(time2-time1)
  # SCBs
  boot_sample[[i]] = boots.all
}
save(boot_sample,file = '../estimation/boot_reg.RData')
