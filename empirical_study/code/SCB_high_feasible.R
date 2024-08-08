library(foreach)
library(doParallel)
library(rstudioapi)
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


## MV
se2 <- function(k,Sigmas){
  r = length(Sigmas)
  Sigma_bar = Sigmas[[1]][[k]]
  for (i in 2:r) {
    Sigma_bar = Sigma_bar+Sigmas[[r]][[k]]
  }
  Sigma_bar = Sigma_bar/r
  se2 = diag((Sigmas[[1]][[k]]-Sigma_bar) %*% (Sigmas[[1]][[k]]-Sigma_bar))
  for (i in 2:r) {
    se2= se2 + (Sigmas[[i]][[k]]-Sigma_bar) %*% (Sigmas[[i]][[k]]-Sigma_bar)
  }
  se2 = sum(diag(se2))/(r-1)
  return(se2)
}

MV_choice <- function(residuals,candidates){# candidates in increasing order
  Sigma_candidate = list()
  for (i in 1:length(candidates)) {
    Sigma_candidate[[i]] = id_structure_var(residuals=residuals,m=candidates[i])
  }
  MV = c()
  for (i in 1:(length(candidates)) ) {
    Sigma_r = list()
    j=1
    for (r in max(1,i-3):min(length(candidates),i+3) ) {
      Sigma_r[[j]] = Sigma_candidate[[r]]
      j=j+1
    }
    se2s = c()
    for (k in max(candidates):n) {
      se2s = c(se2s,se2(k,Sigma_r) )
    }
    MV = c(MV,max(se2s))
  }
  return(candidates[which.min(MV)])
}


id_structure_var = function(residuals,m){ # residuals n*p matrix
  residuals = as.matrix(residuals)
  n = nrow(residuals)
  X = (1:n)/n
  
  delta = list()
  delta[[1]] = residuals[1,] %*% t(residuals[1,])/(m+1)
  for (i in 2:n) {
    Qs = colSums(residuals[max(1,i-m):i,]) 
    delta[[i]] = Qs %*% t(Qs)/(m+1)
  } 
  return(delta)
}
#### functions for bootstrap ####
# Generate p*n matrix 
Gaussian_gen_high_j = function(Sigma_j){ # Sigma_j is covariance
  p = nrow(Sigma_j)
  V = mvrnorm(n=1,mu=rep(0,p),Sigma = Sigma_j)
  return(V)
}
Gaussian_gen_high = function(Sigma){
  lapply(Sigma, Gaussian_gen_high_j)
}

bootstrap_high = function(Kds, n, N, hr, V, Ws){
  Gk = t(Kds) %*% V
  max(abs(t(Ws) %*% Gk),na.rm = NA)
}


#### import data ####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load('../estimation/high-dim_2023.RData')

monotones = export$monotone
residuals = export$residual
supports = export$support
hrs.opt = export$hr.opt
Weights_ls = export$weight
hds = export$hds
remove(export)

B=2000
n = nrow(monotones)
X = (1:n)/n
N = 4000
m_star = MV_choice(residuals = residuals,candidates = 10:250  )
hat_Sigma_var = id_structure_var(residuals = residuals,m = m_star )
Kds_stat_p = list()
for (p in 1:ncol(monotones)) {
  hr = hrs.opt[p]
  Kds_stat_p[[p]] = K_star_tilde(n=n,N=N,hr=hr)
}




# This may take 30 min
cl <- makeCluster(10) # 调用核心数量
registerDoParallel(cl)
# 启动并行计算
time1 <- Sys.time() 
boots.all <- foreach(b = 1:B , .packages = c('locpol','MASS','expm'), .errorhandling = "remove") %dopar% {
  set.seed(b)
  Vs = Gaussian_gen_high(hat_Sigma_var)
  boots = rep(NA,ncol(monotones)) 
  
  for (p in 1:ncol(monotones)) {
    hr = hrs.opt[p]
    hd = hds[p]
    lower.all = max(monotones[1,],hd*log(1/hd))
    upper.all = min(monotones[n,],1-hd*log(1/hd))
    Kds_stat = Kds_stat_p[[p]]
    index <- monotones[,p] >= lower.all & monotones[,p] <= upper.all
    boots[p] = max(bootstrap_high(Kds=Kds_stat,n=n,N=N,hr=hr,
                                  V = unlist(lapply(Vs, "[[",p)),Ws=Weights_ls[[p]][,index]))
  }
  c(max(boots,na.rm = T))
}
time2 = Sys.time()
print(time2-time1)
# 在计算结束后关闭集群
stopImplicitCluster()
stopCluster(cl)

SCB90 = quantile(unlist(lapply(boots.all, "[[",1)),probs=0.9)
SCB95 = quantile(unlist(lapply(boots.all, "[[",1)),probs=0.95)
print(c(SCB90,SCB95))

write.csv(unlist(boots.all),file='../estimation/bootstrap_samples_2023.csv')


