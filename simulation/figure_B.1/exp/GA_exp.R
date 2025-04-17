library(locpol)
library(foreach)
library(doParallel)
library(Matrix)
library(PLRModels)
library(expm)
library(MASS)
library(mlrv)

#### functions ####
# locally stationay
LS_process1 = function(n,burn = 10000,sigma=1){
  out = rep(NA,n)
  eta = rnorm(burn+n,sd=sigma)
  a = function(t) 0.5-(t-0.5)^2
  H = function(t,i){
    as = a(t)^(0:(burn))
    out.H = as*eta[(i+burn):i]/4
    sum(out.H)
  }
  for (i in 1:n) out[i] = H(i/n,i)
  out
}

long_run_var = function(t,sigma=1){# t is a value
  a = function(t) 0.5-(t-0.5)^2
  out = 0.25/(1-a(t))
  out*sigma
}


# put inverse m_I^-1(t) in the lower bound and exactly compute the integration based on special Kernel function
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
## Ts is vector
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
S_l <- function(Ts,l,hr){ # 
  u1 = (1-Ts)/hr; l1 = -Ts/hr
  u1[u1 >1] = 1
  l1[l1 < -1] = -1
  out = 0.75*(u1^(l+1)/(l+1)-u1^(l+3)/(l+3)) - 0.75*(l1^(l+1)/(l+1)-l1^(l+3)/(l+3))
}

## kernel star involved in SCB
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



#### simulation and estimation once####

# G(t,F_j) = b(t)G(t,F_{j-1})+\xi_j
# b is time-varying coefficient
# \xi_j is independent on Sigma_e
# b(t) = 0.15(0.9+0.1sin(2\pi t))
G_local_stat = function(t,j,Sigma_e_sqrt,burn){ # burn should be j-1000
  if(j <= burn){
    return(0)
  }else{
    error = rexp(nrow(Sigma_e_sqrt))-1
    error = Sigma_e_sqrt %*% error
    out = b_tv(t)*G_local_stat(t,j-1,Sigma_e,burn=burn)+ error
    return(out)
  }
}

LS_process_high = function(n,Burn=500,Sigma_e_sqrt){
  outputs = matrix(NA,nrow = n,ncol = nrow(Sigma_e_sqrt))
  for (i in 1:n) outputs[i,] = G_local_stat(i/n,i,Sigma_e_sqrt,burn = i-Burn)
  return(outputs)
}
b_tv = function(t){
  0.15*(0.9+0.1*sin(2*pi*t))
}
# true long-run covariance
lrv_high = function(n,Sigma_e){
  outputs = list()
  for (i in 1:n) {
    outputs[[i]] = Sigma_e / (1-b_tv(i/n))^2
  }
  outputs
}
GCV.selector = function(Ys){
  p = ncol(Ys)
  n = nrow(Ys)
  hrs.opt = rep(NA,p)
  candidate = matrix(NA,p,2)
  for (j in 1:p) {
    candidate[j,] = rule_of_thumb(Ys[,j],as.matrix(X) )   
    candidate[j,2] = candidate[j,2] * n^(1/6-1/7) # between cn^(-1/4) and cn^(-1/7) origin is cn^(-1/6), here we modify it
    hrs.opt[j] = np.gcv(cbind(Ys[,j],X),
                        h.seq = seq(candidate[j,1],candidate[j,2],length.out=10), estimator = 'LLP')$h.opt
  }
  cbind(candidate[,1],hrs.opt,candidate[,2])
}


#try = simulation_estimation_high(hr=0.1,hd=0.001,p=9,Sigma_e = Sigma_e)
#for (i in 1:p) {
#  plot(X,try$sample[,i],pch=16,col='gray80',cex=0.5)
#  lines(X,try$ll[,i],col='blue',lty=5)
#  lines(try$monotone[,i],try$support[,i],col='red')
#  lines(X,m_targets(X,high_dim_A)[,i])
#}


#### SCB high-dim ####

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
    for (r in max(1,i-3):min(length(candidates),i+3) ) {
      Sigma_r[[r]] = Sigma_candidate[[r]]
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

# weight( not weight star, use known m')
Weights = function(m_N,m_t,dm_t,hd){
  W = matrix(NA,nrow = N,ncol = length(m_t))
  W_star = W
  hd_valid = min(outer(m_N[2:(N+1)],m_t,FUN = function(a,b) abs(a-b))) # avoid hd is too small s.t. W=0 
  #(this actually would not happen in theory since support is a dense interval)
  hd_valid = max(hd,hd_valid)
  for (t in 1:length(m_t)) W[,t] = dm_t[t]*Kernel_ep( (m_N[2:(N+1)]-m_t[t])/hd_valid )/(N*hd_valid)
  return(Matrix(W,sparse = TRUE) ) # most entries of W is in fact 0, save memory with package Matrix
}

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
  abs(t(Ws) %*% Gk)
}


###########################verify GA###################################

Vt_high = function(Kds, n, N, hr, V, Ws){ # GA if V is gaussian r.v.  stochastic expansion if V is error process
  Gk = t(Kds) %*% V
  abs(t(Ws) %*% Gk)
}


######## setting params 
p=27
n=500
high_dim_A = sqrt(1+(1:p)/p)
m_targets = function(X,high_dim_A){
  p = length(high_dim_A)
  ms = matrix(NA,nrow=length(X),ncol = p)
  for (i in 1:floor((p/3)) ) ms[,i] = high_dim_A[i]*(0.5*X^2+X)
  for (i in (1+floor((p/3))):floor((2*p/3)) ) ms[,i] = high_dim_A[i]*2*log(X+1)
  for (i in (1+floor((2*p/3))):p ) ms[,i] = high_dim_A[i]*exp(X)
  return(ms)
}
# the derivative
dm_targets = function(X,high_dim_A){
  p = length(high_dim_A)
  ms = matrix(NA,nrow=length(X),ncol = p)
  for (i in 1:floor((p/3)) ) ms[,i] = high_dim_A[i]*(X+1)
  for (i in (1+floor((p/3))):floor((2*p/3)) ) ms[,i] = high_dim_A[i]*2/(X+1)
  for (i in (1+floor((2*p/3))):p ) ms[,i] = high_dim_A[i]*exp(X)
  return(ms)
}

Sigma_block = function(sigma){
  Sigma_e = matrix(NA,p,p)
  for (i in 1:p) {
    if (i <= p/2){
      for (j in 1:p) {
        if(j <= p/2){
          Sigma_e[i,j]=sigma*(-0.95)^(abs(i-j))
        }else{
          Sigma_e[i,j]=0
        }
      }
    }else{
      Sigma_e[i,]=0
    }
    Sigma_e[i,i] = sigma
  }
  Sigma_e = (t(Sigma_e)+Sigma_e)/2
  Sigma_e
}
sigma=1
Sigma_e = Sigma_block(sigma)
Sigma = lrv_high(n,Sigma_e)

######################### simulation by parallel ###############################
simulations = 400

time1 = Sys.time() # this will take around 1.6 hour
# parallel computation
cl <- makeCluster(10) # used cpu cores
registerDoParallel(cl)
# start computation
GAs <- foreach(k = 1:simulations , .packages = c('locpol','MASS','expm','Matrix'),
               .errorhandling = "remove") %dopar% {
                 set.seed(k)
                 hr = 0.25
                 hd = hr^(4/3)
                 N = 4000
                 xevals = seq(0,1,by=1/N)
                 grids = (1:n)/n
                 LS = LS_process_high(n,Sigma_e_sqrt=sqrtm(Sigma_e))
                 Weights_ls = list()
                 for (j in 1:p) {
                   m_N = m_targets(xevals,high_dim_A)[,j]
                   m_t = m_targets(grids,high_dim_A)[,j]
                   dm_t = dm_targets(grids,high_dim_A)[,j]
                   Weights_ls[[j]] = Weights(m_N,m_t,dm_t,hd)
                 }
                 Vs = Gaussian_gen_high(Sigma)
                 Kds_stat = K_star_tilde(n=n,N=N,hr=hr)
                 VV.hd = EE.hd = rep(NA,p)
                 for (i in 1:p) {
                   index.hd <- grids >= hd*log(1/hd) & grids <= 1-hd*log(1/hd)
                   VV.hd[i] = max( Vt_high(Kds_stat,n,N,hr,unlist(lapply(Vs, "[[",i)),Ws=Weights_ls[[i]][,index.hd]) )
                   EE.hd[i] = max( Vt_high(Kds_stat,n,N,hr,LS[,i],Ws=Weights_ls[[i]][,index.hd]) )
                 }
                 c(VV.hd,EE.hd)
               }
time2 <- Sys.time()
# time cost
print(time2-time1)
# close cluster
stopImplicitCluster()
stopCluster(cl)

# arrange as dataframe
GAs_df = matrix(NA,nrow = simulations,ncol=2)
for (i in 1:2) GAs_df[,i] = unlist(lapply(GAs,"[[",i))

print(GAs_df)
write.csv(GAs_df,'GA_verify_exp.csv')


