library(locpol)
library(foreach)
library(doParallel)
library(Matrix)
library(PLRModels)
library(expm)
library(MASS)
# args containing parameters: 1. index for batch (1-10) 2. model ("a" or "b") 3. sample size 4. dimension 5. candidate 6. bandwidths 7. simulation per batch
args = commandArgs(trailingOnly=TRUE)
#### functions ####

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
inverse_jack <- function(Ts, ll, N, hd){ # Ts support set, ll local linear estimation
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

interpolate = function( ts,xs,ys){
  out = rep(NA,length(ts))
  for (i in 1:length(out)) {
    out[i] = mean(ys[which.min(abs(xs-ts[i]))])
  }
  out
}





#### simulation and estimation once####
# locally stationary
# G(t,F_j) = b(t)G(t,F_{j-1})+\xi_j
# b is time-varying coefficient
# \xi_j is independent on Sigma_e
# b(t) = 0.15(0.9+0.1sin(2\pi t))
b_tv = function(t){
  0.15*(0.9+0.1*sin(2*pi*t))
}
G_local_stat = function(t,j,Sigma_e,burn){ # burn should be j-1000
  if(j <= burn){
    return(0)
  }else{
    out = b_tv(t)*G_local_stat(t,j-1,Sigma_e,burn=burn)+mvrnorm(n=1,mu=rep(0,nrow(Sigma_e)),Sigma = Sigma_e )
    return(out)
  }
}

LS_process_high = function(n,Burn=500,Sigma_e){
  outputs = matrix(NA,nrow = n,ncol = nrow(Sigma_e))
  for (i in 1:n) outputs[i,] = G_local_stat(i/n,i,Sigma_e,burn = i-Burn)
  return(outputs)
}

# piecewise locally stationary
# Xi = G0(ti,Fi) for ti <= 1/3; Xi=G1(ti,Fi) for ti>1/3
# G0(t,F_j) = 0.5G0(t,F_{j-1})+\xi_j
# G1(t,F_j) = -0.5G1(t,F_{j-1})+\xi_j
# \xi_j is independent on Sigma_e
G0_plocal_stat = function(t,j,Sigma_e,burn){ # burn should be j-500
  if(j <= burn){
    return(0)
  }else{
    out = 0.5*G0_plocal_stat(t,j-1,Sigma_e,burn=burn)+mvrnorm(n=1,mu=rep(0,nrow(Sigma_e)),Sigma = Sigma_e )
    return(out)
  }
}
G1_plocal_stat = function(t,j,Sigma_e,burn){ # burn should be j-500
  if(j <= burn){
    return(0)
  }else{
    out = -0.5*G1_plocal_stat(t,j-1,Sigma_e,burn=burn)+mvrnorm(n=1,mu=rep(0,nrow(Sigma_e)),Sigma = Sigma_e )
    return(out)
  }
}
#G_plocal_stat(1,n,Sigma_e,burn = n-500)
PLS_process_high = function(n,Burn=500,Sigma_e){
  outputs = matrix(NA,nrow = n,ncol = nrow(Sigma_e))
  for (i in 1:floor(n/3)) outputs[i,] = G0_plocal_stat(i/n,i,Sigma_e,burn = i-Burn)
  for (i in (floor(n/3)+1):n) outputs[i,] = G1_plocal_stat(i/n,i,Sigma_e,burn = i-Burn)
  return(outputs)
}

simulation_estimation_high_zhou = function(hr,hd,p,Sigma_e){
  #### simulation+estimation ####
  Ys = matrix(NA,nrow = n,ncol = p)
  # local linear
  lls = matrix(NA,nrow = n,ncol = p)
  # residuals between Y and ll
  residuals = matrix(NA,nrow = n,ncol = p)
  # error between monotone estimator and true function
  errors.hr = rep(NA,p)
  # generate data
  LS_p = error_process_high(n,Sigma_e = Sigma_e)
  Ys = m_targets(X,high_dim_A)+LS_p
  
  for (j in 1:p) {
    Y = Ys[,j]
    XY_data <- as.data.frame(cbind(X,Y))
    colnames(XY_data) = c('X','Y')
    fit1 = locpol(Y~X,data=XY_data,deg = 1,bw = hr/sqrt(2),kernel = EpaK, xeval = X )
    fit2 = locpol(Y~X,data=XY_data,deg = 1,bw = hr,kernel = EpaK,xeval = X )
    ll = 2*fit1$lpFit[,2]-fit2$lpFit[,2]
    lls[,j] = ll
    residuals[,j] = ll-Y 
  }
  # error needs considering the support set
  index.hr = X >= hr & X <= 1-hr
  for (j in 1:p) {
    errors.hr[j] = max(abs(m_targets(X[index.hr],high_dim_A)[,j] - lls[index.hr,j]))
  }
  out.list = list(sample=Ys,ll=lls,residual=residuals,error.hr = errors.hr)
  return(out.list)
} 

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


# weight
Weights_star = function(ll,support,hd){
  W = matrix(NA,nrow = N,ncol = length(support))
  W_star = W
  hd_valid = min(outer(ll[2:(N+1)],support,FUN = function(a,b) abs(a-b))) # avoid hd is too small s.t. W=0 
  #(this actually would not happen in theory since support is a dense interval)
  hd_valid = max(hd,hd_valid)
  for (t in 1:length(support)) W[,t] = Kernel_ep( (ll[2:(N+1)]-support[t])/hd_valid )
  for (t in 1:length(support)) W_star[,t] = W[,t]/sum(W[,t])
  return(Matrix(W_star,sparse = TRUE) ) # most entries of W is in fact 0, save memory with package Matrix
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
bootstrap_high = function( n, hr, V){
  Kds = K_star_tilde(n,n,hr)
  Gk = t(Kds) %*% V
  abs(Gk)
}
interpolate = function( ts,xs,ys){
  out = rep(NA,length(ts))
  for (i in 1:length(out)) {
    out[i] = mean(ys[which.min(abs(xs-ts[i]))])
  }
  out
}

# conservative inmproving SCB
## rearrange
SCB_rearrange = function(SCB.zhou.u,SCB.zhou.l,N,hr,xevals,hd=0.01){ # SCB.zhou.u/l is the upper and lower bound
  Ts.u = seq(min(SCB.zhou.u),max(SCB.zhou.u),length.out = length(xevals) )
  upper.zhou.m.inv = inverse_jack(Ts.u,SCB.zhou.u,N,hd)
  Ts.l = seq(min(SCB.zhou.l),max(SCB.zhou.l),length.out=length(xevals))
  lower.zhou.m.inv = inverse_jack(Ts.l,SCB.zhou.l,N,hd)
  u_l.zhou.m = Ts.u- interpolate(ts=upper.zhou.m.inv,xs=lower.zhou.m.inv,ys=Ts.l)
  return(min(u_l.zhou.m[upper.zhou.m.inv > hr & upper.zhou.m.inv < 1-hr])/2)
}
## iso
SCB_iso = function(SCB.zhou.u,SCB.zhou.l,hr){ # SCB.zhou.u/l is the upper and lower bound
  bound.zhou = X > hr & X < 1-hr
  upper.zhou.iso = isoreg(SCB.zhou.u)$yf
  lower.zhou.iso = isoreg(SCB.zhou.l)$yf
  width.zhou.iso = min(upper.zhou.iso[bound.zhou] - lower.zhou.iso[bound.zhou])
  return(width.zhou.iso/2)
}
  
#### Batch simulations ####
######## params
index = as.integer(args[1])
model = args[2]
n <- as.integer(args[3])  # Sample size
p <- as.integer(args[4])   # Dimension
candidate_m <- eval(parse(text=args[5]))  # Candidate vector (e.g., 10:100)
hrs  <- as.numeric(strsplit(args[6], " ")[[1]])  # Bandwidth vector (split string into numeric vector)
simulations <- as.integer(args[7])    # Simulation parameter (numeric value)

X = (1:n)/n
N = 5000
xevals = seq(0,1,by=1/N)
B = 2000


######## setting mean functions
high_dim_A = sqrt(1+(1:p)/p)
m_targets = function(X,high_dim_A){
  p = length(high_dim_A)
  ms = matrix(NA,nrow=length(X),ncol = p)
  for (i in 1:floor((p/3)) ) ms[,i] = high_dim_A[i]*(0.5*X^2+X)
  for (i in (1+floor((p/3))):floor((2*p/3)) ) ms[,i] = high_dim_A[i]*2*log(X+1)
  for (i in (1+floor((2*p/3))):p ) ms[,i] = high_dim_A[i]*exp(X)
  return(ms)
}

# setting errors high dim
if (model =="a") {
   error_process_high = LS_process_high
}else if (model =="b") {
   error_process_high = PLS_process_high
}


# setting errors high dim
# G(t,F_j) = b(t)G(t,F_{j-1})+\xi_j
# b is time-varying coefficient
# \xi_j is independent on Sigma_e
# b(t) = 0.15(0.9+0.1sin(2\pi t))
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




bandwidth_choices = data.frame(hr=hrs)
SCBs_compare = bandwidth_choices
SCBs_compare$max_width90.zhou = NA
SCBs_compare$max_width95.zhou = NA
SCBs_compare$max_width90.zhou.m = NA
SCBs_compare$max_width95.zhou.m = NA
SCBs_compare$max_width90.zhou.iso = NA
SCBs_compare$max_width95.zhou.iso = NA


for (h in 1:length(hrs)) {
  hr = hrs[h]

  time1 = Sys.time()
  # parallel computation
  cl <- makeCluster(10) # used cpu cores
  registerDoParallel(cl)
  # start computation
  SCBs <- foreach(k = 1:simulations , .packages = c('locpol','MASS','expm','Matrix'),
                   .errorhandling = "remove") %dopar% {
                     set.seed(index*1000+k) # to ensure no random seed overlaps for each batch and core 
                     gen_simul = simulation_estimation_high_zhou(hr=hr,p=p,Sigma_e = Sigma_e)
                     m_star = MV_choice(residuals = gen_simul$residual,candidates = candidate_m) 
                     hat_Sigma = id_structure_var(residuals = gen_simul$residual,m = m_star )
                     output = rep(NA,6)
                     
                     max_boots_high.hr = rep(NA,B)
                     
                     for (b in 1:B) {
                       Vs = Gaussian_gen_high(hat_Sigma)
                       boots.hr = rep(NA,p)
                       for (i in 1:p) {
                         index.hr <- X >= hr & X <= 1-hr
                         boots.hr[i] = max(bootstrap_high(n=n,hr=hr,V = unlist(lapply(Vs, "[[",i)))[index.hr] )
                       }
                       max_boots_high.hr[b] = max(boots.hr,na.rm = T)
                     } 
                     # Zhou SCB
                     output[1:2] = quantile(max_boots_high.hr,probs=c(0.9,0.95))
                     upper90.zhou = gen_simul$ll+output[1]
                     lower90.zhou = gen_simul$ll-output[1]
                     upper95.zhou = gen_simul$ll+output[2]
                     lower95.zhou = gen_simul$ll-output[2]
                     
                     # improve Zhou SCB as conservative SCB
                     width90.m = width95.m = width90.iso = width95.iso  = rep(NA,p)
                     for (i in 1:p) {
                       width90.m[i] = SCB_rearrange(upper90.zhou[,i],lower90.zhou[,i],N,hr,xevals)
                       width95.m[i] = SCB_rearrange(upper95.zhou[,i],lower95.zhou[,i],N,hr,xevals)
                       width90.iso[i] = SCB_iso(upper90.zhou[,i],lower90.zhou[,i],hr)
                       width95.iso[i] = SCB_iso(upper95.zhou[,i],lower95.zhou[,i],hr)
                     }
                     output[3:6] = c(max(width90.m),max(width95.m),max(width90.iso),max(width95.iso))
                     output
                     
                     
                     
                   }
  time2 <- Sys.time()

  print(time2-time1)
  
  # close cluster
  stopImplicitCluster()
  stopCluster(cl)
  
  # SCBs dataframe
  SCBs_compare[h,2:7] = Reduce("+",SCBs)/length(SCBs)
  print(SCBs_compare[h,])
}


coverage_name = paste0(index,'batch',B,'boot',model,simulations,'p',p,'sigma',sigma,'zhou')
write.csv(SCBs_compare,file = paste0(coverage_name,'.csv'))




