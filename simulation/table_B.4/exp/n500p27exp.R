library(locpol)
library(foreach)
library(doParallel)
library(Matrix)
library(PLRModels)
library(expm)
library(MASS)
library(mlrv)
args = commandArgs(trailingOnly=TRUE)
args = as.numeric(args)
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

# with many hds
simulation_estimation_high.gcv = function(p,Sigma_e_sqrt){
  #### simulation+estimation ####
  Ys = matrix(NA,nrow = n,ncol = p)
  supports = matrix(NA,nrow = n,ncol = p)
  # local linear
  lls = matrix(NA,nrow = n,ncol = p)
  # residuals between Y and ll
  residuals = matrix(NA,nrow = n,ncol = p)
  # error between monotone estimator and true function
  errors.hr = rep(NA,p)
  errors.hd = rep(NA,p)
  # generate data
  LS_p = LS_process_high(n,Sigma_e = Sigma_e_sqrt)
  Ys = m_targets(X,high_dim_A)+LS_p
  
  hrs_u_o_l = GCV.selector(Ys)
  hr = mean(hrs_u_o_l[,2])
  hd = max(hr^(4/3), (n*hr)^(-1/3))
  # monotone
  monotones = matrix(NA,nrow = n,ncol = p)
  Weights_ls = list()
  for (j in 1:p) {
    Y = Ys[,j]
    XY_data <- as.data.frame(cbind(X,Y))
    colnames(XY_data) = c('X','Y')
    fit1 = locpol(Y~X,data=XY_data,deg = 1,bw = hr/sqrt(2),kernel = EpaK, xeval = xevals )
    fit2 = locpol(Y~X,data=XY_data,deg = 1,bw = hr,kernel = EpaK,xeval = xevals )
    ll = 2*fit1$lpFit[,2]-fit2$lpFit[,2]
    supports[,j] = seq(min(ll[-1]),max(ll[-1]),length.out=n)
    lls[,j] = ll[ (1:n)*floor(N/n) ]
    residuals[,j] = ll[ (1:n)*floor(N/n) ]-Y 
    monotones[,j] = inverse_jack(Ts=supports[,j],ll=ll,N=N,hd=hd,hr=hr)
    Weights_ls[[j]] = Weights_star(ll=ll,support = supports[,j],hd=hd)
  }
  # error needs considering the support set
  lower.all.hd = max(max(monotones[1,]),hd*log(1/hd))
  upper.all.hd = min(min(monotones[nrow(monotones),]),1-hd*log(1/hd))
  for (j in 1:p) {
    index.hd <- monotones[,j] >= lower.all.hd & monotones[,j] <= upper.all.hd
    errors.hd[j] = max(abs(m_targets(monotones[,j],high_dim_A)[index.hd,j] - supports[index.hd,j]))
  }
  out.list = list(sample=Ys,monotone=monotones,support=supports,hr=hr,hd=hd,hr.gcv = hrs_u_o_l,
                  ll=lls,residual=residuals,weight_ls=Weights_ls,error.hd=errors.hd)
  return(out.list)
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
bootstrap_high = function(Kds, n, N, hr, V, Ws){
  Gk = t(Kds) %*% V
  abs(t(Ws) %*% Gk)
}



#### Batch simulations ####
# dims
p = 27
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


######## params
n = 500
X = (1:n)/n
B = 2000
simulations = 40
# setting errors high dim
# G(t,F_j) = b(t)G(t,F_{j-1})+\xi_j
# b is time-varying coefficient
# \xi_j is independent on Sigma_e
# b(t) = 0.15(0.9+0.1sin(2\pi t))

b_tv = function(t){
  0.15*(0.9+0.1*sin(2*pi*t))
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
Sigma_e_sqrt = sqrtm(Sigma_e)

N = 4000 # the computation limit is below 5000, too big N wastes too much memory space
xevals = seq(0,1,by=1/N)


time1 = Sys.time()
# parallel computation
cl <- makeCluster(10) # used cpu cores
registerDoParallel(cl)
# start computation
cover <- foreach(k = 1:simulations , .packages = c('locpol','MASS','expm','Matrix','mlrv','PLRModels'),
                 .errorhandling = "remove") %dopar% {
                   set.seed(args*1000+k) # to ensure no random seed overlaps for each batch and core 
                   gen_simul = simulation_estimation_high.gcv(p=p,Sigma_e_sqrt = Sigma_e_sqrt)
                   hr = gen_simul$hr
                   hd = gen_simul$hd
                   Kds_stat = K_star_tilde(n=n,N=N,hr=hr)
                   m_star = MV_choice(residuals = gen_simul$residual,candidates = 20:100)
                   hat_Sigma = id_structure_var(residuals = gen_simul$residual,m = m_star )
                   output = rep(NA,7)
                   
                   max_boots_high.hd = rep(NA,B)
                   lower.all.hd = max(gen_simul$monotone[1,],hd*log(1/hd))
                   upper.all.hd = min(gen_simul$monotone[n,],1-hd*log(1/hd))
                   for (b in 1:B) {
                     Vs = Gaussian_gen_high(hat_Sigma)
                     boots.hd = rep(NA,p)
                     for (i in 1:p) {
                       index.hd <- gen_simul$monotone[,i] >= lower.all.hd & gen_simul$monotone[,i] <= upper.all.hd
                       boots.hd[i] = max(bootstrap_high(Kds=Kds_stat,n=n,N=N,hr=hr,
                                                        V = unlist(lapply(Vs, "[[",i)),Ws=gen_simul$weight_ls[[i]][,index.hd]))
                     }
                     max_boots_high.hd[b] = max(boots.hd,na.rm = T)
                   }
                   output[3:4] = quantile(max_boots_high.hd,probs=c(0.9,0.95))
                   output[1:2] = max( gen_simul$error.hd ) < output[3:4]
                   output[5:7] = as.numeric(apply(gen_simul$hr.gcv,2,mean))
                   output
                 }
time2 <- Sys.time()

print(time2-time1)

# close cluster
stopImplicitCluster()
stopCluster(cl)

# coverage
coverage_df = data.frame(matrix(unlist(cover),nrow=length(cover),ncol = 7,byrow = T))
colnames(coverage_df) = c('is.cover90','is.cover95','SCB90','SCB95','hr.lower','hr.opt','hr.upper')

# bandwidth selection

coverage_name = paste0(args,'GCV',B,'boot',simulations,'LS',p,'sigma',sigma,'n',n)
write.csv(coverage_df,file = paste0(coverage_name,'.csv'))




