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


LS_process1 = function(n,burn = 10000,sigma=1){
  out = rep(NA,n)
  eta = rnorm(burn+n,sd=sigma)
  a = function(t) 0.5-(t-0.5)^2
  H = function(t,i){
    as = a(t)^(0:(burn))
    out.H = as*eta[(i+burn):i]
    sum(out.H)
  }
  for (i in 1:n) out[i] = H(i/n,i)
  out
}

LS_process_X = function(n,burn = 10000,sigma=1){
  out = rep(NA,n)
  eta = rnorm(burn+n,sd=sigma)
  a = function(t) 0.25 + t/2
  H = function(t,i){
    as = a(t)^(0:(burn))
    out.H = as*eta[(i+burn):i]
    sum(out.H)
  }
  for (i in 1:n) out[i] = H(i/n,i)
  out
}


#### functions for bootstrap ####
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
# estimate long-run covariance matrix

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


###################GCV fot regression#####################
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


#### simulation + estimation ####

# m function
m_target_1 = function(x){
  0.5*x^2 + x
}

m_target_2 = function(x){
  exp(x)
}

# simulation once
simulation_estiamtion_reg_GCV = function(sigma){
  # simulation+estimation #
  LS1 = LS_process1(n = n,sigma = sigma)
  LSX = LS_process_X(n = n,sigma = sigma)
  
  
  Ys = m_target_1(X)+ m_target_2(X)*LSX+LS1
  XX = cbind(rep(1,n),LSX)
  hr_u_o_l = GCV.selector(Ys,XX)
  hr = hr_u_o_l[2]
  hd = min(max(hr^(4/3), (n*hr)^(-1/3),hr^2*n^0.01),hr/n^0.005 ) # optimum hd hd^2=R_n/h_d, meanwhile ensure hd is between (hr^2,hr)
  
  fit1 = sapply(xevals, FUN = tv_beta,y_dat = Ys,x_dat = XX, hr=hr/sqrt(2))
  fit2 = sapply(xevals, FUN = tv_beta,y_dat = Ys,x_dat = XX, hr=hr)
  ll = t(2*fit1)-t(fit2)
  y.hat = ll[(1:n)*floor(N/n),1]+ll[(1:n)*floor(N/n),2] * LSX
  residuals = y.hat-Ys
  residuals_m = XX*residuals
  # long-run variance estimation
  m_star = MV_choice(residuals = residuals_m,candidates = 8:100)
  hat_Sigma = id_structure_var(residuals = residuals_m,m=m_star)
  
  
  hat_M_ts = list()
  for (j in 1:n) {
    hat_M_ts[[j]] = hat_M(X[j],x_dat = XX,hr=hr)
  }
  # monotone
  # monotone estimation
  ll.m1 = seq( min(ll[-1,1]), max(ll[-1,1]), length.out=n)
  ll.m2 = seq( min(ll[-1,2]), max(ll[-1,2]), length.out=n)
  monotone.m1 = inverse_jack(Ts=ll.m1,ll=ll[,1],N=N,hd=hd,hr=hr)
  monotone.m2 = inverse_jack(Ts=ll.m2,ll=ll[,2],N=N,hd=hd,hr=hr)
  weights.m1 = Weights_star(ll=ll[,1],support = ll.m1,hd=hd)
  weights.m2 = Weights_star(ll=ll[,2],support = ll.m2,hd=hd)
  # errors
  errors = matrix(NA,nrow = n,ncol = 2)
  errors[,1] = abs(m_target_1(monotone.m1)-ll.m1)
  errors[,2] = abs(m_target_2(monotone.m2)-ll.m2)
  
  
  list(sample = cbind(Ys,XX), support = cbind(ll.m1,ll.m2), hr=hr,hd=hd,hr.gcv = hr_u_o_l,
       monotone = cbind(monotone.m1,monotone.m2),ll = ll,error = errors,
       weight =  list(weight.m1= weights.m1,weight.m2 = weights.m2),
       Sigma = hat_Sigma, M = hat_M_ts)
} 

#time1 = Sys.time()
#try = simulation_estiamtion_reg(hr=hr,hd=hd,sigma = sigma)
#time2 = Sys.time()
#print(time2-time1)

#plot(X,try$sample[,1],pch=16,col='gray80',cex=0.5)
#plot(xevals,try$ll[,1],'l')
#lines(try$monotone[,1],try$support[,1],col='blue')
#lines(X,m_target_1(X),'l')
#plot(xevals,try$ll[,2],'l',col='red')
#lines(try$monotone[,2],try$support[,2],col='blue')
#lines(X,m_target_2(X),'l')

#lines(try$monotone[,1],seq(min(try$ll[,1]),max(try$ll[,1]),length.out=n),col='red')
#lines(X,m_target(X))


#### SCB+simulation ####

n = 500

N = 4000
X = (1:n)/n
xevals = seq(0,1,by=1/N)
sigma = 1

B = 2000
simulations = 40



time1 = Sys.time()
# parallel computation
cl <- makeCluster(10) # used cpu cores
registerDoParallel(cl)
# start computing
time1 <- Sys.time() 
cover <- foreach(k = 1:simulations , .packages = c('locpol','MASS','expm','mlrv'), .errorhandling = "remove") %dopar% {
  set.seed(args*1000+k)
  gen_simul = simulation_estiamtion_reg_GCV(sigma = sigma)
  hr = gen_simul$hr
  hd = gen_simul$hd
  Kds_stat = K_star_tilde(n=n,N=N,hr=hr)
  output = rep(NA,7)
  Sigma.m12 = list()
  for (j in 1:n) {
    Sigma.m12[[j]] = solve(gen_simul$M[[j]])  %*% ( (gen_simul$Sigma[[j]]) ) %*% solve(gen_simul$M[[j]])
  }
  # Bootstrap
  max_boots_high = rep(NA,B)
  lower.all = max(gen_simul$monotone[1,],hd*log(1/hd))
  upper.all = min(gen_simul$monotone[n,],1-hd*log(1/hd))
  for (b in 1:B) {
    ### m12
    Vs = Gaussian_gen_high(Sigma = Sigma.m12)
    boots = rep(NA,2)
    index = matrix(NA,nrow = n,ncol = 2)
    for (i in 1:2) {
      index[,i] <- gen_simul$monotone[,i] >= lower.all & gen_simul$monotone[,i] <= upper.all
      boots[i] = max(bootstrap_high(Kds=Kds_stat,n=n,N=N,hr=hr,
                                    V = unlist(lapply(Vs, "[[",i)),Ws=gen_simul$weight[[i]][,index[,i] ]),na.rm=T)
    }
    max_boots_high[b] = max(boots,na.rm = T)
  }
  
  #### cover-m12 ####
  output[3:4] = quantile(max_boots_high,probs=c(0.9,0.95) )
  output[1:2] = max( max(gen_simul$error[index[,1],1]),max(gen_simul$error[index[,2],2]) ) < output[3:4]
  output[5:7] = gen_simul$hr.gcv
  output
}
time2 <- Sys.time()
print(time2-time1)

# close cluster
stopImplicitCluster()
stopCluster(cl)
# output coverage
coverage_df = data.frame(matrix(unlist(cover),nrow=length(cover),ncol = 7,byrow = T))
colnames(coverage_df) = c('is.cover90','is.cover95','SCB90','SCB95','hr.lower','hr.opt','hr.upper')

print(coverage_df)

name = paste0(args,'GCV',B,'boot',simulations,'LS_simul_SCB_reg_est','n',n,'.csv')
write.csv(coverage_df,file=name)









