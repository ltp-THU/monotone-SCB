library(locpol)
library(foreach)
library(doParallel)
library(Matrix)
library(PLRModels)
library(expm)
library(MASS)
library(mlrv)
# args containing parameters: 
#1. index for batch (1-10) 
#2. model ("a" or "b") 
#3. sample size 
#4. dimension 
#5. candidate 
#6. bandwidths 
#7. GCV or not 
#8. simulation per batch
#9. penalization C1
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





# with many hds
simulation_estimation_high = function(hr,hd,p,Sigma_e){
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
  LS_p = error_process_high(n,Sigma_e = Sigma_e)
  Ys = m_targets(X,high_dim_A)+LS_p

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
  lower.all.hr = max(max(monotones[1,]),hr)
  upper.all.hr = min(min(monotones[nrow(monotones),]),1-hr)
  lower.all.hd = max(max(monotones[1,]),hd*log(1/hd))
  upper.all.hd = min(min(monotones[nrow(monotones),]),1-hd*log(1/hd))
  for (j in 1:p) {
    index.hr <- monotones[,j] >= lower.all.hr & monotones[,j] <= upper.all.hr
    index.hd <- monotones[,j] >= lower.all.hd & monotones[,j] <= upper.all.hd
    errors.hr[j] = max(abs(m_targets(monotones[,j],high_dim_A)[index.hr,j] - supports[index.hr,j]))
    errors.hd[j] = max(abs(m_targets(monotones[,j],high_dim_A)[index.hd,j] - supports[index.hd,j]))
  }
  out.list = list(sample=Ys,monotone=monotones,support=supports,
                  ll=lls,residual=residuals,weight_ls=Weights_ls,error.hr = errors.hr,error.hd=errors.hd)
  return(out.list)
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
simulation_estimation_high.gcv = function(p,Sigma_e){
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
  LS_p = error_process_high(n,Sigma_e = Sigma_e)
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
######## params
index = as.integer(args[1])
model = args[2]
n <- as.integer(args[3])  # Sample size
p <- as.integer(args[4])   # Dimension
candidate_m <- eval(parse(text=args[5]))  # Candidate vector (e.g., 10:100)
hrs  <- as.numeric(strsplit(args[6], " ")[[1]])  # Bandwidth vector (split string into numeric vector)
GCV <- as.logical(args[7])           # GCV parameter (TRUE or FALSE)
simulations <- as.integer(args[8])    # Simulation parameter (numeric value)

X = (1:n)/n
B = 2000

N = 4000 # the computation limit is below 5000, too big N wastes too much memory space
xevals = seq(0,1,by=1/N)

######## setting mean functions
high_dim_A = sqrt(1+(1:p)/p)
C1 = as.numeric(args[9])
C3 = C1^(19/8)
penal = function(X) C1*X+C3*X^2+C3*X^3
m_targets = function(X,high_dim_A){
  p = length(high_dim_A)
  ms = matrix(NA,nrow=length(X),ncol = p)
  for (i in 1:floor((p/3)) ){
    ms[,i] = (X<1/3)*10*exp(1/3)/9+(X>=1/3)*( (X^2-4*X/3+13/9)*exp(X) )
    ms[,i] = high_dim_A[i]*(ms[,i]+penal(X))
  } 
  for (i in (1+floor((p/3))):floor((2*p/3)) ){
    ms[,i] = (X<1/3)*(1+(3*(X-1/3))^3)+(X>=1/3&X<2/3)*1+(X>=2/3)*(1+(3*(X-2/3))^3)
    ms[,i] = high_dim_A[i]*(ms[,i]+penal(X))
  }
  for (i in (1+floor((2*p/3))):p ){
    ms[,i] = (X<2/3)*2*(cos( 3*pi*(X-2/3)/4))+(X>=2/3)*2
    ms[,i] = high_dim_A[i]*(ms[,i]+penal(X))
  } 
  return(ms)
}

# a help function to check monotonicity
check_monotonicity_increase <- function(data) {# matrix or data frame with p columns
  # Function to check if a single vector is monotonically increasing
  is_monotone_increasing <- function(vec) {
    all(diff(vec) >= 0)
  }
  monotone_results <- apply(data, 2, is_monotone_increasing)
  all_monotone <- all(monotone_results)
  return(all_monotone)
}

#### setting errors high dim
if (model =="a") {
   error_process_high = LS_process_high
}else if (model =="b") {
   error_process_high = PLS_process_high
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






######## parallel simulation ########
if (GCV) {
  time1 = Sys.time()
  # parallel computation
  cl <- makeCluster(10) # used cpu cores
  registerDoParallel(cl)
  # start computation
  cover <- foreach(k = 1:simulations , .packages = c('locpol','MASS','expm','Matrix','mlrv','PLRModels'),
                  .errorhandling = "remove") %dopar% {
                    set.seed(index*1000+k) # to ensure no random seed overlaps for each batch and core 
                    gen_simul = simulation_estimation_high.gcv(p=p,Sigma_e = Sigma_e)
                    hr = gen_simul$hr
                    hd = gen_simul$hd
                    Kds_stat = K_star_tilde(n=n,N=N,hr=hr)
                    m_star = MV_choice(residuals = gen_simul$residual,candidates = candidate_m)
                    hat_Sigma = id_structure_var(residuals = gen_simul$residual,m = m_star )
                    output = rep(NA,8)
                    # monotone shape
                    if (C1==0){
                      output[8] = TRUE
                    }else{
                      hat_y = matrix(NA, nrow = n, ncol = p)
                      for (i in 1:p) {
                        hat_y[,i] = gen_simul$support[,i] - penal(gen_simul$monotone[,i]) 
                        hat_y[,i] = hat_y[,i] + log10(n)*C3*gen_simul$monotone[,i]/sqrt(hd*n*hr)
                      }
                      all_monotone = check_monotonicity_increase(hat_y)
                      output[8] = all_monotone
                    }
                    
                    # SCB coverage
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
                    # error with correction
                    errors.hd = rep(NA,p)
                    for (i in 1:p) {
                      index.hd <- gen_simul$monotone[,i] >= lower.all.hd & gen_simul$monotone[,i] <= upper.all.hd
                      resid_temp = m_targets(gen_simul$monotone[,i],high_dim_A)[index.hd,i] - gen_simul$support[index.hd,i]
                      resid_temp = resid_temp-log10(n)*C3*gen_simul$monotone[,i]/sqrt(hd*n*hr)
                      errors.hd[i] = max(abs(resid_temp))
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
  coverage_df = data.frame(matrix(unlist(cover),nrow=length(cover),ncol = 8,byrow = T))
  colnames(coverage_df) = c('is.cover90','is.cover95','SCB90','SCB95','hr.lower','hr.opt','hr.upper','shape')
  coverage_name = paste0(index,'GCV',B,'boot',model,n,'p',p,'sigma',sigma)
  write.csv(coverage_df,file = paste0(coverage_name,'.csv'))
}else {
    hd.opt = c()
  for (i in hrs) {
    hd.opt = c(hd.opt, min(max(i^(4/3), (n*i)^(-1/3),i^2*n^0.01),i/n^0.005 )) # optimum hd hd^2=R_n/h_d, meanwhile ensure hd is between (hr^2,hr)
  }
  hds = hd.opt
  bandwidth_choices = data.frame(hr=hrs,hd=hds)
  coverage_df = data.frame( matrix(NA,nrow = length(hds),ncol = 8)  )
  colnames(coverage_df) = c('is.cover90.hr','is.cover95.hr','is.cover90.hd','is.cover95.hd',
                      'SCB90.hr','SCB95.hr','SCB90.hd','SCB95.hd')
  coverage_df = cbind(bandwidth_choices,coverage_df)

  for (h in 1:length(hrs)) {
    hr = hrs[h]
    hd = hds[h]
    
    Kds_stat = K_star_tilde(n=n,N=N,hr=hr)
    
    time1 = Sys.time()
    # parallel computation
    cl <- makeCluster(10) # used cpu cores
    registerDoParallel(cl)
    # start computation
    cover <- foreach(k = 1:simulations , .packages = c('locpol','MASS','expm','Matrix'),
                    .errorhandling = "remove") %dopar% {
                      set.seed(index*1000+k) # to ensure no random seed overlaps for each batch and core 
                      gen_simul = simulation_estimation_high(hr=hr,hd=hd,p=p,Sigma_e = Sigma_e)
                      m_star = MV_choice(residuals = gen_simul$residual,candidates = candidate_m) 
                      hat_Sigma = id_structure_var(residuals = gen_simul$residual,m = m_star )
                      output = rep(NA,8)
                      
                      max_boots_high.hr = rep(NA,B)
                      max_boots_high.hd = rep(NA,B)
                      lower.all.hr = max(gen_simul$monotone[1,],hr)
                      upper.all.hr = min(gen_simul$monotone[n,],1-hr)
                      lower.all.hd = max(gen_simul$monotone[1,],hd*log(1/hd))
                      upper.all.hd = min(gen_simul$monotone[n,],1-hd*log(1/hd))
                      for (b in 1:B) {
                        Vs = Gaussian_gen_high(hat_Sigma)
                        boots.hr = rep(NA,p)
                        boots.hd = rep(NA,p)
                        for (i in 1:p) {
                          index.hr <- gen_simul$monotone[,i] >= lower.all.hr & gen_simul$monotone[,i] <= upper.all.hr
                          index.hd <- gen_simul$monotone[,i] >= lower.all.hd & gen_simul$monotone[,i] <= upper.all.hd
                          boots.hr[i] = max(bootstrap_high(Kds=Kds_stat,n=n,N=N,hr=hr,
                                                            V = unlist(lapply(Vs, "[[",i)),Ws=gen_simul$weight_ls[[i]][,index.hr]))
                          boots.hd[i] = max(bootstrap_high(Kds=Kds_stat,n=n,N=N,hr=hr,
                                                            V = unlist(lapply(Vs, "[[",i)),Ws=gen_simul$weight_ls[[i]][,index.hd]))
                        }
                        max_boots_high.hr[b] = max(boots.hr,na.rm = T)
                        max_boots_high.hd[b] = max(boots.hd,na.rm = T)
                      } 
                        output[5:6] = quantile(max_boots_high.hr,probs=c(0.9,0.95))
                        output[7:8] = quantile(max_boots_high.hd,probs=c(0.9,0.95))
                        output[1:2] = max( gen_simul$error.hr ) < output[5:6]
                        output[3:4] = max( gen_simul$error.hd ) < output[7:8]
                      
                      output
                    }
    time2 <- Sys.time()
    # time cost
    print(time2-time1)
    
    # close cluster
    stopImplicitCluster()
    stopCluster(cl)
    
    # coverage dataframe
    coverage_df[h,3:10] = Reduce("+",cover)/length(cover)
    print(coverage_df[h,])
  }
  coverage_name = paste0(index,'batch',B,'boot',model,n,'p',p,'sigma',sigma)
  write.csv(coverage_df,file = paste0(coverage_name,'.csv'))
}







