library(sm)
library(Matrix)
library(foreach)
library(doParallel)
################ functions ############################3
# CHETVERIKOV test
test_stat = function(xdata=X,Y,Q_list,sigma_list,S){ # Q is a list containing all choices
  test.out = c()
  Y_matrix = outer(Y,Y,FUN=function(a,b) a-b)
  X_sign_matrix = t(outer(X,X,FUN=function(a,b) as.numeric(a-b>0) ))
  for (i in 1:nrow(S)) {
    bs = 0.5 * sum( Y_matrix * X_sign_matrix * Q_list[[i]] )
    sigma = sigma_list[[ which(Hn==S[i,2]) ]]
    V = sum(sigma^2 * (X_sign_matrix %*% Q_list[[i]])^2)
    test.out = c(test.out,bs/sqrt(V))
  }
  return(max(test.out))
}
Kernel_ep <- function(u){
  support <- as.numeric(abs(u)<=1)
  0.75*(1-u^2)*support
}
# Q for CHETVERIKOV test
Q_cheterikov0_list = function(xdata=X,S){ # S are matrix containing all the choices of Sn
  out = list()
  for (i in 1:nrow(S)) {
    out[[i]] = Matrix(Kernel_ep((xdata-S[i,1])/S[i,2]) %*% t(Kernel_ep((xdata-S[i,1])/S[i,2])), sparse = T)
  }
  out
}

Q_cheterikov1_list = function(xdata=X,S){ # S are matrix containing all the choices of Sn
  out = list()
  for (i in 1:nrow(S)) {
    X_diff = outer(xdata,xdata, FUN = function(a,b) abs(a-b))
    out[[i]] =  Matrix(X_diff * (Kernel_ep((xdata-S[i,1])/S[i,2]) %*% t(Kernel_ep((xdata-S[i,1])/S[i,2]))),sparse = T)
  }
  out
}
# GSV-test
GSV_test = function(s,X,Y){ # s can be vector (t,h) of time point and bandwidth
  n = length(X)
  b = 0
  V = 0
  for (i in 1:n) {
    Vi = 0
    for (j in 1:n) {
      signpart = as.numeric( X[j]-X[i] >0 )*as.numeric(Y[i]-Y[j]>0)
      b = b + signpart*Kernel_ep((X[i]-s[1])/s[2])*Kernel_ep((X[j]-s[1])/s[2])
    }
  }
  return(b)
}


# sigma estimation
sigma_est = function(X,Y,h){
  n = length(X)
  sigma = rep(NA,n)
  for (i in 1:n) {
    Ji <- which(abs(X-X[i]) <= h)
    Ji <- Ji[-length(Ji)]
    sigma[i] = sum((Y[Ji+1]-Y[Ji])^2)/(2*length(Ji))
  }
  return(sqrt(sigma) )
}

# bootstrap
test_stat_boot = function(xdata=X,ydata,Q_list,sigma_list,S){
  n = length(xdata)
  b = 0
  V = 0
  X_sign_matrix = t(outer(X,X,FUN=function(a,b) as.numeric(a-b>0) ))
  noise = rnorm(n)
  test.out = c()
  for (i in 1:nrow(S)) {
    sigma = sigma_list[[ which(Hn==S[i,2]) ]]
    sigma_matrix = outer(sigma*noise ,sigma*noise ,FUN=function(a,b) a-b)
    bs = 0.5 * sum( sigma_matrix * X_sign_matrix * Q_list[[i]] )
    V =  sum(sigma^2 * ( X_sign_matrix %*% Q_list[[i]])^2)
    test.out = c(test.out,bs/sqrt(V))
  }
  return(max(test.out))
}

###################### test ##############################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ts_all.trunc = read.csv('../estimation/data_processed.csv')
ts_all.trunc = ts_all.trunc[,-1]
n = nrow(ts_all.trunc)
X = (1:n)/n


### Bowman test on heathrow data
set.seed(0)
mtest = sm.monotonicity(x=X,y=ts_all.trunc$heathrow,type = "continuous",display="none")


### Chetverikov test on heathrow data
hmax = 0.5
hmin = 0.4*hmax*(log(n)/n)^(1/3)
Hn = c(hmax)
l=1
while (min(Hn) >= hmin) {
  h = hmax * (0.5)^l
  Hn = c(Hn,h)
  l = l+1
}

Sn = as.matrix(expand.grid(X,Hn))
Q_list0 = Q_cheterikov0_list(X,Sn)

sigma_ls = list()
for (h in 1:length(Hn)) sigma_ls[[h]] = sigma_est(X,ts_all.trunc$heathrow,Hn[h]) 

# about 2min
test_stat = test_stat(X,ts_all.trunc$heathrow,Q_list=Q_list0,sigma_list = sigma_ls,Sn) 

# about 1 hour
B=200
cl <- makeCluster(10) # 调用核心数量
registerDoParallel(cl)
# 启动并行计算
time1 <- Sys.time() 
boots <- foreach(b = 1:B , .packages = c('Matrix'), .errorhandling = "remove") %dopar% {
  set.seed(b)
  test_stat_boot(X,ts_all.trunc$heathrow,Q_list=Q_list0,sigma_list=sigma_ls,Sn) 
}
time2 = Sys.time()
print(time2-time1)
# 在计算结束后关闭集群
stopImplicitCluster()
stopCluster(cl)

q.boots = quantile(unlist(boots),probs = seq(0.01,1,0.01) )
p_chet = which.min(abs(q.boots-test_stat))/100

print(paste0('P-values in page 2: ', p_chet, '(Chetverikov test); ', mtest$p, '(Bowman test)'))

