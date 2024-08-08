library(stringr)
#### import estimation results####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
ts_all.trunc = read.csv('../estimation/data_processed.csv')
ts_all.trunc = ts_all.trunc[,-1]
names = colnames(ts_all.trunc)
load('../estimation/high-dim_2023.RData')
bootstrap_sample = read.csv('../estimation/bootstrap_samples_2023.csv')
monotones = export$monotone
supports = export$support
hrs = export$hr.opt
hds = export$hds
remove(export)

n = nrow(monotones)
X = (1:n)/n
N = 4000

############ SCB for testing quadratic trend ##################3
# ordered by latitude
location = read.csv('../raw_data/locations.csv',encoding = 'UTF-8')
# match location
lat = rep(NA,length(names))
for (i in 1:length(names) ) {
  lat[i] = location$lat[location$name==names[i]]
}


SCB = c(quantile(bootstrap_sample[,2],0.9), quantile(bootstrap_sample[,2],0.95) )
# find p-value
check_p_values_quadra = function(q_fitting,m_fitting){
  qs = seq(0,1,0.001)
  quantiles = quantile(bootstrap_sample[,2],qs)
  1-qs[which.min( abs( max(abs( q_fitting-m_fitting )) - quantiles ) )]
}
# 二次函数拟合
p_values = c()
sites.order = order(lat)

# save png
png('../figures/Figure 2.png',width = 900,height = 500)
par(mfrow=c(3,9),mar=c(0.5,2,2.5,0.5),omi=c(0.1,0.1,0.1,0.1))
for (j in sites.order) {
  index.m = monotones[,j] >= hds[j]*log(1/hds[j]) & monotones[,j] <= 1-hds[j]*log(1/hds[j])
  index.x = X >= hds[j]*log(1/hds[j]) & X <= 1-hds[j]*log(1/hds[j])
  plot(X[index.x],ts_all.trunc[index.x,j],ann=F,xaxt="n",yaxt="n", ylim = c(6,18),
       pch=16,col='gray80',cex=0.5)
  #axis(1,cex.axis=0.8,tcl=-0.2,mgp=c(3,0.2,0))
  axis(2,cex.axis=0.8,tcl=-0.2,mgp=c(3,0.2,0))
  #title(ylab='max_temp.',line=0.8)
  title(main= str_to_title(names[j]),line=0.8)
  #lines(monotones[,j],supports[,j],col='red')
  lines(monotones[index.m,j],supports[index.m,j]+SCB[1],lty='dotted',col='red')
  lines(monotones[index.m,j],supports[index.m,j]-SCB[1],lty='dotted',col='red')
  lines(monotones[index.m,j],supports[index.m,j]+SCB[2],lty='dashed',col='red')
  lines(monotones[index.m,j],supports[index.m,j]-SCB[2],lty='dashed',col='red')
  # 二次函数
  y = ts_all.trunc[index.x,j]
  xx = X[index.x]
  quadratic_model <- nls(y ~ a * xx^2 + b * xx + c, algorithm = 'port',
                         start = list(a = 1, b = 1, c = 1),
                         lower = c(a=0,b=0))
  a = coef(quadratic_model)[1]
  b = coef(quadratic_model)[2]
  c = coef(quadratic_model)[3]
  lines(monotones[index.m,j],a*monotones[index.m,j]^2+b*monotones[index.m,j]+c,col='blue')
  p_values = c(p_values,check_p_values_quadra(supports[index.m,j],a*monotones[index.m,j]^2+b*monotones[index.m,j]+c))
  if(max(a*monotones[index.m,j]^2+b*monotones[index.m,j]+c - supports[index.m,j]-SCB[2])>0) print(names[j])
}
dev.off()


print(paste0('P-value in Figure 2: ', round(min(p_values),3) ) )




########### p-value for testing 0.5 growth delta ################

max_increase = function(delta=0.5,monotones,supports){
  if(delta == 1){
    out = c()
    for (j in 1:ncol(monotones)) {
      index.m = monotones[,j] >= hds[j]*log(1/hds[j]) & monotones[,j] <= 1-hds[j]*log(1/hds[j])
      monotone = monotones[index.m,j]
      support = supports[index.m,j]
      diff = max(support)-min(support)
      out = c(out,diff)
    }
    return(out)
  }else{
    out = c()
    for (j in 1:ncol(monotones)) {
      diff = c() 
      index.m = monotones[,j] >= hds[j]*log(1/hds[j]) & monotones[,j] <= 1-hds[j]*log(1/hds[j])
      monotone = monotones[index.m,j]
      support = supports[index.m,j]
      for (t in monotone[monotone>delta] ) {
        diff = c(diff,support[which.min(abs(monotone-(t) ))]-support[which.min(abs(monotone-(t-delta) ))]  )
      }
      out = c(out,max(diff))
    }
    return(out)
  }
}

check_p_value = function(delta=0.5,C=0.5,bootstrap){
  qs = seq(0,1,0.001)
  quantiles = quantile(bootstrap,probs=qs)
  round(1-qs[which.min( abs(max(max_increase(delta,monotones,supports)-C)/2-quantiles))],3)
}


p_value_all = check_p_value(delta = 1,C=0.5,bootstrap_sample[,2]) 
p_value_delta = check_p_value(delta = 0.75,C=0.5,bootstrap_sample[,2]) 
print(paste0('p-value for (8.2): ',p_value_all))
print(paste0('p-value for (8.4): ',p_value_delta))











