#########  figure B.3 and its p-values ###########
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load('../estimation/regression_sun.RData')
load('../estimation/boot_reg.RData')
n = length(sun_reg_estimations[[1]]$monotone)
X = (1:n)/n
N = 4000
xevals = seq(0,1,by=1/N)
# one_station: plot zhou's SCB, iso SCB, monotone SCB
SCB_design = function(boot_sample,q){ # generate upper&lower bounds for zhou's SCB, iso SCB, monotone SCB single q
  SCB_zhou_ul = list()
  SCB_iso_ul = list()
  SCB_mon_ul = list()
  for (i in 1:2) {
    data_ls = sun_reg_estimations[[i]]
    # SCB zhou    
    SCB.zhou = data_ls$Sigma * quantile(boot_sample[[i]][,3],q)
    hr = data_ls$hr.gcv
    bound = X>hr & X<1-hr
    bound.N = xevals > hr & xevals < 1-hr
    bound.m = data_ls$monotone > hr & data_ls$monotone < 1-hr
    ll.n = data_ls$ll[ (1:n)*floor(N/n) ,2]
    SCB_zhou_ul[[i]] = cbind(ll.n[bound]+SCB.zhou[bound],ll.n[bound]-SCB.zhou[bound])
    # SCB iso
    y.i.u = isoreg( x = X[bound],y = ll.n[bound]+SCB.zhou[bound])
    y.i.l = isoreg( x = X[bound],y = ll.n[bound]-SCB.zhou[bound])
    SCB_iso_ul[[i]] = cbind(y.i.u$yf , y.i.l$yf)
    # SCB monotone
    SCB.b = quantile(boot_sample[[i]][,2],q)
    SCB_mon_ul[[i]] = cbind(data_ls$support[bound.m]+ SCB.b,data_ls$support[bound.m]- SCB.b )
  }
  return(list(zhou=SCB_zhou_ul,iso=SCB_iso_ul,mon=SCB_mon_ul))
}

SCBs_90 = SCB_design(boot_sample,0.9)
SCBs_95 = SCB_design(boot_sample,0.95)

names = c('Heathrow','Lerwick')
# plot Figure B.3(a)
png('../figures/Figure B.3(a).png',width =400,height = 450,res=100 )
par(mfrow=c(2,1), mai=c(0.5,0.8,0.5,0.1),omi=c(0.1,0.1,0.1,0.1))
for (i in 1:2) {
  ymax = max( SCBs_95$zhou[[i]],SCBs_95$iso[[i]],SCBs_95$mon[[i]] )
  ymin = min( SCBs_95$zhou[[i]],SCBs_95$iso[[i]],SCBs_95$mon[[i]] )
  data_ls = sun_reg_estimations[[i]]
  hr = data_ls$hr.gcv
  bound = X>hr & X<1-hr
  bound.N = xevals > hr & xevals < 1-hr
  ll.n = data_ls$ll[ (1:n)*floor(N/n) ,2]
  
  # Zhou's SCB
  year_evals = seq(1979,2023,length.out=length(xevals))
  year_X = seq(1979,2023,length.out=length(X))
  plot(year_X[bound],ll.n[bound],'l',xaxt="n",yaxt="n",ann=F,ylim = c(ymin*1.1,ymax*1.1))
  axis(1,cex.axis=0.8,tcl=-0.2,mgp=c(3,0.2,0))#,at=c(1990,1995,2000,2005,2010))
  axis(2,cex.axis=0.8,tcl=-0.2,mgp=c(3,0.2,0))
  title(ylab='sunshine coefficient',line=2)
  title(xlab='Year',line=1.5)
  title(main= names[i] ,line=0.8)
  abline(h=0)
  # SCB90
  lines(year_X[bound],SCBs_90$zhou[[i]][,1],lty='dotted')
  lines(year_X[bound],SCBs_90$zhou[[i]][,2],lty='dotted')
  #SCB95
  lines(year_X[bound],SCBs_95$zhou[[i]][,1],lty='dashed')
  lines(year_X[bound],SCBs_95$zhou[[i]][,2],lty='dashed')
}
dev.off()
# plot Figure B.3(b)
png('../figures/Figure B.3(b).png',width =400,height = 450,res=100 )
par(mfrow=c(2,1), mai=c(0.5,0.8,0.5,0.1),omi=c(0.1,0.1,0.1,0.1))
for (i in 1:2) {
  ymax = max( SCBs_95$zhou[[i]],SCBs_95$iso[[i]],SCBs_95$mon[[i]] )
  ymin = min( SCBs_95$zhou[[i]],SCBs_95$iso[[i]],SCBs_95$mon[[i]] )
  data_ls = sun_reg_estimations[[i]]
  hr = data_ls$hr.gcv
  bound = X>hr & X<1-hr
  bound.N = xevals > hr & xevals < 1-hr
  ll.n = data_ls$ll[ (1:n)*floor(N/n) ,2]
  y.i = isoreg( x = xevals[bound.N],y = data_ls$ll[bound.N,2])
  # improved SCB
  year_evals = seq(1979,2023,length.out=length(xevals))
  year_X = seq(1979,2023,length.out=length(X))
  
  plot(year_evals[bound.N],y.i$yf,'l',xaxt="n",yaxt="n",ann=F,ylim =c(ymin*1.1,ymax*1.1),col='blue')
  axis(1,cex.axis=0.8,tcl=-0.2,mgp=c(3,0.2,0))#,at=c(1990,1995,2000,2005,2010))
  axis(2,cex.axis=0.8,tcl=-0.2,mgp=c(3,0.2,0))
  title(ylab='sunshine coefficient',line=2)
  title(xlab='Year',line=1.5)
  title(main= names[i] ,line=0.8)
  abline(h=0)
  # SCB90
  lines(year_X[bound],SCBs_90$iso[[i]][,1],lty='dotted',col='blue')
  lines(year_X[bound],SCBs_90$iso[[i]][,2],lty='dotted',col='blue')
  #SCB95
  lines(year_X[bound],SCBs_95$iso[[i]][,1],lty='dashed',col='blue')
  lines(year_X[bound],SCBs_95$iso[[i]][,2],lty='dashed',col='blue')
}
dev.off()
# plot Figure B.3(c)
png('../figures/Figure B.3(c).png',width =400,height = 450,res=100 )
par(mfrow=c(2,1), mai=c(0.5,0.8,0.5,0.1),omi=c(0.1,0.1,0.1,0.1))
for (i in 1:2) {
  ymax = max( SCBs_95$zhou[[i]],SCBs_95$iso[[i]],SCBs_95$mon[[i]] )
  ymin = min( SCBs_95$zhou[[i]],SCBs_95$iso[[i]],SCBs_95$mon[[i]] )
  data_ls = sun_reg_estimations[[i]]
  hr = data_ls$hr.gcv
  bound.m = data_ls$monotone > hr & data_ls$monotone < 1-hr
  ll.n = data_ls$ll[ (1:n)*floor(N/n) ,2]
  # improved SCB
  year_monotone = 1979 + (2023-1979) * data_ls$monotone
  plot(year_monotone[bound.m],data_ls$support[bound.m],col='red','l',ylim=c(1.1*ymin,1.1*ymax),
       xaxt="n",yaxt="n",ann=F)
  axis(1,cex.axis=0.8,tcl=-0.2,mgp=c(3,0.2,0))#,at=c(1990,1995,2000,2005,2010))
  axis(2,cex.axis=0.8,tcl=-0.2,mgp=c(3,0.2,0))
  title(ylab='sunshine coefficient',line=2)
  title(xlab='Year',line=1.5)
  title(main= names[i] ,line=0.8)
  abline(h=0)
  #lines(data_ls$monotone[bound.m],data_ls$support[bound.m],col='red')
  # SCB90
  lines(year_monotone[bound.m],SCBs_90$mon[[i]][,1],col='red',lty="dotted")
  lines(year_monotone[bound.m],SCBs_90$mon[[i]][,2],col='red',lty="dotted")
  #SCB95
  lines(year_monotone[bound.m],SCBs_95$mon[[i]][,1],col='red',lty="dashed")
  lines(year_monotone[bound.m],SCBs_95$mon[[i]][,2],col='red',lty="dashed")
}
dev.off()


######## compute p-values ########
lower_qs = lapply(seq(0,1,0.001), FUN=SCB_design, boot_sample=boot_sample)

check_p_values = function(station_index=1,qs=lower_qs){
  min_lower = matrix(NA,nrow = length(qs),ncol=3)
  for (k in 1:length(qs)) {
    min_lower[k,1] = min(qs[[k]]$zhou[[station_index]][,2])
    min_lower[k,2] = min(qs[[k]]$iso[[station_index]][,2])
    min_lower[k,3] = min(qs[[k]]$mon[[station_index]][,2])
  }
  
  p.zhou = 1 - (which.min(abs(min_lower[,1]))-1)/1000
  p.iso = 1 - (which.min(abs(min_lower[,2]))-1)/1000
  p.mon = 1 - (which.min(abs(min_lower[,3]))-1)/1000
  return(c(p.zhou,p.iso,p.mon))
}

check_p_values(1,lower_qs)
check_p_values(2,lower_qs)
