
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load('../estimation/high-dim_2023.RData')
ts_all.trunc = read.csv('../estimation/data_processed.csv')
ts_all.trunc = ts_all.trunc[,-1]






fastlocalt=function(yy,bb,T=1/length(yy),endj=length(yy))###T should be in 1/n, and enj>=floor(nn*bb)
{nn=length(yy)

a1=rep(0,nn)
a2=rep(0,nn)
a3=rep(0,nn)
a4=rep(0,nn)
a5=rep(0,nn)
a6=rep(0,nn)

b1=rep(0,nn)
b2=rep(0,nn)
b3=rep(0,nn)
b4=rep(0,nn)
b5=rep(0,nn)
b6=rep(0,nn)

TT=rep(0,nn)
Tstart=nn*T

TT[1]=T

a1[1]=sum(yy*((1:nn)>=1)*((1:nn)<=floor(Tstart+nn*bb)) )
a2[1]=sum((1:nn)*yy*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
a3[1]=sum((1:nn)^2*yy*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
a4[1]=sum((1:nn)^3*yy*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
#a5[1]=sum((1:nn)^4*yy*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
#a6[1]=sum((1:nn)^5*yy*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )


b1[1]=sum(((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
b2[1]=sum((1:nn)*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
b3[1]=sum((1:nn)^2*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
b4[1]=sum((1:nn)^3*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
b5[1]=sum((1:nn)^4*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
#b6[1]=sum((1:nn)^5*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )

#Lend=rep(-1,500)
#Rend=rep(-1,500)
Lend=-1*floor((nn*bb)-Tstart)
Rend=floor(Tstart+(nn*bb))+1
Endj=endj
for(i in 2: Endj)########The following is for Updating#######
{lend1=(Lend>=1)
rend1= (Rend<=Endj)
lend=Lend*lend1+(Lend<1)
rend=Rend*rend1+(Rend>Endj)
a1[i]=a1[i-1]-yy[lend]*(1)*lend1+yy[rend]*(1)*rend1
a2[i]=a2[i-1]-yy[lend]*(lend)*lend1+yy[rend]*(rend)*rend1
a3[i]=a3[i-1]-yy[lend]*(lend^2)*lend1+yy[rend]*(rend^2)*rend1
a4[i]=a4[i-1]-yy[lend]*(lend^3)*lend1+yy[rend]*(rend^3)*rend1
#     a5[i]=a5[i-1]-yy[lend]*(lend^4)*lend1+yy[rend]*(rend^4)*rend1

b1[i]=b1[i-1]-(1)*lend1+(1)*rend1
b2[i]=b2[i-1]-(lend)*lend1+(rend)*rend1
b3[i]=b3[i-1]-(lend^2)*lend1+(rend^2)*rend1
b4[i]=b4[i-1]-(lend^3)*lend1+(rend^3)*rend1
b5[i]=b5[i-1]-(lend^4)*lend1+(rend^4)*rend1

T=T+1/nn
Rend=Rend+1
Lend=Lend+1
TT[i]=T
# Lend[i]=lend
# Rend[i]=rend
}

Y_1=a1
Y_2=bb^(-1)*(a2/nn-TT*a1)
Y_3=bb^(-2)*(a3/nn^2-2*a2*TT/nn+TT^2*a1)
Y_4=bb^(-3)*(a4/nn^3-3*a3*TT/nn^2+3*a2*TT^2/nn-TT^3*a1)
#  Y_5=bb^(-4)*(a5/nn^4-4*a4*TT/nn^3+6*a3*TT^2/nn^2-4*a2/nn*TT^3+TT^4*a1)



R_1=b1
R_2=bb^(-1)*(b2/nn-TT*b1)
R_3=bb^(-2)*(b3/nn^2-2*b2*TT/nn+TT^2*b1)
R_4=bb^(-3)*(b4/nn^3-3*b3*TT/nn^2+3*b2*TT^2/nn-TT^3*b1)
R_5=bb^(-4)*(b5/nn^4-4*b4*TT/nn^3+6*b3*TT^2/nn^2-4*b2/nn*TT^3+TT^4*b1)


K1=3/4*(R_1-R_3)/nn/bb

K2=3/4*(R_2-R_4)/nn/bb

K3=3/4*(R_3-R_5)/nn/bb


Y1=3/4*(Y_1-Y_3)/nn/bb

Y2=3/4*(Y_2-Y_4)/nn/bb

norm=K1*K3-K2^2
coeff1=(K3*Y1-K2*Y2)/norm
coeff2=(K1*Y2-K2*Y1)/norm
return(list(mean=coeff1[1:endj],derivative=coeff2[1:endj],time=TT[1:endj]))

}


band_matrix=function(nn,bb)
{Ma_temp=function(t)
{XX=cbind(rep(1,nn),(1:nn)/nn-t/nn)
kt=EPK(((1:nn)/nn-t/nn)/bb)
KT=diag(kt)

uuu=(solve(t(XX*kt)%*%XX)%*%t(XX*kt))[1,t]
return(uuu)}
DIAG=unlist(lapply(1:nn, Ma_temp))

return(sum(DIAG))
}

band_matrix_fast=function(nn,bb,Endj=nn)##Endj should be smaller than nn###
{Ma_temp=function(t)
{XX=cbind(rep(1,Endj),(1:Endj)/nn-t/nn)
kt=EPK(((1:Endj)/nn-t/nn)/bb)
x1=kt;y1=rep(1,Endj);y2=((1:Endj)/nn-t/nn)
x2=kt*y2
x11=sum(x1*y1);x12=sum(x1*y2);x21=sum(x2*y1);x22=sum(x2*y2)
scale=1/(x11*x22-x12*x21)
y11=x22;y22=x11;y12=-1*x12;y21=-1*x21
uuu=(matrix(c(y11,y12)*scale,nrow=1)%*%t(XX*kt))[t]
return(uuu)}

DIAG=unlist(lapply(1:Endj, Ma_temp))

return(sum(DIAG))
}


Estimatrix_constau=function(yy,l,tau,NConst=length(yy))
{ Endj=length(yy)
Nvector=rep(0, (NConst-Endj))
yy=c(yy,Nvector)
error=yy[1:Endj]-mean(yy[1:Endj])
Error=matrix(0,nrow=(l+1),ncol=Endj)
Error[1,]=error^2
if(l>=1)
{
  for(ss in 1:l)
  {u2=Endj-ss;u1=ss+1
  Error[u1,1:u2]=error[1:u2]*error[u1:Endj]
  }
}
#length(error[1:(nn-ss)])###Why??
#length(error[(ss+l):nn])
AA=matrix(0,ncol=Endj,nrow=Endj)
#error_compliment=rep(0,(length(yy)))
errortemp=rep(0,length(yy));errortemp[1:Endj]=error;error=errortemp;
diag(AA)=(fastlocalt(error^2,bb=tau,endj=Endj)$mean)[1:Endj]
if(l>=1)
{
  for (ss in 1:l)
  {uu2=Endj-ss;uu1=ss+1
  error_temp=error[1:uu2]*error[uu1:Endj]
  u1=rep(0,Endj)
  u2=rep(0,Endj)
  r1=floor(ss/2)
  r2=ceiling(ss/2)
  r3=r1+1
  u1[(1+r1):(Endj-r2)]=error_temp
  u2[(1+r2):(Endj-r1)]=error_temp
  Tstart=(1+(2+ss)/2-floor((2+ss)/2))/NConst
  gammaplus=fastlocalt(c(u1,Nvector),bb=tau,T=Tstart)$mean[r3:(Endj-ss+r3-1)]
  gammaminus=fastlocalt(c(u2,Nvector),bb=tau,T=Tstart)$mean[r3:(Endj-ss+r3-1)]
  
  AA[row(AA) == (col(AA)-ss)]=(gammaplus+gammaminus)/2
  #AA[row(AA) == (col(AA)+ss)]=(gammaplus+gammaminus)/2
  
  }
}
AA=AA+t(AA)
diag(AA)=diag(AA/2)
return(AA)
}


Estimatrix=function(yy,l,tau=seq(0.1,0.35,by=0.02),NConst=length(yy),adjtau=rep(0,9),QQtau=0)##adjtau is to adjust tau for lags
{ Tau=rep(0,(l+1))
Endj=length(yy)
Nvector=rep(0, (NConst-Endj))
error=yy[1:Endj]-mean(yy[1:Endj])
Error=matrix(0,nrow=(l+1),ncol=Endj)
Error[1,]=error^2
QQ=lapply(tau,function(x) band_matrix_fast(nn=NConst,bb=x,Endj=Endj))
yyy=lapply(tau,function(x) {sum((fastlocalt(yy=c(Error[1,],Nvector),bb=x,endj=Endj)$mean-Error[1,])^2)})
gcvv=unlist(yyy)/(Endj*(1-unlist(QQ)/Endj)^2)
tau_used=tau[order(gcvv)[1]]
Tau[1]=tau_used
if(l>=1)
{
  for(ss in 1:l)
  {u2=Endj-ss;u1=ss+1
  Error[u1,1:u2]=error[1:u2]*error[u1:Endj]
  }
}
#length(error[1:(nn-ss)])###Why??
#length(error[(ss+l):nn])
AA=matrix(0,ncol=Endj,nrow=Endj)
error=c(error,Nvector)
diag(AA)=(fastlocalt(error^2,bb=tau_used,endj=Endj)$mean)[1:Endj]
if(max(abs(QQtau)==0)){
  lindex=l
  QQtau=matrix(0,nrow=length(tau),ncol=lindex)
  
  for (vv in 1:lindex){ tau_ss=pmax(pmin(tau+adjtau[vv],1),0)
  QQtau[,vv]=unlist(lapply(tau_ss,function(x) band_matrix_fast(nn=NConst,bb=x,Endj=Endj)));
  }
}
if(l>=1)
{
  for (ss in 1:l)
  {uu2=Endj-ss;uu1=ss+1
  error_temp=error[1:uu2]*error[uu1:Endj]
  u1=rep(0,Endj)
  u2=rep(0,Endj)
  r1=floor(ss/2)
  r2=ceiling(ss/2)
  r3=r1+1
  u1[(1+r1):(Endj-r2)]=error_temp
  u2[(1+r2):(Endj-r1)]=error_temp
  Tstart=(1+(2+ss)/2-floor((2+ss)/2))/NConst
  tau_ss=pmax(pmin(tau+adjtau[ss],1),0)
  yyy=lapply(tau_ss,function(x) {sum((fastlocalt(yy=c(u1,Nvector),bb=x,endj=Endj)$mean[1:Endj]-u1)^2)})
  gcvv=unlist(yyy)/(Endj*(1-unlist(QQtau[,ss])/Endj)^2)
  tau_used=tau_ss[order(gcvv)[1]]
  Tau[(ss+1)]=tau_used
  gammaplus=fastlocalt(c(u1,Nvector),bb=tau_used,T=Tstart)$mean[r3:(Endj-ss+r3-1)]
  gammaminus=fastlocalt(c(u2,Nvector),bb=tau_used,T=Tstart)$mean[r3:(Endj-ss+r3-1)]
  
  AA[row(AA) == (col(AA)- ss)]=(gammaplus+gammaminus)/2
  #AA[row(AA) == (col(AA)+ss)]=(gammaplus+gammaminus)/2
  
  }
}

AA=AA+t(AA)
diag(AA)=diag(AA/2)

return(list(AA=AA,tau=Tau))
}

# kernel
EPK <- function(u){
  support <- as.numeric(abs(u)<=1)
  0.75*(1-u^2)*support
}

# heathrow is j=12
n = nrow(export$residual)
X = (1:n)/n
hr = export$hr.opt[12]
com_fact = export$residual$heathrow[X>= hr & X<= 1-hr]

vtest1=Estimatrix(com_fact,l=1)
vtest2=Estimatrix(com_fact,l=2)
vtest3=Estimatrix(com_fact,l=3)
c1=0
c2=0
c3=0
TT =length(com_fact)
for (i in 1: (TT-3) ){
  c1[i]=vtest1$AA[i,i+1]/sqrt(vtest1$AA[i,i]*vtest1$AA[i+1,i+1])
  c2[i]=vtest2$AA[i,i+2]/sqrt(vtest1$AA[i,i]*vtest1$AA[i+2,i+2])
  c3[i]=vtest3$AA[i,i+3]/sqrt(vtest1$AA[i,i]*vtest1$AA[i+3,i+3])
}
years = seq(1979,2023,length.out=TT-3)

hd = export$hds[12]
index.m = export$monotone$heathrow>= hd*log(1/hd) & export$monotone$heathrow <= 1-hd*log(1/hd)
# figure 1 (a)
png('../figures/Figure 1(a).png',width=360,height = 440,res=100)
par(mfrow=c(1,1),mar=c(2,4,2.5,2),omi=c(0.1,0.1,0.1,0.1))
years.x = seq(1979,2023,length.out=n)
years.m = 1979 + (2023-1979) * export$monotone$heathrow
plot(years.x,ts_all.trunc$heathrow,ann=F,xaxt="n",yaxt="n",pch=16,col='gray80',cex=0.5)
axis(1,cex.axis=0.8,tcl=-0.2,mgp=c(3,0.2,0))
axis(2,cex.axis=0.8,tcl=-0.2,mgp=c(3,0.2,0))
lines(years.m[index.m],export$support$heathrow[index.m],col='red',lwd=2)
title(xlab='Year',line=-1.2)
dev.off()

# figure 1 (b)
png(filename = '../figures/Figure 1(b).png',width = 550, height =  487,res=100)
par(mfrow=c(3,1),mar=c(2,4,2.5,2),omi=c(0.1,0.1,0.1,0.1))
plot(years, c1,'l',xaxt="n",ann=F)
legend(x=1979,y=0.3, legend = "1st order autocorrelation",bty='n',cex = 1.25)
plot(years, c2,'l',xaxt="n",ann=F,ylim = c(-0.2,0.2))
legend(x=1979,y=0.15, legend = "2nd order autocorrelation",bty='n',cex=1.25)
plot(years, c3,'l',ann=F,ylim = c(-0.2,0.2))
legend(x=1979,y=0.2, legend = "3rd order autocorrelation",bty='n',cex=1.25)
title(xlab='Year',line=-1.2)
dev.off()
















