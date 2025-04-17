setwd(dirname(rstudioapi::getSourceEditorContext()$path))
exp = read.csv('exp/GA_verify_exp.csv')
lognorm = read.csv('lognorm/GA_verify_lognorm.csv')

# draw qq-plot to see how the 100 quantiles of two distributions match each other

exp.V.q = quantile(exp$V1,probs = seq(0.01,1,0.01))
exp.E.q = quantile(exp$V2,probs = seq(0.01,1,0.01))

lognorm.V.q = quantile(lognorm$V1,probs = seq(0.01,1,0.01))
lognorm.E.q = quantile(lognorm$V2,probs = seq(0.01,1,0.01))

png('Figure B.1.png',width=778,height= 452)
par(mfrow=c(1,2))
qqplot(exp.V.q,exp.E.q,ylab = 'error process',xlab = 'Gaussian process',main='Exponential')
abline(0,1)
qqplot(lognorm.V.q,lognorm.E.q,ylab = 'error process',xlab = 'Gaussian process',main='Log-normal')
abline(0,1)
dev.off()
