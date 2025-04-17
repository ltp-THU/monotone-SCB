################# reproduce Table 2 ########################

summary_table2_hrs = function(model,n,p){
  df_names = df_names1 = df_names2 = df_names3 = rep(NA,10)
  if(n==1000 & p==27){
    for(i in 1:10) {
      df_names1[i]=paste0(  getwd(),'/../',model,n,'p',p,'sigma1/hr0.05-0.1/',i,'batch2000boot',model,n,'p',p,'sigma1.csv')
      df_names2[i]=paste0(  getwd(),'/../',model,n,'p',p,'sigma1/hr0.15-0.2/',i,'batch2000boot',model,n,'p',p,'sigma1.csv')
      df_names3[i]=paste0(  getwd(),'/../',model,n,'p',p,'sigma1/hr0.25-0.35/',i,'batch2000boot',model,n,'p',p,'sigma1.csv')
    }
    df1 = read.csv( df_names1[1],na.strings = "#NAME?")
    df2 = read.csv( df_names2[1],na.strings = "#NAME?")
    df3 = read.csv( df_names3[1],na.strings = "#NAME?")
    for (j in 2:10) {
      df1 = df1 + read.csv(df_names1[j],na.strings = "#NAME?")
      df2 = df2 + read.csv(df_names2[j],na.strings = "#NAME?")
      df3 = df3 + read.csv(df_names3[j],na.strings = "#NAME?")
    }
    df= rbind(df1/10,df2/10,df3/10)
    return(df[c('hr','is.cover90.hd','is.cover95.hd','SCB90.hd','SCB95.hd')])
  }else{
    for(i in 1:10) df_names[i]=paste0(  getwd(),'/../',model,n,'p',p,'sigma1/',i,'batch2000boot',model,n,'p',p,'sigma1.csv')
    df = read.csv( df_names[1],na.strings = "#NAME?")
    for (j in 2:10) df = df + read.csv(df_names[j],na.strings = "#NAME?")
    df = df/10
    
    #filtered_df <- df[abs(df$hr-0.075)>0.01,]
    return(df[c('hr','is.cover90.hd','is.cover95.hd','SCB90.hd','SCB95.hd')])
  }
}

summary_table2_GCV = function(model,n,p){
  df_names = rep(NA,10)
  if(model=='a') type = 'LS'
  if(model=='b') type = 'PLS'
  
  if(n==300 & p==9){
    df_names[1]=paste0(  getwd(),'/../GCV/',model,n,'p',p,'gcv/',1,'GCV2000boot',model,n,'p',p,'sigma1.csv')
    df = read.csv( df_names[1],na.strings = "#NAME?")
    return(apply(df[c('is.cover90','is.cover95','SCB90','SCB95')], 2, mean))
  }else{
    for(i in 1:10) df_names[i]=paste0(  getwd(),'/../GCV/',model,n,'p',p,'gcv/',i,'GCV2000boot',model,n,'p',p,'sigma1.csv')
    df = read.csv( df_names[1],na.strings = "#NAME?")
    for (j in 2:10) df = rbind(df, read.csv(df_names[j],na.strings = "#NAME?"))
    return(apply(df[c('is.cover90','is.cover95','SCB90','SCB95')], 2, mean))
  }
}


hrs = c('GCV',seq(0.05,0.35,0.05))

table2_a = table2_b = NA
for (n in c(300,500,1000)) {
  table2_a_n = table2_b_n = data.frame(hr=hrs)
  for (p in c(9,18,27)) {
    df_a = rbind(summary_table2_GCV('a',n,p)[c('is.cover90','is.cover95')], summary_table2_hrs('a',n,p)[c('is.cover90.hd','is.cover95.hd')])
    table2_a_n = cbind(table2_a_n,df_a)
    df_b = rbind(summary_table2_GCV('b',n,p)[c('is.cover90','is.cover95')], summary_table2_hrs('b',n,p)[c('is.cover90.hd','is.cover95.hd')])
    table2_b_n = cbind(table2_b_n,df_b)
  }
  table2_a = rbind(table2_a,table2_a_n)
  table2_b = rbind(table2_b,table2_b_n)
}

cbind(table2_a[-1,],table2_b[-1,-1])

write.csv(cbind(table2_a[-1,],table2_b[-1,-1]),file='table2.csv')


#################### reproduce Table B.2 ###########################
tableB1_a = tableB1_b = NA
for (n in c(300,500,1000)) {
  tableB1_a_n = tableB1_b_n = data.frame(hr=hrs)
  for (p in c(9,18,27)) {
    df_a = rbind(summary_table2_GCV('a',n,p)[c('SCB90','SCB95')], summary_table2_hrs('a',n,p)[c('SCB90.hd','SCB95.hd')])
    tableB1_a_n = cbind(tableB1_a_n,df_a)
    df_b = rbind(summary_table2_GCV('b',n,p)[c('SCB90','SCB95')], summary_table2_hrs('b',n,p)[c('SCB90.hd','SCB95.hd')])
    tableB1_b_n = cbind(tableB1_b_n,df_b)
  }
  tableB1_a = rbind(tableB1_a,tableB1_a_n)
  tableB1_b = rbind(tableB1_b,tableB1_b_n)
}

cbind(tableB1_a[-1,],tableB1_b[-1,-1])

write.csv(cbind(tableB1_a[-1,],tableB1_b[-1,-1]),file='tableB.2.csv')

##################### reproduce Figure B.2 ########################
library(patchwork)
library(ggplot2)
library(tidyr)

summary_hist = function(model,n,p){
  df_names = rep(NA,10)
  if(model=='a') type = 'LS'
  if(model=='b') type = 'PLS'
  if(n==300 & p==9){
    df_names[1]=paste0(  getwd(),'/../GCV/',model,n,'p',p,'gcv/',1,'GCV2000boot',model,n,'p',p,'sigma1.csv')
    df = read.csv( df_names[1],na.strings = "#NAME?")
    return(df)
  }else{
    for(i in 1:10) df_names[i]=paste0(  getwd(),'/../GCV/',model,n,'p',p,'gcv/',i,'GCV2000boot',model,n,'p',p,'sigma1.csv')
    df = read.csv( df_names[1],na.strings = "#NAME?")
    for (j in 2:10) df = rbind(df, read.csv(df_names[j],na.strings = "#NAME?"))
    return(df)
  }
}

df_LS <- summary_hist('a',1000,27)
df_PLS <- summary_hist('b',1000,27)

df_long_LS <- gather(df_LS[,6:8], key = "variable", value = "value")
df_long_PLS <- gather(df_PLS[,6:8], key = "variable", value = "value")

P1 <- ggplot(df_long_LS, aes(x = value, fill = variable, color = variable)) +
  geom_histogram(aes(y = ..density..),alpha=0.5, position = "identity") +labs(title = "model (a), n=1000, p=27", x = "bandwidth", y = "density")+
  theme_minimal() + xlim(0,0.35)+theme(legend.position = "none")
#+geom_vline(xintercept =  c(0.05,500^(-0.2)),show.legend = T,color=c('red','blue'),linetype='dashed')+theme(legend.position = "none")

P2 <- ggplot(df_long_PLS, aes(x = value, fill = variable, color = variable)) +
  geom_histogram(aes(y = ..density..),alpha=0.5, position = "identity") +labs(title = "model (b), n=1000, p=27", x = "bandwidth", y = "density")+
  theme_minimal() + xlim(0,0.35)#+geom_vline(xintercept =  c(0.05,500^(-0.2)),show.legend = T,color=c('red','blue'),linetype='dashed')

png('Figure B.2.png',width=963, height=362)
P1+P2
dev.off()





















