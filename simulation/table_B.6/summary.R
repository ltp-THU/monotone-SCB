
summary_tableB6_hrs = function(model,n){
  df_names = rep(NA,10)
  if(model=='c') type = 'LS'
  if(model=='d') type = 'PLS'
  
  for(i in 1:10) df_names[i]=paste0(model,n,'/',i,'batch2000boot40',type,'_simul_SCB_reg_estn',n,'.csv')
  df = read.csv( df_names[1],na.strings = "#NAME?")
  for (j in 2:10) df = df + read.csv(df_names[j],na.strings = "#NAME?")
  df = df/10
  return(df[c('coverage90','coverage95')])
}


summary_tableB6_GCV = function(model,n){
  df_names = rep(NA,10)
  if(model=='c') type = 'LS'
  if(model=='d') type = 'PLS'
  
  for(i in 1:10) df_names[i]=paste0('GCV/',model,n,'gcv/',i,'GCV2000boot40',type,'_simul_SCB_reg_estn',n,'.csv')
  df = read.csv( df_names[1],na.strings = "#NAME?")
  for (j in 2:10) df = rbind(df, read.csv(df_names[j],na.strings = "#NAME?"))
  return(apply(df[c('is.cover90','is.cover95')], 2, mean))
}




hrs = c('GCV',seq(0.05,0.35,0.05))

tableB6_c_n = tableB6_d_n = hrs
for (n in c(300,500,1000)) {
  df_c = rbind(summary_tableB6_GCV('c',n), summary_tableB6_hrs('c',n))
  tableB6_c_n = cbind(tableB6_c_n,df_c)
  df_d = rbind(summary_tableB6_GCV('d',n), summary_tableB6_hrs('d',n))
  tableB6_d_n = cbind(tableB6_d_n,df_d)
}


write.csv(cbind(tableB6_c_n,tableB6_d_n[,-1]),file='tableB.6.csv')



