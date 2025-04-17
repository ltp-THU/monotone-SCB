

summary_batch = function(model='exp'){
  df_names = rep(NA,10)
  for(i in 1:10){
    df_names[i]=paste0(model,'/',i,'GCV2000boot40LS27sigma1n500.csv')
  } 
  
  files = dir()
  df = read.csv( df_names[1],na.strings = "#NAME?")
  for (j in 2:10) {
    df = rbind(df, read.csv(df_names[j],na.strings = "#NAME?"))
  }
  return(apply(df[c('is.cover90','is.cover95','SCB90','SCB95')], 2, mean))
}

# table B.4
write.csv(rbind(summary_batch('exp'),summary_batch('lognorm')),file = 'tableB.4.csv')

