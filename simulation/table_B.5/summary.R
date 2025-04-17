summary_table = function(model,penal){
  df_names = rep(NA,10)
  for (i in 1:10) {
    df_names[i]=paste0(  getwd(),'/',model,500,'p',27,'penal',penal,'/',i,'GCV2000boot',model,500,'p',27,'sigma1.csv')
  }
  df = read.csv( df_names[1],na.strings = "#NAME?")
  for (j in 2:10) df = rbind(df, read.csv(df_names[j],na.strings = "#NAME?"))
  return(apply(df[c('is.cover90','is.cover95','SCB90','SCB95','shape')], 2, mean))
}

penals = c(0.3,0.2,0.1,0)
models = c('a','b')
dfa = matrix(NA,4,5)
dfb = matrix(NA,4,5)
for (i in 1:4) {
  dfa[i,] = summary_table('a',penals[i])
  dfb[i,] = summary_table('b',penals[i])
}

table_B5 = cbind(dfa,dfb)

# table B.5
write.csv(table_B5,file = 'tableB.5.csv')



