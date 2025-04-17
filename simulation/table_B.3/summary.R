
summary_tableB3 = function(model,n,p){
  df_names = df_names1 = df_names2 = df_names3 = rep(NA,10)
  
  if(n==1000 & p==27){
    for(i in 1:10) {
      df_names1[i]=paste0(  model,n,'p',p,'zhou/hr0.05-0.1/',i,'batch2000boot', model,'40p',p,'sigma1zhou.csv')
      df_names2[i]=paste0(  model,n,'p',p,'zhou/hr0.15-0.2/',i,'batch2000boot', model,'40p',p,'sigma1zhou.csv')
      df_names3[i]=paste0(  model,n,'p',p,'zhou/hr0.25-0.35/',i,'batch2000boot',model,'40p',p,'sigma1zhou.csv')
      
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
    return(df[c('max_width90.zhou','max_width95.zhou')])
  }else{
    for(i in 1:10) df_names[i]=paste0(model,n,'p',p,'zhou/',i,'batch2000boot', model,'40p',p,'sigma1zhou.csv')
    df = read.csv( df_names[1],na.strings = "#NAME?")[,1:4]
    for (j in 2:10) df = df + read.csv(df_names[j],na.strings = "#NAME?")[,1:4]
    df = df/10
    #filtered_df <- df[abs(df$hr-0.075)>0.01,]
    return(df[c('max_width90.zhou','max_width95.zhou')])
  }
}

hrs = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35)
tableB3_a = tableB3_b = NA
for (n in c(300,500,1000)) {
  tableB3_a_n = tableB3_b_n = data.frame(hr=hrs)
  for (p in c(9,18,27)) {
    df_a = summary_tableB3('a',n,p)
    df_b = summary_tableB3('b',n,p)
    tableB3_a_n = cbind(tableB3_a_n,df_a)
    tableB3_b_n = cbind(tableB3_b_n,df_b)
  }
  tableB3_a = rbind(tableB3_a,tableB3_a_n)
  tableB3_b = rbind(tableB3_b,tableB3_b_n)
}

cbind(tableB3_a[-1,],tableB3_b[-1,-1])

write.csv(cbind(tableB3_a[-1,],tableB3_b[-1,-1]),file='tableB.3.csv')

