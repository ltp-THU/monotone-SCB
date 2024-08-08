process = function(ts){
  # missing check
  ts = ts[!is.na(ts$yyyy),]
  years = seq(min(ts$yyyy), min( c(max(ts$yyyy),2024) ),by=1 )
  df = matrix(NA,nrow = 12*length(years),ncol = 3)
  for (i in 1:length(years)){
    df[((i-1)*12+1):(i*12) , 1] = years[i]
    df[ ((i-1)*12+1):(i*12) , 2] = seq(1,12,1)
  }
  for (i in 1:nrow(df)){
    if(length(ts[ (ts$yyyy==df[i,1] & ts$mm==df[i,2]) ,3])==1){
      df[i,3] = ts[ (ts$yyyy==df[i,1] & ts$mm==df[i,2]) ,3]
    }
  } 
  
  # interpolation
  df_inter = data.frame(year=df[,1],month=df[,2],tmax=df[,3])
  for (m in which(is.na(df[,3])) ) {
    month_m = df[m,2]
    year_m = df[m,1]
    ts_m = ts(df[df[,2]==month_m,3],start=years[1],end=years[length(years)])
    df_inter[m,3] = approx(x=time(ts_m),y=ts_m,xout = year_m)$y
  }
  
  # eliminate seasonal effect
  # 若第一年不满，去除
  while ( sum(!is.na(df_inter[df_inter$year==years[1],3]))<12 ) {
    df_inter = df_inter[df_inter$year>years[1],]
    years = years[-1]
  }
  # 若最后一年不满，去除
  while ( sum(!is.na(df_inter[df_inter$year==years[length(years)],3]))<12 ) {
    df_inter = df_inter[df_inter$year<years[length(years)],]
    years = years[-length(years)]
  }
  
  ts_data = ts(df_inter$tmax,start=c(years[1],1),end=c(years[length(years)],12),frequency=12)
  
  stl_data = stl(ts_data,s.window = "periodic")
  seasonal_adj_data = ts_data - stl_data$time.series[,1]
  return( list(ts_data=ts_data,stl_data=stl_data$time.series,seasonal_adj=seasonal_adj_data) )
}

#### deseasonalize and interpolate missing ####

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
files = dir('../raw_data')
# not include locations.csv
files = files[files!='locations.csv']

ts_all = list()

for (i in 1:length(files)) {
  raw_data = read.table(paste0('../raw_data/',files[i]),skip = 7,col.names=c('yyyy', 'mm', 'tmax','tmin','af','rain','sun')
                        ,fill = T)
  raw_data = apply(raw_data, 2, as.numeric)
  
  ts = raw_data[,1:3]
  ts = as.data.frame(ts)
  
  stl_ts = process(ts)
  
  ts_all[[i]] = stl_ts$seasonal_adj
}

starts = c()
ends = c()
for (i in 1:length(files)) {
  starts = c(starts,start(ts_all[[i]])[1] )
  ends = c(ends,end(ts_all[[i]])[1] )
}

# choose time series ranged from 1979 to 2022
start_year = 1979
end_year = 2023

ts_all.trunc = ts_all[[1]][time(ts_all[[1]]) < end_year+1 & time(ts_all[[1]]) >= start_year] #time() compute yyyy/m as yyyy+(m-1)/12

for (i in which(ends >= end_year & starts<=start_year) ) {
  ts_all.trunc = cbind(ts_all.trunc,  ts_all[[i]][time(ts_all[[i]]) <= end_year & time(ts_all[[i]]) >= start_year])
}

ts_all.trunc = data.frame(ts_all.trunc[,-1])


names = files[which(ends>=end_year & starts<=start_year)]
for (i in 1:length(names)) {
  names[i] = substr(names[i],start = 1,stop = nchar(names[i])-8)
}
colnames(ts_all.trunc) = names


write.csv(ts_all.trunc,file = '../estimation/data_processed.csv')

