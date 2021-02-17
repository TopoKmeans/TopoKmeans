# Data generation #

# daily #
data <- covid19us::get_states_daily(state = "all", date = "2021-01-11")
data<-data[c(3,4,5,8,9,10,11,12,13,14,15,20,21,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,42,43,44,45,46,47,49,50,51,52,53)]
N<-dim(data)[2]

data[is.na(data)] = 0
cs<-colSums(data)
for (i in 1:N){
  if (cs[i] >0){
    data[,i]<-data[,i]/cs[i]
  }
}

write.csv(data, file = "covid19us_01_11.csv")

# 7-day window #
for (ii in c(10:16)) {
  tmp_date <- paste("2021","01",ii,sep = "-")
  data <- covid19us::get_states_daily(state = "all", date = tmp_date)
  
  data<-data[c(3,4,5,8,9,10,11,12,13,14,15,20,21,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,42,43,44,45,46,47,49,50,51,52,53)]
  N<-dim(data)[2]
  
  data[is.na(data)] = 0
  cs<-colSums(data)
  for (i in 1:N){
    if (cs[i] >0){
      data[,i]<-data[,i]/cs[i]
    }
  }
  if (ii == 10) {
    final_data <- data
  }else{
    final_data <- data+ final_data
  }
}
avg_final_data <- final_data/7
write.csv(avg_final_data, file = "avg_final_data_01_10_01_16.csv")
