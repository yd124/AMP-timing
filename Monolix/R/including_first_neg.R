rm(list = ls())

data <- read.csv(paste0("~/Documents/GitHub/AMP-timing/Monolix/Data/wpd_all_10wks.csv"), header=TRUE,stringsAsFactors=FALSE)
datacopy = data

IDs = unique(data[,"id_var"])
datacopy = cbind(datacopy,rep(0,dim(datacopy)[1]))
colnames(datacopy)[7]="Cens"


datanew=c()
for(ID in IDs){
  
  data_ind = datacopy[which(ID==data[,"id_var"]),]

  row_temp = data_ind[which.min(as.numeric(data_ind[,"days"])),]
  row_temp["days"]=as.numeric(row_temp["days"])-3.5
  row_temp["log10VL"]=log10(15)
  row_temp["Cens"]=1
  
  data_ind = rbind(row_temp,data_ind)
  datanew = rbind(datanew,data_ind)

}

write.csv(datanew,
          file=paste0("~/Documents/GitHub/AMP-timing/Monolix/Data/wpd_all_10wks_1st_neg.csv"), 
          row.names = FALSE)

