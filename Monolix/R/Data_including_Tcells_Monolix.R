rm(list = ls())

data <- read.csv(paste0("~/Documents/GitHub/AMP-timing/Monolix/Data/DBRout_RV217.csv"), header=TRUE,stringsAsFactors=FALSE)
datacopy = data.frame(data)

IDs = unique(data[,"ID"])
colnames(datacopy)[7]="y"

datacopy = cbind(datacopy,rep(1,dim(datacopy)[1]))
colnames(datacopy)[8]="ytype"

datanew=c()
for(ID in IDs){
  
  data_ind = datacopy[which(ID==datacopy$ID),]

  # deleting day 0
  if (data_ind$VL[which(data_ind$days==0)]>1){
    data_ind_V = data_ind
  }
  else
  {
    data_ind_V = data_ind[which(data_ind$days>0),]
  }
  data_ind_V$y = log10(data_ind_V$VL)
  
  # obtain CD4s and compute log10CD4s
  data_ind_CD4s = data_ind[which(!is.na(data_ind$CD4)),]
  data_ind_CD4s$y = log10(data_ind_CD4s$CD4)
  data_ind_CD4s$ytype = rep(2, length(data_ind_CD4s$ytype))
  
  # obtain CD8s and compute log10CD8s
  data_ind_CD8s = data_ind[which(!is.na(data_ind$CD8)),]
  data_ind_CD8s$y = log10(data_ind_CD8s$CD8)
  data_ind_CD8s$ytype = rep(3, length(data_ind_CD8s$ytype))
  
  # add to individual rows

  datanew = rbind(datanew,data_ind_V,data_ind_CD4s,data_ind_CD8s)
  
}

write.csv(datanew,
          file=paste0("~/Documents/GitHub/AMP-timing/Monolix/Data/RV217_VL_Tcells_Monolix.csv"), 
          row.names = FALSE)

