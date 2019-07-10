rm(list = ls())

data <- read.csv(paste0("~/Documents/GitHub/AMP-timing/Monolix/Data/DBRout_RV217.csv"), header=TRUE,stringsAsFactors=FALSE)
datacopy = data.frame(data)

IDs = unique(data[,"ID"])
colnames(datacopy)[7]="y"

datacopy = cbind(datacopy,rep(1,dim(datacopy)[1]))
colnames(datacopy)[8]="ytype"

datacopy = cbind(datacopy,rep("VL",dim(datacopy)[1]))
colnames(datacopy)[9]="ytypeID"

datacopy = cbind(datacopy,rep(0,dim(datacopy)[1]))
colnames(datacopy)[10]="Cens"


datacopy$Cens[which(datacopy$VL<=30)]=1

datanew=c()
for(ID in IDs){
  
  data_ind = datacopy[which(ID==datacopy$ID),]
  
  if (data_ind$VL[1]==1){
    data_ind = data_ind[which(data_ind$days>0),]
  }
  
  data_ind$days=data_ind$days-min(data_ind$days)
  data_ind$y = log10(data_ind$VL)
  
  # obtain CD4s and compute log10CD4s
  data_ind_CD4s = data_ind[which(!is.na(data_ind$CD4)),]
  data_ind_CD4s$y = log10(data_ind_CD4s$CD4)
  data_ind_CD4s$ytype = rep(2, length(data_ind_CD4s$ytype))
  data_ind_CD4s$ytypeID = rep("CD4s", length(data_ind_CD4s$ytype))
  
  # obtain CD8s and compute log10CD8s
  data_ind_CD8s = data_ind[which(!is.na(data_ind$CD8)),]
  data_ind_CD8s$y = log10(data_ind_CD8s$CD8)
  data_ind_CD8s$ytype = rep(3, length(data_ind_CD8s$ytype))
  data_ind_CD8s$ytypeID = rep("CD8s", length(data_ind_CD8s$ytype))
 
   # add to individual rows
  datanew = rbind(datanew,data_ind,data_ind_CD4s,data_ind_CD8s)
  
}

write.csv(datanew,
          file=paste0("~/Documents/GitHub/AMP-timing/Monolix/Data/RV217_VL_Tcells_Monolix.csv"), 
          row.names = FALSE)

