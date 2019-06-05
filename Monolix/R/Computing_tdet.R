rm(list = ls())


Project = "/Model_4/1_Project_model_4_v2_first_neg"

predictions <- read.table(paste0("~/Documents/GitHub/AMP-timing/Monolix/Projects",Project,"/predictions.txt"), sep=",",header=TRUE,stringsAsFactors=FALSE)
 
IDs = unique(predictions[,"id"])

tdet=c()
for(ID in IDs){
  
  Vdata = 10^as.numeric( predictions[which(ID==predictions[,"id"]),"log10VL"])
  tVdata = as.numeric( predictions[which(ID==predictions[,"id"]),"time"])

  tdet = c(tdet,tVdata[which.min(abs(Vdata-15))])
  
}
names(tdet)=IDs
hist(tdet)

write.table(tdet,
          file=paste0("~/Documents/GitHub/AMP-timing/Monolix/Projects",Project,"/tdet_estimates.txt"), 
          row.names = FALSE)

