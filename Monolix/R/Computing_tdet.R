rm(list = ls())

# -Fabian Cardozo
rm(list = ls())

### MODELS
HIV_acute_model <- function(t,x,params){
  with(as.list(x),{   
    # Set current state values
    
    S  = x[1]
    AU  = x[2]
    AP  = x[3]
    P  = x[4]
    E  = x[5]
    V  = x[6]
    
    
    ddt_S = aS - dS*S - Bt0*S*V                      
    ddt_AU = (1-tau)*Bt0*S*V - dI*AU - k*E*AU 
    ddt_AP = tau*Bt0*S*V - dI*AP - k*E*AP 
    
    ddt_P = aE + w*(1-f)*P*(AP+AU)/(1+(AP+AU)/I50) -dp*P 
    ddt_E = w*f*P*(AP+AU)/(1+(AP+AU)/I50) - dE*E
    
    ddt_V = p*AP - g*V - Bt0*S*V  
    
    der <- c(ddt_S,ddt_AU,ddt_AP,ddt_P,ddt_E,ddt_V)
    #print(x)
    list(der)
  })       
}


###########################
# MAIN PROGRAM
###########################

#load R library for ordinary differential equation solvers
library(deSolve)
library(gplots)
require(lattice)      
require(latticeExtra) 
library(plotrix)

Project = "/Model_4/1_Project_model_4_v2_first_neg"

predictions <- read.table(paste0("~/Documents/GitHub/AMP-timing/Monolix/Projects",Project,"/ChartsData/IndividualFits/log10VL_fits.txt"), sep=",",header=TRUE,stringsAsFactors=FALSE)
estimates <- read.table(paste0("~/Documents/GitHub/AMP-timing/Monolix/Projects",Project,"/IndividualParameters/estimatedIndividualParameters.txt"), header=TRUE,stringsAsFactors=FALSE,sep=",")
data <- read.csv(paste0("~/Documents/GitHub/AMP-timing/Monolix/Data/wpd_all_10wks_1st_neg.csv"), header=TRUE,stringsAsFactors=FALSE)

IDs = unique(predictions[,"id_var"])


tdetp=c()
tdet=c()
for(ID in IDs){
  
  Vdata = 10^as.numeric(data[which(ID==data[,"id_var"]),"log10VL"])
  tVdata = as.numeric(data[which(ID==data[,"id_var"]),"days"])
  
  
  Vdatap = 10^as.numeric( predictions[which(ID==predictions[,"id_var"]),"indivPredMode"])
  tVdatap = as.numeric( predictions[which(ID==predictions[,"id_var"]),"time"])

  tdetp = c(tdetp,tVdatap[which.min(abs(Vdatap-30))])
  
  
  index_e = which(ID==IDs)
  
  initT =estimates[index_e,"initT_mode"]
  aS =estimates[index_e,"aS_mode"]
  dS =estimates[index_e,"dS_mode"]
  tau =estimates[index_e,"tau_mode"]
  lBt0 =estimates[index_e,"lBt0_mode"]
  lp =estimates[index_e,"lp_mode"]
  dI =estimates[index_e,"dI_mode"]
  lk =estimates[index_e,"lk_mode"]
  f =estimates[index_e,"f_mode"]
  dp =estimates[index_e,"dp_mode"]
  lw =estimates[index_e,"lw_mode"]
  dE =estimates[index_e,"dE_mode"]
  lI50 =estimates[index_e,"lI50_mode"]
  
  
  lae = lk
  g = 23        
  
  Bt0 = 10^(lBt0) 
  p = 10^(lp)
  k = 10^(lk)        
  w = 10^(lw)     
  I50 = 10^(lI50) 
  aE = 10^(lae) 
 
  t0 = initT
  S_0 = aS/dS                             
  
  I_0 = 1e-6
  
  V_0 = p*I_0/g                   
  AU_0 = tau*I_0
  AP_0 = (1-tau)*I_0
  P_0 = aE/dp
  E_0 = 0      
  

  tend = 70
  init.x <- c(S_0,AU_0,AP_0,P_0,E_0,V_0)
  out <- as.data.frame(lsodar(init.x,seq(t0,tend,by=0.1),HIV_acute_model,params)) 
  
  t.out = out[,1]
  
  S  = out[,2]
  AU  = out[,3]
  AP  = out[,4]
  P  = out[,5]
  E  = out[,6]
  V  = out[,7]
  
  
  tdet = c(tdet,t.out[which.min(abs(V*1e3-30))])
  
}
names(tdet)=IDs
hist(tdet)

write.table(tdet,
          file=paste0("~/Documents/GitHub/AMP-timing/Monolix/Projects",Project,"/tdet_estimates.txt"), 
          row.names = FALSE)

