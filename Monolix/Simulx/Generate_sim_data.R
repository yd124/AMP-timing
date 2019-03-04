library("mlxR")

setwd("~/Dropbox/RV217/AMP-timing/Monolix/Simulx")
project.file = "../Projects/Model_4/1_Project_model_4_v12_lod_15_cors.mlxtran"

N <- 1000


Vout  <- list(name = 'Vout',time = seq(0, 70, by=7))

outp <- c("initT", "aS","dS", "tau", "lBt0","lp", "dI", "lk", "f","dp", "lw", "dE", "lI50","lod")

res1  <- simulx(project = project.file,
                group   = list(size = N),
                output    = list(Vout,outp),
                fim     = "lin",
                result.file = "sim_data_lod_15.csv")