### This script will tidy Dan's simulated sparse datasets into monolix friendly datasets
#sink(paste0("combLog", make.names(Sys.Date()), ".txt"), type = c("output", "message"))
##First read in the arguments listed at the command line
args=commandArgs(trailingOnly = TRUE)

if(length(args)==0){
    stop("No inputs supplied.")
}

sparse_input = args[1]
if(length(args) > 1) output_fname = args[2] else output_fname = "sparsedata_mono.csv"
if(length(args) > 2) truetime_input = args[3] else truetime_input = NULL

cat("Processing sim data on: ", paste(Sys.Date()), "\n")
cat("Input simulation data: ", paste(sparse_input), "\n")
if(is.null(truetime_input)) true_str = "non given" else true_str = truetime_input
cat("True time data: ", paste(true_str), "\n")
cat("Output filename: ", paste(output_fname), "\n")

if(file.exists(output_fname)) stop("won't overwrite")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  suppressWarnings(library(stringr))
})

tidy_sim = function(sparse_fname, truetime_fname = NULL){
  sparse_data_in = read.csv(sparse_fname, stringsAsFactors = F) 

  sparse_data =  sparse_data_in %>%
    gather(weeks, log10_VL, -X) %>%
    rename(simID = X) %>%
    arrange(simID) %>%
    mutate(
      simID = simID + 1,
      weeks = str_remove_all(weeks, "X"))
  
  # error check
  for(i in 1:nrow(sparse_data_in)){
    if(!all(subset(sparse_data, simID == i)$log10_VL == sparse_data_in[i,-1])) stop("transpose error")
  }
  
  if(!is.null(truetime_fname)){
     time_dat = read.csv(truetime_fname, stringsAsFactors = F) %>%
      rename(simID = X, truetime = X0) %>%
      mutate(simID = simID + 1)
      sparse_data_out = left_join(sparse_data, time_dat, by = "simID")
  } else sparse_data_out = sparse_data
  if(nrow(sparse_data) != nrow(sparse_data_out)) stop("merging error with true times")
  
  sparse_data_out
}

write.csv( tidy_sim(sparse_input, truetime_input), output_fname, row.names = F)

