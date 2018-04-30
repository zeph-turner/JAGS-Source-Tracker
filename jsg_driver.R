#Main method

#Source jsg_functions
setwd("~/GitHub/JAGS-Source-Tracker")
fileroot <- getwd()
source(paste0(fileroot, "/jsg_functions.R"))

#Run setup functions
setup()
fileNameRoot <- paste0(fileroot, "/")

#Loop used for simulation studies of JAGS ST. 
r2=list()
nRuns <- 1
for(i in 1:nRuns){
  print(paste0("Run ", i, " started: "))
  this <- create_data(1, 5, 20, TRUE, 0.2, FALSE, 50, 45)
  write_otu_table(this$otu_table, 1, 5, paste0(fileNameRoot, "table", i, ".csv"))
  t <- this$trueSources
  dump(c("t"), file=paste0(fileNameRoot, "trueSources", i, ".R"))
  samples <- estimate_proportions(this$otu_table, this$map, 1, 5, 1, chains=4, 
                                  adapt=500, burnin=500, samplesPerChain=500,
                                  alpha=10, beta=10)
  r2[[i]] <- calculate_r2(this$trueSources, samples, 1, 5, 1)
  dump(c("samples", "this"), file=paste0(fileNameRoot, "allData", i, ".R"))
  print(paste0("Run ", i, " concluded."))
}

dump(c("r2"), file=paste(fileNameRoot, "r2summary.R", sep=""))


this <- create_data(1, 5, 20, TRUE, 0.2, FALSE, 500, 450)
t <- this$trueSources
samples <- estimate_proportions(this$otu_table, this$map, 5, 1, 1, chains=4, 
                                adapt=500, burnin=500, samplesPerChain=500,
                                alpha=1.0, beta=.1)
r2 <- calculate_r2(this$trueSources, samples, 1, 5, 1)
print_violinplots(1, 5, 1, samples)
