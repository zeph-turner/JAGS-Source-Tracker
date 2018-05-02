#Load required libraries
setup <- function() {
  require(rjags)
  library(reshape2)
  library(mosaic)
  library(dplyr)
}


#Returns otu_table, map, trueSources
#Note: draw from unknown not working yet - trueSources not correctly updated
#with column for draws from unknown.
create_data <- function(nSink, nSource, nOTU=20, drawFromKnown = TRUE, 
                        sparsity = 0.1, moreThanOne=TRUE, highOTU=500,
                        lowOTU=400){

  map <- data.frame(SampleID=paste("s", 0:(nSource+nSink-1), sep=""),
                    SourceSink=c(rep("source", nSource), rep("sink", nSink)),
                    Env = rep(NA, nSource+nSink))
  #Don't let it coerce to factor
  map$SampleID <- as.character(map$SampleID)
  
  otu_table <- as.data.frame(matrix(0, nrow=nOTU, ncol=(nSource+nSink)))
  colnames(otu_table) <- map$SampleID
  rownames(otu_table) <- paste("OTU", rownames(otu_table), sep="")
  
  #Seed sources
  for(i in 1:nSource){
    if(!grepl("unk", colnames(otu_table)[i])){
      otu_table[i,i] = highOTU;          #Set one OTU to highOTU
      #Draw others from uniform multinomial. Adjusting number drawn changes J-S divergence.
      otu_table[-(i),i] = rmultinom(1, lowOTU*(nrow(otu_table)-1), rep(1, nrow(otu_table)-1))
      
      #Sparisfy
      for(j in 1:nrow(otu_table)){
        if(runif(1) < sparsity) {        #Adjust sparsity by changing this condition
          otu_table[j, i] <- 0
        }
      }
    }
  }
  
  #Seed sinks by drawing from sources and record from which sources they were drawn
  trueSources = matrix(rep(0, nSink*nSource), nrow=nSink, ncol=nSource)
  if(drawFromKnown){
    for(i in (nSource+1):(nSource+nSink)){
      k <- resample(1:nSource, 1)                               #Pick source
      otu_table[,i] <- rmultinom(1, 200, unlist(otu_table[,k])) #Draw from source into sample
      trueSources[(i-nSource), k] <- trueSources[(i-nSource), k] + 1                #Record draw
      j = nSource-1
      print(paste("Sample ", i, " drawn from ", k))
      while(resample(1:2, 1)>1 && j > 0 && moreThanOne){ 
        #Add some number of additional sources
        k = resample(1:(nSource), 1)
        otu_table[,i] = otu_table[,i] + rmultinom(1, 200, unlist(otu_table[,k]))
        j<-j-1
        trueSources[(i-nSource), k] <- trueSources[(i-nSource), k] + 1
        print(paste("Sample ", i, " drawn from ", k))
      }
    }
  } else {
    for(i in (nSource+1):(nSource+nSink)){
      otu_table[,i] <- rmultinom(1, 200, runif(nrow(otu_table), 1, 50))
      trueSources[(i-nSource),nSource] <- 1
    }
  }
  
  #Normalize true sources matrix for later model evaluation
  for(i in 1:nrow(trueSources)){
    trueSources[i,] <-  trueSources[i,]/sqrt(sum(trueSources[i,])^2)
  }
  
  returnObj <- list("otu_table"=otu_table, "map"=map, "trueSources"=trueSources)
  return(returnObj)
}

#Returns void
#If run right after create_data, unkn_count should be 0 (we haven't 
#added unknown columns to the otu table yet). If after estimate_proportions,
#unkn_count should be the same count used in estimate_proportions.
write_otu_table <- function(otu_table, nSink, nSource, filename, unkn_count=0) {
  #Write OTU table without unknown rows
  write.table(otu_table[,c(1:(nSource-unkn_count), (nSource+1):(nSource+nSink))],
              file=filename)
}

#Returns codaSamples
#Inputs: Map with *no* unknown row, column SourceSink with entries "source" or "sink",
#and column SampleID matching OTU table.
#OTU table with OTUs as rows and samples as columns.
#nSource = number of source samples (NOT including unknown)
#nSink = number of sink samples
#unkn_count = number of unknown categorical distributions to generate
#chains = number of restarts of the Gibbs sampler to run
#adapt = number of adaptation iterations to fine-tune the Gibbs sampler
#burnin = number of full cycles through all parameter components to run before sampling posterior distibution
#samplesperchain = number of posterior parameter samples to take per chain
#alpha = prior count of sources in each sink
#beta = prior count of taxa in each source
estimate_proportions <- function(otu_table, map, nSink, nSource, unkn_count=0,
                                 chains=4, adapt=500, burnin=500, samplesPerChain=500,
                                 alpha=0.1, beta=10){
  #Data manipulation to *reduce* number of errors
  otu_table <- as.data.frame(otu_table)
  map <- select(map, SampleID, SourceSink)
  
  #Create unknown columns and rows
  if(unkn_count > 0){
    for(i in 1:unkn_count){
      #Create new row in map for unkn_count
      map[nrow(map)+1,] <- c(paste("unk", i, sep=""), "source")
      
      #Add column to OTU table and create column name
      otu_table[,ncol(otu_table)+1] <- rep(0, nrow(otu_table))
      colnames(otu_table)[ncol(otu_table)] <- paste("unk", i, sep="")
    }
  }
  
  #Save various counts that will be useful later
  #nSource is now UPDATED to INCLUDE unknown sources!
  nOtu = nrow(otu_table)
  nSource = nrow(map[map$SourceSink=="source",])
  nSamp = nrow(map[map$SourceSink=="sink",])
  nTax=nrow(otu_table)
  
  #Initialize vector to store output samples of alpha
  codaSamples=c()
  
  #Get source and sink column names from map
  sourceNames=as.vector(map[map$SourceSink=="source","SampleID"]) 
  sourceNames <- unlist(sourceNames)
  #print(sourceNames)
  sinkNames = as.vector(map[map$SourceSink=="sink","SampleID"])
  sinkNames <- unlist(sinkNames)
  
  #Rearrange OTU table so that otu_table[,1:nSource] is sources and
  # otu_table[,nSource+1:nSource+nSink] is sinks.
  otu_table <- otu_table[,c(sourceNames, sinkNames)]
  
  print("Source guide: ")
  for(i in 1:length(sourceNames)){
    print(paste0(i, ": ", sourceNames[i]))
  }
  
  #Specify Bayesian model (same model for each sink)
  modelString=paste("
  model {
  for(j in 1:nSource) {
    A[j, 1:nTax] ~ ddirch(rep(", alpha, ", nTax))
  }
  for(k in 1:nAllSourceObservations) {
    observation[k] ~ dcat(A[sampleNumber[k], 1:nTax])
  }
  for(i in 1:nSeq) {
    z[i] ~ dcat(B)
    x[i] ~ dcat(A[z[i], 1:nTax])
  }
  B ~ ddirch(rep(", beta, ", nSource))
  }
  ")
  #Initialize lists of observations for each sample
  otu_lists <- list()
  
  #Expand OTU table into lists of observations and store in otu_lists
  #If sample s1 had 1 observation of OTU 1, 2 of OTU 2, and 3 of OTU 3, then
  #otu_lists[[s1]] will be (1, 2, 2, 3, 3, 3).
  for(i in 1:ncol(otu_table)) {
    vect = c()
    for(j in 1:nrow(otu_table)) {
      vect= c(vect, rep(j, otu_table[j,i]))
    }
    otu_lists[[i]] <- vect
    #otu_lists[i] is the list corresponding to 
    #otu_table[,i]; it contains both x and z
  }
  
  #Melt lists
  melted <- melt(otu_lists)
  
  #Keep only sources, because we'll pass the sink list separately.
  melted <- filter(melted, L1 %in% 1:nSource)
  
  #Split melted dataset into separate columns because we can't do this inside JAGS
  sampleNumber <- melted$L1
  observation <- melted$value
  
  #Count number of total observed sequences among all sources
  nAllSourceObservations <- length(sampleNumber)
  
  for(j in 1:nSamp){
    #Current source predictions (start random)
    
    #Specify which sample to guess source
    #i = sample ID
    thisSeq = nSource+j
    print(paste("Calculating for ", thisSeq))
    #Put sequences from sink sample into x
    x = otu_lists[[thisSeq]]
    #Number of sequences to assign source
    nSeq = length(otu_lists[[thisSeq]])
    #Initial source guesses
    z = rep("", nSeq)
    
    dataList = list(    # Put the information into a list.
      nSeq=nSeq, #Number of observed(?) sequences
      x=x,       #vector of sample observations
      nTax=nTax,       #number of taxa
      nSource=nSource, #number of sources incl. unknown
      sampleNumber=sampleNumber, #which sample we're on
      observation=observation,   #observed source data
      nAllSourceObservations = nAllSourceObservations #num of observed sequences from sources
    )
    
    writeLines( modelString , con="TEMPmodel.txt" )
    
    # Run the chains:
    jagsModel = jags.model( file="TEMPmodel.txt" , data=dataList, 
                            n.chains=chains , n.adapt=adapt )
    update( jagsModel , n.iter=burnin )
    print(paste("Saving ", thisSeq, " in ", (thisSeq-nSource)))
    codaSamples[thisSeq - nSource] = coda.samples( jagsModel , variable.names=c("B") ,
                                                   n.iter=samplesPerChain )
  }
  
  return(codaSamples)
}

calculate_r2 <- function(trueSources, sampleAlphas, nSink, nSource, unkn_count){
  r2 <- vector(mode="numeric", length=nSink)
  trueSources <- cbind(trueSources, matrix(0, nrow=nrow(trueSources), ncol=unkn_count))
  for(i in 1:nSink) {
    colMeansMat <- matrix(colMeans(matrix(sampleAlphas[[i]], ncol=nSource+unkn_count)))
    trueSourcesMat <- matrix(trueSources[i,])
    ssRes <- sum((trueSourcesMat-colMeansMat)^2)
    ssTot <- sum((trueSourcesMat-mean(trueSourcesMat))^2)
    r2[i] <- 1-(ssRes/ssTot)
  }
  return(r2)
}

print_violinplots <- function(nSink, nSource, nUnknown, codaSamples) {
  for(i in 1:nSink){
    samplepost <- matrix(codaSamples[[i]], ncol=(nSource+nUnknown))
    samplepost.df <- data.frame(samplepost)
    print(ggplot(data=melt(samplepost.df, measure.vars=colnames(samplepost.df))) + 
            geom_violin(aes(y=value, x=variable), scale="width") +
            ggtitle(paste("Estimated source proportions for sample ", i)) +
            xlab("Source") + ylab("Proportion"))
  }
}

print_boxplots <- function(nSink, nSource, nUnknown, codaSamples) {
  for(i in 1:nSink){
    samplepost <- matrix(codaSamples[[i]], ncol=(nSource+nUnknown))
    samplepost.df <- data.frame(samplepost)
    print(ggplot(data=melt(samplepost.df, measure.vars=colnames(samplepost.df))) + 
            geom_boxplot(aes(y=value, x=variable)) +
            ggtitle(paste("Estimated source proportions for sample ", i)) +
            xlab("Source") + ylab("Proportion"))
  }
}
