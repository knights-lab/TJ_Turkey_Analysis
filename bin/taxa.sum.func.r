#Parses taxa_table to make proper taxa summary plots
#inputs are the taxa table and the subset you need for plotting


make_taxa_sums <- function(otu, union, cutoff){
  
  #Keep only samples from that subset
  subset_otu <- otu[,union, drop=FALSE]
  subset_otu <- as.data.frame(t(subset_otu))
  
  #if less than cutoff, set to 0 and change name to other
  for(i in 1:ncol(subset_otu)){
    for(k in 1:nrow(subset_otu)){
      if(subset_otu[k,i] < cutoff){
        subset_otu[k,i] <- 0
      }
    }
  }
  subset_otu$Other <- 0
  #Set other to be the counts subtracted earlier
  for(v in 1:nrow(subset_otu)){
    subset_otu[v,"Other"] <- (1 - rowSums(subset_otu[v,]))
  }
  
  #Keep only taxa present (remove 0 counts)
  subset_otu <- subset_otu[,colSums(subset_otu) > 0]
  
  #Find taxa present at highest numbers
  subset_otu_na <- subset_otu[,1:ncol(subset_otu)-1]
  subset_otu_na[subset_otu_na == 0] <- NA
  max_abund <- colMeans(subset_otu_na, na.rm=TRUE)
  names(max_abund) <- colnames(subset_otu_na)
  max_abund <- c(max_abund, 0)
  names(max_abund)[length(max_abund)] <- "Other"
  max_abund <- sort(max_abund, decreasing=TRUE)
  
  #add sample IDs to otu table and melt table
  subset_otu <- subset_otu[,names(max_abund)]
  subset_otu$SampleID <- rownames(subset_otu)
  
  otu_made <- melt(subset_otu, id.vars = "SampleID", variable.name = "Taxa", value.name = "Relative_Abundance")
  
  #Merge metadata
  new_map <- mapping[union,]
  otu_made <- merge(otu_made, new_map, by="SampleID")
  
  #Assign taxa levels (order in bar plot)
  otu_made$Taxa <- factor(otu_made$Taxa, levels = c(names(max_abund)), ordered=TRUE)
  
  otu_made <- otu_made[order(otu_made$Taxa),]
  return(otu_made)
}


### Same but for fungal!
make_taxa_sums2 <- function(otu, union, cutoff){
  
  #Keep only samples from that subset
  subset_otu <- otu[,union, drop=FALSE]
  subset_otu <- as.data.frame(t(subset_otu))
  
  #if less than cutoff, set to 0 and change name to other
  for(i in 1:ncol(subset_otu)){
    for(k in 1:nrow(subset_otu)){
      if(subset_otu[k,i] < cutoff){
        subset_otu[k,i] <- 0
      }
    }
  }
  subset_otu$Other <- 0
  #Set other to be the counts subtracted earlier
  for(v in 1:nrow(subset_otu)){
    subset_otu[v,"Other"] <- (1 - rowSums(subset_otu[v,]))
  }
  
  #Keep only taxa present (remove 0 counts)
  subset_otu <- subset_otu[,colSums(subset_otu) > 0]
  
  #Find taxa present at highest numbers
  subset_otu_na <- subset_otu[,1:ncol(subset_otu)-1]
  subset_otu_na[subset_otu_na == 0] <- NA
  max_abund <- colMeans(subset_otu_na, na.rm=TRUE)
  names(max_abund) <- colnames(subset_otu_na)
  if("Other" %in% colnames(subset_otu)){
    max_abund <- c(max_abund, 0)
    names(max_abund)[length(max_abund)] <- "Other"
  }
  max_abund <- sort(max_abund, decreasing=TRUE)
  #add sample IDs to otu table and melt table
  subset_otu <- subset_otu[,names(max_abund)]
  subset_otu$SampleID <- rownames(subset_otu)
  
  otu_made <- melt(subset_otu, id.vars = "SampleID", variable.name = "Taxa", value.name = "Relative_Abundance")
  
  #Merge metadata
  otu_made <- merge(otu_made, f_map, by="SampleID")
  
  #Assign taxa levels (order in bar plot)
  otu_made$Taxa <- factor(otu_made$Taxa, levels = c(names(max_abund)), ordered=TRUE)
  
  otu_made <- otu_made[order(otu_made$Taxa),]
  return(otu_made)
}
