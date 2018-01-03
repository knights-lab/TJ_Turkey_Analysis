#Heat map of transcripts/pathways and OTUs

#Make the names of the samples the same for transcripts and otu table
samples2 <- rownames(mapping[mapping$SampleID %in% samples & mapping$transcriptomic %in% colnames(ileum_transcript),])
samples3 <- mapping[samples2,"transcriptomic"]

x <- as.matrix(t(CLR_otutable))
rownames(x) <- taxonomy3$V2
Cor_out2 <- cor(t(x))
Cor_out3 <- abs(Cor_out2) > 0.95
new_taxatable <- x

for(i in 1:nrow(Cor_out3)){
  for(j in i:nrow(Cor_out3)){
    if(i==j){
      next
    }
    if(is.na(Cor_out3[j,i])){
      next
    }
    if(Cor_out3[j,i]){
      new_taxatable[j,] <- new_taxatable[j,] + new_taxatable[i,]
      new_taxatable[i,] <- rep(0, ncol(new_taxatable))
    }
  }
}

#Collapse
taxa_table2 <- new_taxatable[rowSums(new_taxatable) > 0, ]
x <- taxa_table2


x <- t(x[,samples2])
# replace rownames with unique ones
otu_names <- as.character(colnames(x))
taxa <- c()
for (i in 1:length(otu_names)){
  genus <- strsplit(otu_names[i], "; ", fixed=T)[[1]][6]
  family <- strsplit(otu_names[i], "; ", fixed=T)[[1]][5]
  order <- strsplit(otu_names[i], "; ", fixed=T)[[1]][4]
  class <- strsplit(otu_names[i], "; ", fixed=T)[[1]][3]
  phylum <- strsplit(otu_names[i], "; ", fixed=T)[[1]][2]
  if(family== "f__"){
    taxon <- gsub("o__", "", order)
    taxon <- paste("order:", taxon)
  } else {
    if(genus == "g__"){
      taxon <- gsub("f__", "", family)
      taxon <- paste("family:", taxon)
    } else {
      taxon <- gsub("g__", "", genus)
    }
  }
  taxa <- c(taxa, taxon)
}
taxa <- data.table(taxa)
taxa <- taxa[, Index := 1:.N , taxa ]
taxa$taxa <- paste(taxa$taxa, taxa$Index, sep="_")
taxa$taxa <- gsub("_1", "", taxa$taxa)
colnames(x) <- taxa$taxa


y <- t(ileum_transcript[,samples3])


#correlations <- associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
# Or, alternatively, the same output is also available in a handy table format
correlation.table <- associate(x, y, method = "pearson", mode = "table", p.adj.threshold = 0.05, n.signif = 1)

# replace rownames with unique ones
#otu_names <- as.character(rownames(x))
# otu_names <- as.character(correlation.table$X1)
# taxa <- c()
# for (i in 1:length(otu_names)){
#   genus <- strsplit(otu_names[i], "; ", fixed=T)[[1]][6]
#   family <- strsplit(otu_names[i], "; ", fixed=T)[[1]][5]
#   order <- strsplit(otu_names[i], "; ", fixed=T)[[1]][4]
#   class <- strsplit(otu_names[i], "; ", fixed=T)[[1]][3]
#   phylum <- strsplit(otu_names[i], "; ", fixed=T)[[1]][2]
#   if(family== "f__"){
#     taxon <- gsub("o__", "", order)
#     taxon <- paste("order:", taxon)
#   } else {
#     if(genus == "g__"){
#       taxon <- gsub("f__", "", family)
#       taxon <- paste("family:", taxon)
#     } else {
#       taxon <- gsub("g__", "", genus)
#     }
#   }
#   taxa <- c(taxa, taxon)
# }
# taxa <- data.table(taxa)
# taxa <- taxa[, Index := 1:.N , taxa ]
# taxa$taxa <- paste(taxa$taxa, taxa$Index, sep="_")
# taxa$taxa <- gsub("_1", "", taxa$taxa)
# #rownames(x) <- taxa$taxa
# correlation.table$X1 <- taxa$taxa

cortable2 <- correlation.table[correlation.table$p.adj < 0.05,]

name <- "Taxa_Trans_Heatmap.pdf"
fp <- paste(out_dir, name, sep="/")
pdf(fp, height=22,width=11)
heat(cortable2, "X1", "X2", fill = "Correlation") 
dev.off()

######################################################################
# For pathways
#Heat map of transcripts/pathways and OTUs

#Make the names of the samples the same for transcripts and otu table
samples2 <- rownames(mapping[mapping$SampleID %in% samples & mapping$transcriptomic %in% colnames(ileum_transcript),])
samples3 <- mapping[samples2,"transcriptomic"]

x <- as.matrix(t(CLR_otutable))
rownames(x) <- taxonomy3$V2
Cor_out2 <- cor(t(x))
Cor_out3 <- abs(Cor_out2) > 0.95
new_taxatable <- x

for(i in 1:nrow(Cor_out3)){
  for(j in i:nrow(Cor_out3)){
    if(i==j){
      next
    }
    if(is.na(Cor_out3[j,i])){
      next
    }
    if(Cor_out3[j,i]){
      new_taxatable[j,] <- new_taxatable[j,] + new_taxatable[i,]
      new_taxatable[i,] <- rep(0, ncol(new_taxatable))
    }
  }
}

#Collapse
taxa_table2 <- new_taxatable[rowSums(new_taxatable) > 0, ]
x <- taxa_table2


x <- t(x[,samples2])
# replace rownames with unique ones
otu_names <- as.character(colnames(x))
taxa <- c()
for (i in 1:length(otu_names)){
  genus <- strsplit(otu_names[i], "; ", fixed=T)[[1]][6]
  family <- strsplit(otu_names[i], "; ", fixed=T)[[1]][5]
  order <- strsplit(otu_names[i], "; ", fixed=T)[[1]][4]
  class <- strsplit(otu_names[i], "; ", fixed=T)[[1]][3]
  phylum <- strsplit(otu_names[i], "; ", fixed=T)[[1]][2]
  if(family== "f__"){
    taxon <- gsub("o__", "", order)
    taxon <- paste("order:", taxon)
  } else {
    if(genus == "g__"){
      taxon <- gsub("f__", "", family)
      taxon <- paste("family:", taxon)
    } else {
      taxon <- gsub("g__", "", genus)
    }
  }
  taxa <- c(taxa, taxon)
}
taxa <- data.table(taxa)
taxa <- taxa[, Index := 1:.N , taxa ]
taxa$taxa <- paste(taxa$taxa, taxa$Index, sep="_")
taxa$taxa <- gsub("_1", "", taxa$taxa)
colnames(x) <- taxa$taxa


y <- reactome_samples[samples3,]


#correlations <- associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
# Or, alternatively, the same output is also available in a handy table format
correlation.table2 <- associate(x, y, method = "pearson", mode = "table", p.adj.threshold = 0.5, n.signif = 1)

# replace rownames with unique ones
#otu_names <- as.character(rownames(x))
# otu_names <- as.character(correlation.table$X1)
# taxa <- c()
# for (i in 1:length(otu_names)){
#   genus <- strsplit(otu_names[i], "; ", fixed=T)[[1]][6]
#   family <- strsplit(otu_names[i], "; ", fixed=T)[[1]][5]
#   order <- strsplit(otu_names[i], "; ", fixed=T)[[1]][4]
#   class <- strsplit(otu_names[i], "; ", fixed=T)[[1]][3]
#   phylum <- strsplit(otu_names[i], "; ", fixed=T)[[1]][2]
#   if(family== "f__"){
#     taxon <- gsub("o__", "", order)
#     taxon <- paste("order:", taxon)
#   } else {
#     if(genus == "g__"){
#       taxon <- gsub("f__", "", family)
#       taxon <- paste("family:", taxon)
#     } else {
#       taxon <- gsub("g__", "", genus)
#     }
#   }
#   taxa <- c(taxa, taxon)
# }
# taxa <- data.table(taxa)
# taxa <- taxa[, Index := 1:.N , taxa ]
# taxa$taxa <- paste(taxa$taxa, taxa$Index, sep="_")
# taxa$taxa <- gsub("_1", "", taxa$taxa)
# #rownames(x) <- taxa$taxa
# correlation.table$X1 <- taxa$taxa

cortable3 <- correlation.table2[correlation.table2$p.adj < 0.25,]

name <- "Taxa_Pathways_Heatmap.pdf"
fp <- paste(out_dir, name, sep="/")
pdf(fp, height=8,width=8)
heat(cortable3, "X1", "X2", fill = "Correlation") 
dev.off()


