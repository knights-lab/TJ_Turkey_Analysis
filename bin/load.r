library("vegan")
library("RColorBrewer")
library("ggplot2")
library("reshape2")
library("plyr")
library('beeswarm')
library('ape')
library("grid")
library("gridExtra")
library("cowplot")
library("biom")
library("robCompositions")
library("microbiome")
library("dplyr")
library("tidyr")
library("cluster")

source('bin/output_dir.r')
source('bin/diff.test.r')
source("bin/taxa.sum.func.r")

theme_set(theme_bw())

#Find the data
data_dir <- "data/"
otu_fp <- paste(data_dir, "Turkey_1000_Taxa.txt", sep='')
map_fp <- paste(data_dir, "Turkey_mapping.txt", sep='')
alpha_fp <- paste(data_dir, "alpha_diversity10k.txt", sep='')
ileum_transfp <- paste(data_dir, "Ileum_Normalized_Transcriptomics.txt", sep='')
bray_fp <- paste(data_dir, "beta_diversity10k/bray_curtis_Turkey_10k_Taxa.txt", sep='')
wuni_fp <- paste(data_dir, "beta_diversity10k/weighted_unifrac_Turkey_10k_Taxa.txt", sep='')
uni_fp <- paste(data_dir, "beta_diversity10k/unifrac_Turkey_10k_Taxa.txt", sep='')

####Load OTU table and metadata####
#Table is pre-filtered to keep samples with >= 1000 counts
#this is 551 samples out of 603 listed in the mapping file
metadata <- read.table(map_fp, 
                       sep = "\t", 
                       comment="", 
                       header=T, 
                       as.is=T, 
                       check.names=F)
colnames(metadata)[1] <- "SampleID"
rownames(metadata) <- metadata$SampleID

#otu table is otus as rows, samples as columns, column 553 is the taxonomy (header= taxonomy)
otu_table <- read.table(otu_fp, 
                        sep="\t", 
                        comment="", 
                        header=T, 
                        skip=1, 
                        as.is=T, 
                        check.names=F)
rownames(otu_table) <- otu_table[,1]

#remove and store taxonomy
remove <- otu_table[grep("Zea",otu_table$taxonomy), "taxonomy"]
remove1 <- otu_table[grep("Streptophyta",otu_table$taxonomy), "taxonomy"]
remove2 <- otu_table[grep("Chloroplast",otu_table$taxonomy), "taxonomy"]

remove <- c(remove, remove1)
otu_table <- otu_table[!otu_table$taxonomy %in% remove,]
taxonomy <- as.data.frame(cbind(otu_table[,1], otu_table[,ncol(otu_table)]))
rownames(taxonomy) <- taxonomy$V1
otu_table <- otu_table[,2:(ncol(otu_table)-1)]
colnames(otu_table)[550] <- "TJPBXFMB11Inoc" 

####Filter####
#drop samples below 1500 counts
otu_table <- otu_table[, colSums(otu_table) > 1500]
#keep otus that occur in > 1 sample (9215 OTUs to 606)
#Change OTUs less than 1/10 millionth of read depth to 0
otu_table2 <- otu_table
otu_table2[otu_table2 < sum(colSums(otu_table2))/10000000] <- 0

#Change singletons to 0 
otu_table2[otu_table2 < 2] <- 0 # This table has 9215 OTUs

#Filter the OTU table to keep OTUs in at least 5% of samples, and whose count is at least 5
otu_table3 <- otu_table2[rowSums(otu_table2 > 0) > (0.05*ncol(otu_table2)),] #this table has 606 OTUs

otutable4 <- otu_table3

####Transform####
#Center log-ratio tranform the data for diff. taxa testing
#Ref: Palarea-Albaladejo J, et al. 2014. JOURNAL OF CHEMOMETRICS. A bootstrap estimation scheme for chemical compositional data with nondetects. 28;7:585–599.
#Ref: Gloor GB, et al. 2016. ANNALS OF EPIDEMIOLOGY. It's all relative: analyzing microbiome data as compositions. 26;5:322-329.
CLR_otutable <- otu_table3
CLR_otutable[CLR_otutable == 0] <- 0.65 #Convert any 0 to 0.65 to allow for CLR transform

CLR_otutable <- t(CLR_otutable) #convert to samples as rows
CLR_otutable <- cenLR(CLR_otutable)$x.clr  #transform

###Filter the metadata to keep only samples in the new otutable
ids_keep <- intersect(rownames(metadata), colnames(otutable4))
mapping <- metadata[ids_keep,]
mapping$MapDepth <- as.factor(colSums(otu_table[,ids_keep]))

####Add the transcriptomics sample IDs to match up with
mapping$transcriptomic <- NA
mapping$treatment4 <- mapping$Treatment2
mapping$treatment4 <- lapply(mapping$treatment4, sub, pattern="No_Inoc", replacement="NOI",fixed=T)
mapping$treatment4 <- lapply(mapping$treatment4, sub, pattern="TJPbx", replacement="TJP",fixed=T)
mapping$treatment4 <- lapply(mapping$treatment4, sub, pattern="FMB11", replacement="FMB",fixed=T)
mapping$treatment4 <- lapply(mapping$treatment4, sub, pattern="GroGel", replacement="GRG",fixed=T)
mapping$treatment4 <- unlist(mapping$treatment4)

substrRIGHT <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
for(r in 1:nrow(mapping)){
  mapping[r,"transcriptomic"] <- paste(mapping$treatment4[r], mapping$Collection_Day[r], mapping$Cage[r], substrRIGHT(mapping$SampleID[r], 3),sep="")
}

####Add taxonomy at species level to OTU table####
t_table <- as.matrix(otutable4)
taxonomy2 <- taxonomy[intersect(rownames(taxonomy), rownames(t_table)),]
rownames(t_table) <- taxonomy2$V2

#Collapse same taxonomies (from 606 to 99)
taxa_table <- aggregate(t_table[,1:ncol(t_table)], by=list(rownames(t_table)), FUN = sum)
rownames(taxa_table) <- taxa_table[,1]
taxa_table <- taxa_table[,2:ncol(taxa_table)]
mapping$nTaxa <- as.factor(colSums(taxa_table > 0))

#Load alpha_div
alpha_div <- read.table(alpha_fp, sep="\t", comment="", row.names=1, header=TRUE)

#Load Ileum transcriptomics
ileum_transcript <- read.table(ileum_transfp, sep='\t', comment="", row=1, header=T)
#These samples aren't in the map:
#Day1 samples, and two Day 3 samples
print("Transcript samples missing:")
colnames(ileum_transcript)[which(!colnames(ileum_transcript) %in% mapping$transcriptomic)]

####Get sample IDs for testing####
mapping[is.na(mapping)] <- "N_A"
mapping["TJPBXBaseInoc","Tissue"] <- "Base"
mapping["TJPBXFMB11Inoc","Tissue"] <- "FMB11"
mapping["TJPBXTJInoc","Tissue"] <- "TJ"
mapping["TJPBXBaseInoc","Probiotic"] <- "N-A"

Ileum <- rownames(mapping[(mapping$Tissue =="Ileum"),])
Trachea <- rownames(mapping[mapping$Tissue =="Trachea",])
Ceca <- rownames(mapping[mapping$Tissue =="Ceca",])

NoInoc <- row.names(mapping[mapping$Treatment2 == "No_Inoc",])
FMB11 <- row.names(mapping[mapping$Treatment2 == "FMB11",])
BMD <- row.names(mapping[mapping$Treatment2 == "BMD",])
TJPbx <- row.names(mapping[mapping$Treatment2 == "TJPbx",])
GroGel <- row.names(mapping[mapping$Treatment2 == "GroGel",])

Day0 <- row.names(mapping[mapping$Treatment == "D0Pool",])
Ileum <- Ileum[which(! Ileum == Day0)]
Base_Inoc <- row.names(mapping[mapping$Treatment == "Base_Inoc",])
FMB11_Inoc <- row.names(mapping[mapping$Treatment == "FMB11_Inoc",])
TJPbx_Inoc <- row.names(mapping[mapping$Treatment == "TJ_Inoc",])

Antibiotics <- row.names(mapping[mapping$Abx == "Yes",])
Probiotics <- row.names(mapping[mapping$Probiotic == "Yes",])
Controls <- row.names(mapping[mapping$Probiotic == "No",])
Controls <- Controls[which(!Controls %in% Antibiotics)]
Treated <- list(Antibiotics, Probiotics)
names(Treated) <- c("Antibiotics", "Probiotics")

Days_avail <- unique(mapping$Collection_Day)
Days <- list()
for(i in 1:length(Days_avail)){
  working_day <- Days_avail[i]
  Day_samples <- list(rownames(mapping[mapping$Collection_Day == working_day,]))
  Days <- c(Days, Day_samples)
  names(Days)[i] <- working_day
}
keep_days <- c("D03", "D06", "D13")
Days <- Days[names(Days) %in% keep_days]

Bodysites <- list(Trachea, Ileum, Ceca)
names(Bodysites) <- c("Trachea", "Ileum", "Cecum")

Pbx <- list(FMB11, TJPbx)
names(Pbx) <- c("FMB11", "TJPbx")

Treatments <- c(Pbx, list(BMD))
names(Treatments)[3] <- "BMD"

Inputs <- list(FMB11_Inoc, TJPbx_Inoc, Base_Inoc)
names(Inputs) <- c("FMB11_Inoc", "TJPbx_Inoc", "Base_Inoc")

samples_no_con <- c(NoInoc,GroGel, BMD, FMB11, TJPbx)

treats_nocon <- list(BMD, TJPbx, FMB11, GroGel) 
names(treats_nocon) <- c("BMD", "TJPbx", "FMB11", "GroGel")
####Add Taxa quantiles to mapping###
# 
# ranges <- c(0, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# mapping <- cbind(mapping, t(taxa_table))
# taxa <- colnames(t(taxa_table))
# for(i in 1:nrow(mapping)){
#   sample_id <- rownames(mapping)[i]
#   for(k in 1:length(taxa)){
#     taxon <- taxa[k]
#     for(m in 1:length(ranges)){
#       if((mapping[sample_id, taxon] >= ranges[m]) && (mapping[sample_id, taxon] < ranges[(m+1)])){
#         mapping[sample_id, taxon] <- ranges[m]
#         next
#       } 
#     }
#   }
# }
# colnames(mapping) <- gsub(" ", "_", colnames(mapping))

####Set Colors####
cols <- brewer.pal(8,'Dark2')
cols2 <- colorRampPalette(cols)
