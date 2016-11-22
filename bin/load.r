library("vegan")
library("RColorBrewer")
library("vegan")
library("ggplot2")
library("reshape2")
library("plyr")
library('beeswarm')
library('ape')
library("grid")
library("gridExtra")

source('bin/output_dir.r')

#Find the data
data_dir <- "data/"
otu_fp <- paste(data_dir, "Turkey_1000_Taxa.txt", sep='')
map_fp <- paste(data_dir, "Turkey_mapping.txt", sep='')
alpha_fp <- paste(data_dir, "Alpha_Div1000.txt", sep='')

####Load OTU table and metadata####
#Table is pre-filtered to keep samples with >= 1000 counts
#this is 551 samples out of 603 listed in the mapping file
metadata <- read.table(map_fp, sep = "\t", comment="", header=T, as.is=T, check.names=F)
colnames(metadata)[1] <- "SampleID"
rownames(metadata) <- metadata$SampleID

#otu table is otus as rows, samples as columns, column 553 is the taxonomy (header= taxonomy)
otu_table <- read.table(otu_fp, sep="\t", comment="", header=T, skip=1, as.is=T, check.names=F)
rownames(otu_table) <- otu_table[,1]

#remove and store taxonomy
taxonomy <- as.data.frame(cbind(otu_table[,1], otu_table[,ncol(otu_table)]))
rownames(taxonomy) <- taxonomy$V1
otu_table <- otu_table[,2:(ncol(otu_table)-1)]
colnames(otu_table)[550] <- "TJPBXFMB11Inoc" 

####Filter and normalize####

#keep otus that occur in > 1 sample (9215 OTUs to 5212)
otutable2 <- otu_table[rowSums(otu_table > 0) > 1,]

#square root transform
otutable3 <- sqrt(otutable2)

#convert to relative abundance
otutable4 <- sweep(otutable3,2,colSums(otutable3),`/`)  #OTU table is 142 taxa and 258 samples

####Filter the metadata to keep only samples in the new otutable
ids_keep <- intersect(rownames(metadata), colnames(otutable4))
mapping <- metadata[ids_keep,]
mapping$MapDepth <- as.factor(colSums(otutable2))

####Add taxonomy at species level to OTU table####
t_table <- as.matrix(otutable4)
taxonomy2 <- taxonomy[intersect(rownames(taxonomy), rownames(t_table)),]
rownames(t_table) <- taxonomy2$V2

#Collapse same taxonomies (from 5212 to 767)
taxa_table <- aggregate(t_table[,1:ncol(t_table)], by=list(rownames(t_table)), FUN = sum)
rownames(taxa_table) <- taxa_table[,1]
taxa_table <- taxa_table[,2:ncol(taxa_table)]
mapping$nTaxa <- as.factor(colSums(taxa_table > 0))

#Load alpha_div
alpha_div <- read.table(alpha_fp, sep="\t", comment="", row.names=1, header=TRUE)

####Get sample IDs for testing####
mapping[is.na(mapping)] <- "N_A"
mapping["TJPBXBaseInoc","Tissue"] <- "Base"
mapping["TJPBXFMB11Inoc","Tissue"] <- "FMB11"
mapping["TJPBXTJInoc","Tissue"] <- "TJ"
mapping["TJPBXBaseInoc","Probiotic"] <- "N-A"

Ileum <- rownames(mapping[(mapping$Tissue =="Ileum"),])
Trachea <- rownames(mapping[mapping$Tissue =="Trachea",])
Ceca <- rownames(mapping[mapping$Tissue =="Ceca",])

NoInoc <- row.names(mapping[mapping$Treatment == "1_NoInoc",])
FMB11 <- row.names(mapping[mapping$Treatment == "3_FMB11GroGel",])
BMD <- row.names(mapping[mapping$Treatment == "5_BMD",])
TJPbx <- row.names(mapping[mapping$Treatment == "4_TJPbxGroGel",])
GroGel <- row.names(mapping[mapping$Treatment == "2_GroGelNoPbx",])

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
Days <- Days[2:(length(Days)-1)]

Bodysites <- list(Ileum, Trachea, Ceca)
names(Bodysites) <- c("Ileum", "Trachea", "Ceca")

Pbx <- list(FMB11, TJPbx)
names(Pbx) <- c("FMB11", "TJPbx")

Treatments <- c(Pbx, list(BMD))
names(Treatments)[3] <- "BMD"

Inputs <- list(FMB11_Inoc, TJPbx_Inoc, Base_Inoc)
names(Inputs) <- c("FMB11_Inoc", "TJPbx_Inoc", "Base_Inoc")

samples_no_con <- c(NoInoc,GroGel, BMD, FMB11, TJPbx)

####Add Taxa quantiles to mapping###

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

cols <- brewer.pal(7,'Paired')
cols2 <- colorRampPalette(cols)
