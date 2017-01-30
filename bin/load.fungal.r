
#Paths to funal tables - alpha and beta generated with RA tables in QIIME
#OTU table picked with ninja
fotu_fp <- paste(data_dir, "Fungal_otutable.biom", sep='')
fmap_fp <- paste(data_dir, "Turkey_mapping.txt", sep='')
falpha_fp <- paste(data_dir, "Alpha_fungal_RA.txt", sep='')
fbray_fp <- paste(data_dir, "bray_curtis_Turkey_fungal_RA.txt", sep='')

alpha_div2 <- read.table(falpha_fp, sep="\t", comment="", row.names=1, header=TRUE)
fbray <- read.table(fbray_fp, sep="\t", header=T, row=1)

####Load OTU table and metadata####

#otu table is OTU ID (rows) x sample ID (columns)
#dimensions = 381 otus and  244 samples
fungal_otu <- as.matrix(biom_data(read_biom(fotu_fp)))

ftaxonomy <- observation_metadata(read_biom(fotu_fp))
ftaxonomy$taxonomy7 <- sapply(strsplit(ftaxonomy$taxonomy7, split='__', fixed=TRUE), function(x) (x[2]))

fungal_map <- read.table(fmap_fp, sep = "\t", comment="", header=T, as.is=T, check.names=F)
colnames(fungal_map)[1] <- "SampleID"
rownames(fungal_map) <- fungal_map$SampleID

#Remove positive and negative controls
controls <- c("Negative", "Positive", "Positive.KAPA.PROG")
fungal_otu1 <- fungal_otu[, ! colnames(otutable) %in% controls]
#Store the positive and negative values
control_table <- cbind(fungal_otu[,"Negative"], fungal_otu[,"Positive"], fungal_otu[,"Positive.KAPA.PROG"])
colnames(control_table) <- controls

#keep singletons that occur only in more than one sample (381 to 205)
fungal_otu2 <- fungal_otu1[rowSums(fungal_otu1 > 0) > 1,]
#keep samples with at least 50 sequence counts (241 to 222), and
fungal_otu3 <- fungal_otu2[,colSums(fungal_otu2) > 50]

#convert to relative abundance
fungal_RA <- sweep(fungal_otu3,2,colSums(fungal_otu3),`/`)

####Filter the metadata to keep only samples in the new otutable
ids_keep <- intersect(rownames(fungal_map), colnames(fungal_RA))
f_map <- fungal_map[ids_keep,]
#There is a sample in the otu table thats not in the map: "TJPBX.D0I"
fungal_RA <- fungal_RA[,ids_keep]
fungal_otu3 <- fungal_otu3[, ids_keep]
f_map$MapDepth <- as.factor(colSums(fungal_otu3))

#Print this table as the normalized OTU table, multiplied by 10000 for beta div problems
fungal_RA <- round(fungal_RA * 100000)
sink("data/Turkey_fungal_RA.txt")
cat("#OTUID")
write.table(fungal_RA, 
            sep="\t", #tell R to make is tab-delimited
            quote=F, #tell R not to put quotes
            col.names=NA) #formats the column headers properly
sink()

#There are no controls left, after filtering and all samples are ileum

#Store treatments, etc.
NoInoc_f <- row.names(f_map[f_map$Treatment == "1_NoInoc",])
FMB11_f <- row.names(f_map[f_map$Treatment == "3_FMB11GroGel",])
BMD_f <- row.names(f_map[f_map$Treatment == "5_BMD",])
TJPbx_f <- row.names(f_map[f_map$Treatment == "4_TJPbxGroGel",])
GroGel_f <- row.names(f_map[f_map$Treatment == "2_GroGelNoPbx",])

Antibiotics_f <- row.names(f_map[f_map$Abx == "Yes",])
Probiotics_f <- row.names(f_map[f_map$Probiotic == "Yes",])
Controls_f <- row.names(f_map[f_map$Probiotic == "No",])
Controls_f <- Controls_f[which(!Controls_f %in% Antibiotics_f)]
Treated_f <- list(Antibiotics_f, Probiotics_f)
names(Treated_f) <- c("Antibiotics", "Probiotics")

Days_avail_f <- unique(f_map$Collection_Day)
Days_f <- list()
for(i in 1:length(Days_avail_f)){
  working_day <- Days_avail_f[i]
  Day_samples <- list(rownames(f_map[f_map$Collection_Day == working_day,]))
  Days_f <- c(Days_f, Day_samples)
  names(Days_f)[i] <- working_day
}
keep_days <- c("D03", "D06", "D13")
Days_f <- Days_f[names(Days_f) %in% keep_days]

Pbx_f <- list(FMB11_f, TJPbx_f)
names(Pbx_f) <- c("FMB11", "TJPbx")

Treatments_f <- c(Pbx_f, list(BMD_f))
names(Treatments_f)[3] <- "BMD"


