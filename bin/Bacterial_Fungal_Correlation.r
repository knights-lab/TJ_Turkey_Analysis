#Fungal-Bacterial Correlation

######################################################################
#procrustes of the fungal and bacterial components
fotu_table <- fungal_otu3
#keep otus that occur in > 1 sample (9215 OTUs to 606)
#Change OTUs less than 1/10 millionth of read depth to 0
fotu_table2 <- fotu_table
fotu_table2[fotu_table2 < sum(colSums(fotu_table2))/10000000] <- 0

#Change singletons to 0 
fotu_table2[fotu_table2 < 2] <- 0 # This table has 9215 OTUs

#Filter the OTU table to keep OTUs in at least 5% of samples, and whose count is at least 5
fotu_table3 <- fotu_table2[rowSums(fotu_table2 > 0) > (0.05*ncol(fotu_table2)),] #this table has 606 OTUs
fotutable4 <- fotu_table3

####Transform####
#Center log-ratio tranform the data for diff. taxa testing
#Ref: Palarea-Albaladejo J, et al. 2014. JOURNAL OF CHEMOMETRICS. A bootstrap estimation scheme for chemical compositional data with nondetects. 28;7:585–599.
#Ref: Gloor GB, et al. 2016. ANNALS OF EPIDEMIOLOGY. It's all relative: analyzing microbiome data as compositions. 26;5:322-329.
fCLR_otutable <- fotu_table3
fCLR_otutable[fCLR_otutable == 0] <- 0.65 #Convert any 0 to 0.65 to allow for CLR transform

fCLR_otutable <- t(fCLR_otutable) #convert to samples as rows
fCLR_otutable <- cenLR(fCLR_otutable)$x.clr  #transform

day6_f <- fCLR_otutable[Days_f[[2]], ]
keep_f <- fungal_otu3[rowSums(fungal_otu3 > 0) > 5, rownames(day6_f)]
keep_f <- intersect(rownames(keep_f), colnames(day6_f))
day6_f <- day6_f[rowSums(day6_f > 0) > 1, keep_f]
test <- t(day6_f)

#Fungal distances
beta_table <- as.matrix(dist(t(test)))


#Bacteria Ileum day 6
otu_test <- read.table("data/Turkey_10k_Taxa.txt", 
                       sep="\t", 
                       comment="", 
                       header=T, 
                       skip=1, 
                       as.is=T, 
                       check.names=F, row=1)
#remove taxa that are plants
remove <- as.character(taxonomy[grep("zea",taxonomy$V2), "V1"])
remove1 <- as.character(taxonomy[grep("Streptophyta",taxonomy$V2), "V1"])
remove2 <- as.character(taxonomy[grep("Chloroplast",taxonomy$V2), "V1"])
remove_these <- c(remove, remove1, remove2)
otu_test <- otu_test[!rownames(otu_test) %in% remove_these,]

tree = read_tree_greengenes('data/97_otus.tree')
ileum_6_bacteria <- intersect(rownames(mapping[Days[[2]],]), Ileum)

beta_keeps <- intersect(ileum_6_bacteria, colnames(otu_test))
otu_test1 <- otu_test[,beta_keeps]
otu_test1 <- otu_test1[rowSums(otu_test1 > 0) > 1, ]
m <- mapping[beta_keeps,]
#generate phyloseq objects
OTU <- otu_table(otu_test1, taxa_are_rows=T)
sampledata = sample_data(m)
#manditory input for UniFrac
physeq1 = phyloseq(OTU, sampledata, tree)
#distances
uni1 <- as.matrix(UniFrac(physeq1, weighted=FALSE, normalized=TRUE))

#find samples that overlap:
overlap_samples <- intersect(rownames(uni1), rownames(beta_table))

#Generate PCoAs
PCOA_b <- pcoa(uni1[overlap_samples, overlap_samples])$vectors
PCOA_f <- pcoa(beta_table[overlap_samples,overlap_samples])$vectors
var_exp <- pcoa(beta_table[overlap_samples,overlap_samples])$values

pro <- procrustes(PCOA_b, PCOA_f)
pro <- procrustes(PCOA_b[,1:3], PCOA_f[,1:3])
summary(pro)
pro2 <- protest(PCOA_b[,1:3], PCOA_f[,1:3], permutations = 999)
sink(paste(main_fp, "/Figure6_mantel.txt", sep=""))
print(mantel(uni1[overlap_samples, overlap_samples], beta_table[overlap_samples,overlap_samples], method="pearson", permutations=999))
sink()

bact_pro <- data.frame(pro$X)
fungi_pro <- data.frame(pro$Yrot)

rownames(fungi_pro) <- overlap_samples
colnames(fungi_pro) <- colnames(bact_pro)

fungi_pro$treatment <- mapping[overlap_samples, "treatment4"]
fungi_pro$type <- "fungi"
bact_pro$treatment <- mapping[overlap_samples, "treatment4"]
bact_pro$type <- "bacteria"

bact_pro$sample <- c(1:35)
fungi_pro$sample <- c(1:35)

print(paste( "P=", pro2$signif))
print(paste( "M=", pro2$ss))
p_val <- pro2$signif
m_val <- pro2$ss

PCOA <- rbind(bact_pro, fungi_pro)
PCOA$SampleID <- rownames(PCOA)
colnames(PCOA) <- gsub("Axis.", "PC", colnames(PCOA))

cols_pro <- c("#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#E6AB02")
PCOA$treatment <- factor(PCOA$treatment, levels=c("BMD", "NOI", "TJP", "FMB", "GRG"))
pro_b_f <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "treatment", shape="type")) +
  geom_line(aes(x= PC1, y=PC2, group=sample, color=treatment)) +
  scale_shape_manual(values=c(16, 15)) +
  theme(plot.title = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "in"),
        legend.text = element_text(size=5),
        legend.position = 'bottom',
        axis.text = element_text(size=5),
        axis.title = element_text(size=8)) +
  scale_color_manual(values=cols_pro) +
  annotate("text", x=0.2, y=0.25, label= paste("P=", p_val), size=2) +
  annotate("text", x=0.2, y=0.20, label= paste("M2=", round(m_val, digits=3)), size=2) +
  #coord_fixed(xlim = c(-0.32, 0.41), ylim = c(-0.32, 0.41)) +
  labs(x=paste("PC1 (", round(var_exp$Relative_eig[1], digits=3)*100, "%)", sep=""), y= paste("PC2 (", round(var_exp$Relative_eig[2], digits=3)*100, "%)", sep=""))


######################################################################
#Are fungi correlated with bacteria?
f_otu <- data.frame(t(test[, overlap_samples]))
b_otu <- data.frame(CLR_otutable[overlap_samples,!colnames(CLR_otutable) %in% remove_these])
colnames(b_otu) <- gsub("X", "", colnames(b_otu))

#Within BMD test for correlation
f_otu_ABX <- f_otu[rownames(f_otu)%in% BMD,]
b_otu_ABX <- b_otu[rownames(b_otu)%in% BMD,]

keep_b <- otu_table3[,rownames(b_otu_ABX)]
keep_b <- keep_b[rowSums(keep_b)<200 & rowSums(keep_b) > 5,] #334 OTUs

keep_f <- data.frame(fungal_otu3[,rownames(f_otu_ABX)])
keep_f <- keep_f[rowSums(keep_f > 0) > 2,] #334 OTUs

f_otu_ABX <- f_otu_ABX[,colnames(f_otu) %in% rownames(keep_f)]
b_otu_ABX <- b_otu_ABX[,colnames(b_otu) %in% rownames(keep_b)]

colnames(b_otu_ABX) <- taxonomy3[colnames(b_otu_ABX), "V2"]
otu_names <- as.character(colnames(b_otu_ABX))
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
colnames(b_otu_ABX) <- taxa$taxa
colnames(f_otu_ABX) <- ftaxonomy3[colnames(f_otu_ABX), "taxonomy7"]

correlation.tableAB <- associate(f_otu_ABX, b_otu_ABX, method = "pearson", mode = "table", p.adj.threshold = 0.5, n.signif = 1)
within_abx <- heat(correlation.tableAB, "X1", "X2", fill = "Correlation", p.adj.threshold = 0.25, star='p.adj', star.size=3,
                   colours=c("purple4", "#5E3C99", "white", "#f5bf99", "#e66101")) 

#Within TJPbx test for correlation
f_otu_ABX <- f_otu[rownames(f_otu)%in% TJPbx,]
b_otu_ABX <- b_otu[rownames(b_otu)%in% TJPbx,]

keep_b <- otu_table3[,rownames(b_otu_ABX)]
keep_b <- keep_b[rowSums(keep_b)<200 & rowSums(keep_b) > 5,] #334 OTUs

keep_f <- data.frame(fungal_otu3[,rownames(f_otu_ABX)])
keep_f <- keep_f[rowSums(keep_f > 0) > 2,] #334 OTUs

f_otu_ABX <- f_otu_ABX[,colnames(f_otu) %in% rownames(keep_f)]
b_otu_ABX <- b_otu_ABX[,colnames(b_otu) %in% rownames(keep_b)]

# replace rownames with unique ones
colnames(b_otu_ABX) <- taxonomy3[colnames(b_otu_ABX), "V2"]
otu_names <- as.character(colnames(b_otu_ABX))
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
colnames(b_otu_ABX) <- taxa$taxa

colnames(f_otu_ABX) <- ftaxonomy3[colnames(f_otu_ABX), "taxonomy7"]

correlation.tableTP <- associate(f_otu_ABX, b_otu_ABX, method = "pearson", mode = "table", p.adj.threshold = 0.5, n.signif = 1)
within_pbx <- heat(correlation.tableTP, "X1", "X2", fill = "Correlation", p.adj.threshold = 0.25, star='p.adj', star.size=3, 
                   colours=c("purple4", "#5E3C99", "white", "#f5bf99", "#e66101")) 


#Within Control1 test for correlation
f_otu_ABX <- f_otu[rownames(f_otu)%in% NoInoc,]
b_otu_ABX <- b_otu[rownames(b_otu)%in% NoInoc,]

keep_b <- otu_table3[,rownames(b_otu_ABX)]
keep_b <- keep_b[rowSums(keep_b)<200 & rowSums(keep_b) > 5,] #334 OTUs

keep_f <- data.frame(fungal_otu3[,rownames(f_otu_ABX)])
keep_f <- keep_f[rowSums(keep_f > 0) > 2,] #334 OTUs

f_otu_ABX <- f_otu_ABX[,colnames(f_otu) %in% rownames(keep_f)]
b_otu_ABX <- b_otu_ABX[,colnames(b_otu) %in% rownames(keep_b)]

# replace rownames with unique ones
colnames(b_otu_ABX) <- taxonomy3[colnames(b_otu_ABX), "V2"]
otu_names <- as.character(colnames(b_otu_ABX))
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
colnames(b_otu_ABX) <- taxa$taxa

colnames(f_otu_ABX) <- ftaxonomy3[colnames(f_otu_ABX), "taxonomy7"]

#Error - bc none are significant
#correlation.table <- associate(f_otu_ABX, b_otu_ABX, method = "pearson", mode = "table", p.adj.threshold = 0.1, n.signif = 1)
#within_con <- heat(correlation.table, "X1", "X2", fill = "Correlation", p.adj.threshold = 0.1, star='p.adj', star.size=1) 
#Double Check by hand:
x <- f_otu_ABX
y <- b_otu_ABX

dat <- NULL

for(i in colnames(x)){ # OTUs
  for(j in colnames(y)){ # Pathways
    a<-as.numeric(x[,i])
    b<-as.numeric(y[,j])
    tmp <- c(i,j,cor(a,b, method = "pearson", use = "everything"), 
             cor.test(a,b,method = "pearson")$p.value)
    if(is.null(dat)){
      dat <- tmp
    }
    else{
      dat<-rbind(dat,tmp)
    }
  }
}

dat<- data.frame(row.names = NULL, dat, stringsAsFactors = FALSE)

colnames(dat)<-c("X1","X2","Correlation","p.val")
dat$Correlation <- as.numeric(dat$Correlation)
dat$p.adj <- p.adjust(dat$p.val, 'fdr')
dat <- data.frame(dat)
#** None are significant 

#Within Control2 test for correlation - can't do bc there are only two GroGel in B
f_otu_ABX <- f_otu[rownames(f_otu)%in% GroGel,]
b_otu_ABX <- b_otu[rownames(b_otu)%in% GroGel,]

##############
#Compile figure
fungal_bacteria <- plot_grid(pro_b_f, within_abx, within_pbx, ncol=3, rel_widths = c(1,1.8,1.8))
pdf(paste(main_fp, "/Supplemental_figure6.pdf", sep=""), width=20, height=5, useDingbats = F)
fungal_bacteria
dev.off()
