#Transcriptome changes by Treatment
ileum_transcript <- read.table(ileum_transfp, sep='\t', comment="", row=1, header=T)

#load the mapping file
i_map <- read.table("data/transcript_map.txt", sep='\t', header=T, comment="", row=1)
i_map$SampleID <- rownames(i_map)

#samples <- intersect(Ileum,Days[[2]])
#samples2 <- rownames(mapping[mapping$SampleID %in% samples & mapping$transcriptomic %in% colnames(ileum_transcript),])
#samples3 <- mapping[samples2,"transcriptomic"]
keeps <- rownames(i_map[i_map$Collection_Day == "D06",])
trans1 <- ileum_transcript[,keeps]

#colnames(trans1) <- mapping[samples2, "SampleID"]
trans <- vegdist(t(trans1))

PCOA2 <- pcoa(trans)$vectors
var_exp <- pcoa(trans)$values

#Run stats for diff. centroids
beta_dist = as.dist(trans)
ad = adonis(trans ~ i_map[keeps,"Treatment"], data=i_map, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
#Run Stats for diff. dispersion
beta_out <- betadisper(trans, i_map[keeps,"Treatment"])
p_val_disp <- permutest(beta_out)$tab[1, 6]

for(c in 1:ncol(PCOA2)){
  colnames(PCOA2)[c] <- paste("PC",c, sep="")
}

PCOA2 <- cbind(PCOA2, rownames(PCOA2))
colnames(PCOA2)[ncol(PCOA2)] <- "SampleID"
PCOA2 <- as.data.frame(PCOA2)    

PCOA2$PC1 <- as.numeric(levels(PCOA2$PC1))[PCOA2$PC1]
PCOA2$PC2 <- as.numeric(levels(PCOA2$PC2))[PCOA2$PC2]

PCOA2 <- cbind(PCOA2, i_map[keeps,])

PCOA2$Treatment <- factor(PCOA2$Treatment, levels = c("BMD", "TJP", "FMB", "NOI", "GRG"))
total_col <- c("#D95F02", "#7570B3", "#E7298A", "#1B9E77", "#E6AB02")
trans_pcoa <- ggplot(PCOA2) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Treatment"), alpha =0.8) + 
  scale_color_manual(values= total_col) +
  scale_fill_manual(values= total_col) +
  stat_ellipse(aes_string(x = "PC1", y = "PC2", color = "Treatment"), linetype = 2) +
  annotate("text", x=0, y=0.2, label= paste("P=", p_val), size=2) +
  annotate("text", x=0, y=0.21, label= paste("R2=", round(r_sq, digits=3)), size=2) +
  labs(x=paste("PC1 (", round(var_exp$Relative_eig[1], digits=3)*100, "%)", sep=""), y= paste("PC2 (", round(var_exp$Relative_eig[2], digits=3)*100, "%)", sep=""))


######################################################################
#Pathway PCoA
#Load reactome data
#Subset to find genes in turkey RNAseq that overlap
reactome <- read.delim("data/ReactomePathways.gmt", sep="\t", header=F)
reactome <- reactome[reactome$V3 == "Reactome Pathway",] #subset just to pathways
reactome <- as.data.frame(apply(reactome, 2, function(x) as.character(x)), stringsAsFactors = F)
x <- unique(rapply(reactome[,4:44],function(x) unique(x))) #store all possible transcript ids
keeps <- intersect(x,rownames(ileum_transcript))

transcript_sub <- ileum_transcript[keeps,] #subset the transcript data

#for each pathway, and each sample, sum the normalized count for the transcripts in that pathway
reactome_samples <- data.frame(matrix(nrow=nrow(reactome), ncol=ncol(ileum_transcript)))
rownames(reactome_samples) <- reactome[,"V1"]
colnames(reactome_samples) <- colnames(ileum_transcript)
reactome_samples <- t(reactome_samples)

for(i in 1:nrow(reactome)){
  IDs <- rapply(reactome[i,4:44], function(xx) as.character(xx))
  reactome_samples[,i] <- colSums(transcript_sub[IDs,], na.rm = T)
}
keeps <- rownames(i_map[i_map$Collection_Day == "D06",])
path1 <- t(reactome_samples[keeps,])

path2 <- vegdist(t(path1))

PCOA2 <- pcoa(path2)$vectors
var_exp <- pcoa(path2)$values

#Run stats for diff. centroids
beta_dist = as.dist(path2)
ad = adonis(path2 ~ i_map[keeps,"Treatment"], data=i_map, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
#Run Stats for diff. dispersion
beta_out <- betadisper(path2, i_map[keeps,"Treatment"])
p_val_disp <- permutest(beta_out)$tab[1, 6]

for(c in 1:ncol(PCOA2)){
  colnames(PCOA2)[c] <- paste("PC",c, sep="")
}

PCOA2 <- cbind(PCOA2, rownames(PCOA2))
colnames(PCOA2)[ncol(PCOA2)] <- "SampleID"
PCOA2 <- as.data.frame(PCOA2)    

PCOA2$PC1 <- as.numeric(levels(PCOA2$PC1))[PCOA2$PC1]
PCOA2$PC2 <- as.numeric(levels(PCOA2$PC2))[PCOA2$PC2]

PCOA2 <- cbind(PCOA2, i_map[keeps,])

centroids <- aggregate(cbind(PCOA2$PC1,PCOA2$PC2) ~ PCOA2$Treatment,PCOA2,mean)
colnames(centroids) <- c("Treatment", "PC1", "PC2")

#pairwise adonis for differences:
pair_ad_list <- list()
for(t in 1:4){
  for(u in (t+1):5){
    t1 <- as.character(unique(i_map$Treatment)[[t]])
    t2 <- as.character(unique(i_map$Treatment)[[u]])
    sams <- rownames(i_map[i_map$Treatment == t1 | i_map$Treatment == t2,])
    sams <- intersect(rownames(as.matrix(path2)), sams)
    beta_dist = as.dist(as.matrix(path2)[sams,sams])
    ad = adonis(beta_dist ~ i_map[sams,"Treatment"], data=i_map, permutations=999)
    p_val <- ad$aov.tab[1,6]
    r_sq <- ad$aov.tab[1,5]
    pair_ad_list[[paste(t1, "_", t2, sep="")]] <- c("pval", p_val, "r_sq", r_sq)
  }
}

PCOA2$Treatment <- factor(PCOA2$Treatment, levels = c("BMD", "TJP", "FMB", "NOI", "GRG"))
total_col <- c("#D95F02", "#7570B3", "#E7298A", "#1B9E77", "#E6AB02")
pathway_pcoa <- ggplot(PCOA2) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Treatment"), alpha =0.8) + 
  scale_color_manual(values= total_col) +
  scale_fill_manual(values= total_col) +
  stat_ellipse(aes_string(x = "PC1", y = "PC2", color = "Treatment"), linetype = 2) +
  coord_fixed(xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2)) +
  annotate("text", x=0, y=0.2, label= paste("P=", p_val), size=2) +
  annotate("text", x=0, y=0.21, label= paste("R2=", round(r_sq, digits=3)), size=2) +
  labs(x=paste("PC1 (", round(var_exp$Relative_eig[1], digits=3)*100, "%)", sep=""), y= paste("PC2 (", round(var_exp$Relative_eig[2], digits=3)*100, "%)", sep=""))+
  geom_point(data=centroids, aes(x=PC1, y=PC2, fill=Treatment), shape = 23, size=3, alpha=0.8, show.legend = F)

sink(paste(main_fp, "treatment_pairwise_centroids_pathways.txt", sep="/"), append =T)
print(pair_ad_list)
sink()

######################################################################
#Procrustes of transcripts and beta div uni
#Bacteria Ileum day 6
otu_test <- read.table("data/Turkey_10k_Taxa.txt", 
                       sep="\t", 
                       comment="", 
                       header=T, 
                       skip=1, 
                       as.is=T, 
                       check.names=F, row=1)
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

#Take beta table of transcripts and find samples that overlap:
samples2 <- rownames(mapping[rownames(mapping) %in% beta_keeps & mapping$transcriptomic %in% rownames(reactome_samples),])
samples3 <- mapping[samples2,"transcriptomic"]
trans1 <- ileum_transcript[,samples3]
colnames(trans1) <- samples2
trans <- as.matrix(vegdist(t(trans1)))

#Generate PCoAs
PCOA_b <- pcoa(uni1[samples2, samples2])$vectors
PCOA_t <- pcoa(trans[samples2,samples2])$vectors

pro <- procrustes(PCOA_b, PCOA_t)
pro2 <- protest(PCOA_b, PCOA_t, permutations = how(nperm = 999))

bact_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)

rownames(trans_pro) <- samples2
colnames(trans_pro) <- colnames(bact_pro)

trans_pro$treatment <- mapping[samples2, "treatment4"]
trans_pro$type <- "transcript"
bact_pro$treatment <- mapping[samples2, "treatment4"]
bact_pro$type <- "bacteria"

bact_pro$sample <- c(1:13)
trans_pro$sample <- c(1:13)

print(paste( "P=", pro2$signif))

PCOA <- rbind(bact_pro, trans_pro)
PCOA$SampleID <- rownames(PCOA)
colnames(PCOA) <- gsub("Axis.", "PC", colnames(PCOA))

pro_b_t <- ggplot(PCOA) +
  geom_point(size = 4, aes_string(x = "PC1", y = "PC2", color = "treatment", shape="type")) +
  geom_line(aes(x= PC1, y=PC2, group=sample, color=treatment)) +
  theme(plot.title = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "in"),
        legend.text = element_text(size=5),
        legend.position = 'bottom',
        axis.text = element_text(size=5),
        axis.title = element_text(size=8)) +
  scale_color_manual(values=cols2(5)) +
  guides(color=guide_legend(nrow=3))

######################################################################
#Procrustes of pathways and beta div uni
#Bacteria Ileum day 6
otu_test <- read.table("data/Turkey_10k_Taxa.txt", 
                       sep="\t", 
                       comment="", 
                       header=T, 
                       skip=1, 
                       as.is=T, 
                       check.names=F, row=1)
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

#Take beta table of transcripts and find samples that overlap:
samples2 <- rownames(mapping[rownames(mapping) %in% beta_keeps & mapping$transcriptomic %in% rownames(reactome_samples),])
samples3 <- mapping[samples2,"transcriptomic"]
reac1 <- t(reactome_samples[samples3,])
colnames(reac1) <- samples2
reac <- as.matrix(vegdist(t(reac1)))

#Generate PCoAs
PCOA_b <- pcoa(uni1[samples2, samples2])$vectors
PCOA_p <- pcoa(reac[samples2,samples2])$vectors
var_exp <- pcoa(reac[samples2,samples2])$values

pro <- procrustes(PCOA_b, PCOA_p)
pro <- procrustes(PCOA_b[,1:3], PCOA_p[,1:3])
summary(pro)
pro2 <- protest(PCOA_b, PCOA_p, permutations = 999)
sink(paste(main_fp, "/Figure5_mantel.txt", sep=""))
print(mantel(uni1[samples2, samples2], reac[samples2,samples2], method="pearson", permutations=999))
sink()

bact_pro <- data.frame(pro$X)
path_pro <- data.frame(pro$Yrot)

rownames(path_pro) <- samples2
colnames(path_pro) <- colnames(bact_pro)

path_pro$treatment <- mapping[samples2, "treatment4"]
path_pro$type <- "pathway"
bact_pro$treatment <- mapping[samples2, "treatment4"]
bact_pro$type <- "bacteria"

bact_pro$sample <- c(1:13)
path_pro$sample <- c(1:13)

print(paste( "P=", pro2$signif))
print(paste( "M=", pro2$ss))
p_val <- pro2$signif
m_val <- pro2$ss

PCOA <- rbind(bact_pro, path_pro)
PCOA$SampleID <- rownames(PCOA)
colnames(PCOA) <- gsub("Axis.", "PC", colnames(PCOA))

PCOA$treatment <- factor(PCOA$treatment, levels=c("BMD", "NOI", "TJP", "FMB", "GRG"))

cols_pro <- c("#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#E6AB02")
pro_b_p <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "treatment", shape="type")) +
  geom_line(aes(x= PC1, y=PC2, group=sample, color=treatment)) +
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
  coord_fixed(xlim = c(-0.32, 0.37), ylim = c(-0.32, 0.37)) +
  labs(x=paste("PC1 (", round(var_exp$Relative_eig[1], digits=3)*100, "%)", sep=""), y= paste("PC2 (", round(var_exp$Relative_eig[2], digits=3)*100, "%)", sep=""))


######################################################################
#Find any differentiated pathways (Abx vs Control, Pbx vs Control)

#Set up tests to run
source("bin/diff.test.r")

keeps <- rownames(i_map[i_map$Collection_Day == "D06",])
trans1 <- reactome_samples[keeps,]

NoInoc_p <- intersect(row.names(i_map[i_map$Treatment== "NOI",]), rownames(trans1))
FMB11_p <- intersect(row.names(i_map[i_map$Treatment == "FMB",]), rownames(trans1))
BMD_p <- intersect(row.names(i_map[i_map$Treatment == "BMD",]), rownames(trans1))
TJPbx_p <- intersect(row.names(i_map[i_map$Treatment == "TJP",]), rownames(trans1))
GroGel_p <- intersect(row.names(i_map[i_map$Treatment == "GRG",]), rownames(trans1))

test.ixs <- list(NoInoc_p, FMB11_p, BMD_p, TJPbx_p, GroGel_p)
names(test.ixs) <- c("NoInoc", "FMB11", "BMD", "TJPbX", "GroGel")

pvals<- c()
ALPHA <- 0.05
diff_list <- list()

for(n in 1:(length(test.ixs)-1)){
  for(m in (n+1):length(test.ixs)){
    test.x <- test.ixs[n]
    test.y <- test.ixs[m]
    working_table <- t(trans1)
    set1 <- intersect(test.x[[1]], colnames(working_table))
    set2 <- intersect(test.y[[1]], colnames(working_table))
    if(length(set1) > 2 && length(set2) > 2){
      full_set <- c(set1, set2)
      test_table <- working_table[, full_set]
      
      test_table <- t(test_table)
      
      #Keep taxa that have a mean rel abundance of at least 0.01
      #test_table <- test_table[,colMeans(test_table)>0.009, drop=F]
      map_test <- i_map[rownames(test_table),]
      difftest <- test.otu.features2(test_table, response=map_test$Treatment, sig.level = 0.10)
      #difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
      difftest$pvals[is.na(difftest$pvals)] <- 1
      cat(paste(names(test.ixs)[n], " vs ", names(test.ixs)[m], "\n"))
      if(any(difftest$pvals <= ALPHA)){
        signif.ix <- which(difftest$pvals <= ALPHA)
        signif.ix <- signif.ix[order(difftest$pvals[signif.ix])]
        #sink("list of diff pathways.txt", append=T)
        cat(paste(names(test.ixs)[n], " vs ", names(test.ixs)[m], "\n"))
        print(table(names(signif.ix)))
        #cat("\n")
        #sink()
        diff_list[[paste(names(test.ixs)[n], names(test.ixs)[m])]] <- names(signif.ix)
      } else {
        cat("Not Signif.")
      }
    }
  }
}

#Only signif = Abx vs No_Inoc
keep_pathways <- unique(diff_list[["NoInoc BMD"]])
x <- data.frame(t(trans1[,keep_pathways]))
rownames(x) <- keep_pathways
x$abx <- 1
x[rownames(x) %in% diff_list[["NoInoc BMD"]], "abx"] <- 0.04

i_table <- data.frame(t(x))
ps <- i_table["abx",]
i_table <- i_table[!rownames(i_table) %in% c("abx"),]
i_map <- mapping[mapping$transcriptomic %in% rownames(i_table),]
rownames(i_map) <- i_map$transcriptomic
keepme <- intersect(rownames(i_map), rownames(i_table))
i_map <- i_map[keepme,]
i_table <- i_table[keepme,]
i_table$Treatment <- i_map$Treatment2
i_table <- melt(i_table)
i_t2 <-aggregate(i_table, by=list(i_table$variable, i_table$Treatment), 
                 FUN=mean, na.rm=TRUE)
ps <- data.frame(t(ps), check.names = F)
ps$Group.1 <- rownames(ps)
i_t2 <- merge(i_t2, ps, by="Group.1")

taxa_order <- i_t2[i_t2$Group.2 == "BMD",]
taxa_order <- taxa_order[order(-taxa_order$value),]
i_t2$Group.1 <- factor(i_t2$Group.1, levels=taxa_order$Group.1)
i_t2$abx <- as.numeric(as.character(i_t2$abx))
i_t2$Group.2 <- factor(i_t2$Group.2, levels = c("BMD", "No_Inoc", "TJPbx", "FMB11", "GroGel"))
pathways_diff <- ggplot(i_t2) + 
  geom_tile(aes(x = Group.2, y = Group.1, fill = value)) +
  scale_fill_gradientn("CLR RA", name="Normalized\nExpression",
                       breaks = seq(from = min(i_t2$value), to = max(i_t2$value), by = 9000), 
                       colours = c("#e2ecd5", "#99be6c", "#6ea32e", "#2e6ea3"), 
                       limits = c(min(i_t2$value),max(i_t2$value))) +
  geom_text(data = i_t2[i_t2$abx < 0.05,], 
            aes(x = 1, y = Group.1, label = "*"), col = "white", size = 4) +
  labs(x= "Treatment", y="Reactome Pathway") +
  scale_x_discrete(labels=c("BMD" = "BMD", "No_Inoc" = "Control", "TJPbx"="T-Pbx", "FMB11"="FMB11", "GroGel"= "Prebx"))
pdf(paste(main_fp, "/Supplemental_Figure7.pdf", sep=""), width=12, height=25)
pathways_diff
dev.off()

######################################################################
#Are any pathways associated with bacteria?

keeps <- rownames(i_map[i_map$Collection_Day == "D06",])

keeps2 <- rownames(mapping[mapping$transcriptomic %in% keeps,])
keeps3 <- mapping$transcriptomic[mapping$transcriptomic %in% keeps]

path1 <- reactome_samples[keeps3,]

x <- t(as.matrix(t(CLR_otutable[keeps2,])))

#Within BMD test for correlation
#pathway_ABX <- path1[rownames(path1) %in% i_map$SampleID[i_map$Treatment == "TJP"],]
#b_otu_ABX <- x[rownames(x) %in% TJPbx,]

pathway_ABX <- path1
b_otu_ABX <- x

keep_b <- otu_table3[,rownames(x)]
keep_b <- keep_b[rowSums(keep_b)<200 & rowSums(keep_b) > 5,] #334 OTUs

b_otu_ABX <- b_otu_ABX[,colnames(b_otu_ABX) %in% rownames(keep_b)]

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

##This takes a while!
correlation.table <- associate(pathway_ABX, b_otu_ABX, method = "pearson", mode = "table", n.signif = 1)
b_pathways_correlations <- heat(correlation.table, "X1", "X2", 
                                fill = "Correlation", p.adj.threshold = 0.25, 
                                star='p.adj', star.size=3, 
                                colours=c("purple4", "#5E3C99", "white", "#f5bf99", "#e66101")) 


######################################################################
#Are any transcripts associated with bacteria?
View(C)
keeps <- rownames(i_map[i_map$Collection_Day == "D06",])

keeps2 <- rownames(mapping[mapping$transcriptomic %in% keeps,])
keeps3 <- mapping$transcriptomic[mapping$transcriptomic %in% keeps]

path1 <- t(ileum_transcript[,keeps3])

x <- t(as.matrix(t(CLR_otutable[keeps2,])))

#Within BMD test for correlation
#pathway_ABX <- path1[rownames(path1) %in% i_map$SampleID[i_map$Treatment == "TJP"],]
#b_otu_ABX <- x[rownames(x) %in% TJPbx,]

pathway_ABX <- path1
b_otu_ABX <- x

keep_b <- otu_table3[,rownames(x)]
keep_b <- keep_b[rowSums(keep_b)<200 & rowSums(keep_b) > 5,] #334 OTUs

b_otu_ABX <- b_otu_ABX[,colnames(b_otu_ABX) %in% rownames(keep_b)]

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

#This takes a while...
#correlation.table2 <- associate(pathway_ABX, b_otu_ABX, method = "pearson", mode = "table", n.signif = 5)
#b_trans_correlations <- heat(correlation.table2, "X1", "X2", fill = "Correlation", p.adj.threshold = 0.05, star='p.adj', star.size=1, colours= c("purple4", 
#                                                                                                                                                 "#5E3C99", "white", "#f5bf99", "#e66101")) 



######################################################################
#Plot Weight
weight_map <- mapping
weight_map$Treatment2 <- factor(weight_map$Treatment2, levels=c("BMD", "No_Inoc", "TJPbx", "FMB11","GroGel"))

colnames(weight_map)[7] <- 3
colnames(weight_map)[8] <- 6
colnames(weight_map)[10] <- 13
weight_map$SampleID <- rownames(weight_map)

weight_map2 <- melt(weight_map, id.vars = c('SampleID', 'Treatment2', "Gain313", "Gain613"), measure.vars = c('3','6','13'))
weight_map2$variable <- as.numeric(as.character(weight_map2$variable))
weight_map2$value <- as.numeric(as.character(weight_map2$value))

colnames(weight_map2)[5] <- "Day"
colnames(weight_map2)[6] <- "Weight"

time_w <- ggplot(weight_map2, aes(x=Day, y=Weight, color=Treatment2, group=Treatment2))+
  geom_jitter(width = 0.35, alpha=0.55) +
  geom_smooth(se=FALSE) +
  scale_color_manual(values=cols_pro) +
  scale_x_continuous(breaks=c(3,6,13)) +
  labs(y="Weight (g)")

#Compare the weights each day
w_map <- mapping[mapping$Treatment2 %in% c("BMD", "No_Inoc", "TJPbx", "FMB11","GroGel"),]
w_map$D3_weight <- as.numeric(w_map$D3_weight)
w_map$D6_weight <- as.numeric(w_map$D6_weight)
w_map$D13_weight <- as.numeric(w_map$D13_weight)
w_map$Treatment2 <- factor(w_map$Treatment2, levels=c("BMD", "No_Inoc", "TJPbx", "FMB11","GroGel"))
w_map2 <- w_map[,c("D3_weight", "D6_weight", "D13_weight", "Treatment2", "SampleID")]
w_map2 <- melt(w_map2, id.vars=c("Treatment2", "SampleID") )

weights <- ggplot(w_map2, aes(x=Treatment2, y=value)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width=0.15, alpha=0.75,aes(color=Treatment2)) +
  scale_color_manual(values=cols_pro) +
  facet_grid(variable~., scales="free") +
  guides(color=F)

#Stats for weights:
pair_weight_list <- list()
columns_weights <- c("D3_weight", "D6_weight", "D13_weight")
for(d in 1:length(columns_weights)){
  this_weight <- columns_weights[d]
  for(t in 1:4){
    for(u in (t+1):5){
      t1 <- as.character(unique(w_map$Treatment2)[[t]])
      t2 <- as.character(unique(w_map$Treatment2)[[u]])
      sams <- rownames(w_map[w_map$Treatment2 == t1 | w_map$Treatment2 == t2,])
      output = t.test(w_map[sams,this_weight] ~ w_map[sams,"Treatment2", drop=T])
      p_val <- output$p.value
      pair_weight_list[[paste(t1, "_", t2, "_", this_weight, sep="")]] <- c("pval", p_val)
    }
  }
}
sink(paste(main_fp, "treatment_pairwise_weights.txt", sep="/"), append =T)
print(pair_weight_list)
sink()

## Change in weight from day 3 to 13
weight_map2$Gain313 <- as.numeric(weight_map2$Gain313)
weight_map2$Gain613 <- as.numeric(weight_map2$Gain613)
weight_map2 <- weight_map2[weight_map2$Treatment2 %in% c("BMD", "No_Inoc", "TJPbx", "FMB11","GroGel"),]
weight <- ggplot(weight_map2, aes(x=Treatment2, y=Gain313)) +
  geom_boxplot() +
  geom_jitter(width=0.15, alpha=0.75, aes(color=Treatment2)) +
  scale_color_manual(values=cols_pro) +
  guides(color=F) +
  labs(y="Change in Weight Day 3 to 13 (g)", x="Treatment")

#Stats for weights:
pair_weight_list <- list()
w_map$Gain313 <- as.numeric(w_map$Gain313)
for(t in 1:4){
  for(u in (t+1):5){
    t1 <- as.character(unique(w_map$Treatment2)[[t]])
    t2 <- as.character(unique(w_map$Treatment2)[[u]])
    sams <- rownames(w_map[w_map$Treatment2 == t1 | w_map$Treatment2 == t2,])
    output = t.test(w_map[sams,"Gain313"] ~ w_map[sams,"Treatment2", drop=T])
    p_val <- output$p.value
    pair_weight_list[[paste(t1, "_", t2, "_Change in weight 3 to 13", sep="")]] <- c("pval", p_val)
  }
}

sink(paste(main_fp, "treatment_pairwise_weights.txt", sep="/"), append =T)
print(pair_weight_list)
sink()


######################################################################
#Compile Figures
path_p <- plot_grid(pathway_pcoa, NULL, ncol=1, rel_heights = c(1,0.5))
pathways_plots <- plot_grid(path_p, b_pathways_correlations, ncol=2)
pdf(paste(main_fp, "/Supplemental_Figure8.pdf", sep=""), height=8, width=8, useDingbats = F)
pathways_plots
dev.off()

#Plot weight over time and weight for each time
right_t <- plot_grid(time_w, NULL, nrow=2, rel_heights = c(1,2))
weight_t <- plot_grid(time_w, weights, ncol=2)
pdf(paste(main_fp, "/Supplemental_Figure9.pdf", sep=""), height=4, width=2, useDingbats = F)
weights
dev.off()

#plot procrustes and change in weight
side <- plot_grid(weight, NULL, nrow=2, rel_heights = c(2,1.3))
pro_path <- plot_grid(pro_b_p, side, ncol=2)
pdf(paste(main_fp, "/Figure5.pdf", sep=""), height=5, width=7, useDingbats = F)
pro_path
dev.off()

count_map <- mapping[!mapping$Gain313 == "N_A",]
nrow(count_map[count_map$Treatment2 == "BMD",])
nrow(count_map[count_map$Treatment2 == "TJPbx",])
nrow(count_map[count_map$Treatment2 == "No_Inoc",])
nrow(count_map[count_map$Treatment2 == "FMB11",])
nrow(count_map[count_map$Treatment2 == "GroGel",])
