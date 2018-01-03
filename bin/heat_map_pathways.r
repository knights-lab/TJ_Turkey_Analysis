#heat map for patways vs OTUs

#Make the names of the samples the same for Pathways and otu table
samples2 <- rownames(mapping[mapping$SampleID %in% samples & mapping$transcriptomic %in% colnames(ileum_transcript),])
samples3 <- mapping[samples2,"transcriptomic"]
x <- as.matrix(t(CLR_otutable[samples2,]))
rownames(x) <- taxonomy3$V2
x <- x[top_50_taxa,] #keep only the top 50 diff taxa (from abx vs no inoc)

# replace rownames with unique ones
otu_names <- as.character(rownames(x))
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
rownames(x) <- taxa$taxa

y <- t(reactome_samples[samples3,top_50_path]) 
colnames(y) <- colnames(x)

x <- as.data.frame(t(x))
y <- as.data.frame(t(y))

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

colnames(dat)<-c("OTUs","Pathways","Correlation","Pvalue")
dat$Correlation <- as.numeric(dat$Correlation)
dat$Pvalue <- as.numeric(dat$Pvalue)
dat$fdr_pval <- (p.adjust(dat$Pvalue, 'fdr'))

dat$Significance<-cut(dat$fdr_pval, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) #converts numeric p-values to factors based on significance level

# prune for complete cases with significant correlations
dat.c<-dat[complete.cases(dat),]
taxkeep<-unique(dat.c$OTUs[dat.c$fdr_pval<1])

taxassign<-taxkeep
dat.w<-dat.c[dat.c$OTUs %in% taxassign,] 

transkeep <- unique(dat.w$Pathways[dat.w$fdr_pval<1])
dat.m<-dat.w[dat.w$Pathways %in% transkeep,]

# Fix labeling
dat.m$Pathways <- gsub(".*L3_", "", dat.m$Pathways)
dat.m$Pathways <- gsub("_", " ", dat.m$Pathways)
dat.m$Pathways <- gsub(".*;s__", "", dat.m$Pathways)
dat.m$Pathways <- gsub("_", " ", dat.m$Pathways)


# order the OTUs
# make wide Pathways v. taxa filled with correlation
library(dplyr)
library(tidyr)
order <- dat.m %>% select(OTUs, Pathways, Correlation)
order <- order %>% spread(OTUs, Correlation)
rownames(order) <- order$Pathways
order <- order[,2:ncol(order)]
otu_order <- hclust((dist(1-cor(order))/2))$order

#can't reorder unless it's a factor first
dat.m$OTUs <- as.factor(dat.m$OTUs)

# reorder the factors in the dataframe to be in the taxaorder
dat.m$OTUs <- factor(dat.m$OTUs, levels(dat.m$OTUs)[otu_order])

# get the order Pathways
order <- dat.m %>% select(OTUs, Pathways, Correlation)
order <- order %>% spread(Pathways, Correlation)
rownames(order) <- order$OTUs
order <- order[,2:ncol(order)]
trans_order <- hclust((dist(1-cor(order))/2))$order

# set Pathways as factor
dat.m$Pathways <- as.factor(dat.m$Pathways)

# reorder the factors in the dataframe to be in the taxaorder
dat.m$Pathways <- factor(dat.m$Pathways, levels(dat.m$Pathways)[trans_order])

#plot the correlations for the collapsed levels
name <- "Pathway_Taxa_Heatmap.pdf"
fp <- paste(out_dir, name, sep="/")
pdf(fp, height=8,width=10)
ggplot(data = dat.m, aes(x=OTUs, y=Pathways, fill=Correlation)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low="#5e3c99", mid="#f7f7f7", high="#e66101", midpoint = 0, limit = c(-1,1), space = "Lab", name = "Spearman") +
  geom_text(data = dat.m, aes(label=Significance), color="black", size=2) +
  #theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1), 
        axis.text.y = element_text(size = 8)) +
  labs(x = "OTUs", y = "Pathways") +
  theme(strip.text.y = element_text(angle = 0, face = "italic"), 
        strip.background = element_rect(color="grey", fill = "white")) +
  coord_fixed()
dev.off()
