#Figure 1
#Succession of the turkey microbiome across body sites (days 3-13)

#Part 1
#Taxa summaries for control birds - one for each bodysite, by day
#Make taxa summaries
taxa_table1 <- taxa_table
tissues <- c(Bodysites[[2]], Bodysites[[3]])
taxa_table1 <- taxa_table1[,tissues]

#Make taxa_table
#Pick which Level and make names accordingly
level=6

names_split <- array(dim=c(nrow(taxa_table1), level))
otu_names <- as.character(rownames(taxa_table1))

for (i in 1:length(otu_names)){
  names_split[i,] <- head(strsplit(otu_names[i], ";", fixed=T)[[1]], n=level)
}

otu_names <- apply(names_split, 1, function(x) paste(x[1:level], sep = "", collapse = ";"))

taxa_table1$taxonomy <- otu_names
samples <- ncol(taxa_table1) - 1
otu <- aggregate(taxa_table1[,1:samples], by=list(taxa_table1$taxonomy), FUN = sum)
names(otu)[1] <- "taxonomy"

#Just for lowest level
holders <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
holder <- holders[level]
otu_names <- as.character(otu$taxonomy)
taxa <- c()
for (i in 1:length(otu_names)){
  taxon <- (strsplit(otu_names[i], holder, fixed=T))[[1]][2]
  if(is.na(taxon)){
    taxon <- (strsplit(otu_names[i], holder, fixed=T))[[1]][1]
    taxon <- (strsplit(taxon, "; "))[[1]][5]
  }
  if(taxon == "f__"){
    taxon <- (strsplit(otu_names[i], holder, fixed=T))[[1]][1]
    taxon <- (strsplit(taxon, "; "))[[1]][4]
  }
  taxa <- c(taxa, taxon)
}

#Add taxonomy as rownames
otu$taxonomy <- taxa
otu[is.na(otu)] <- "Other"
taxa_list <- otu$taxonomy
otu <- as.matrix(otu[,2:ncol(otu)])
rownames(otu) <- taxa_list
otu <- otu[!rownames(otu) == "Other",]

#Assign taxa colors
taxa_cols <- sample(cols2(length(rownames(otu))))
names(taxa_cols) <- rownames(otu)
taxa_cols <- c(taxa_cols, "#d3d3d3")
names(taxa_cols)[nrow(otu)+1] <- "Other"

#convert to RA
otu <- sweep(otu,2,colSums(otu),`/`)

#Make taxa plot for each bodysite on each day, facet by treatment
taxa_dir <- paste(main_fp, "taxa_sum/", sep='/')

cutoff <- 0.10
# Creat one plot for just no-inoc controls, one plot per day

new_cols <- c( "#1B9E77", "#361b9e", "#771b9e", "#E7298A", "#666666", "#E6AB02", "#c978a3", "#A6761D", "#731445", "#7570B3",  
               "#fad4e7", "#D95F02", "#1e66a6","#92a49e", 
               "#66A61E","#a5a2ce" , "#a3c978","#e3e2ef",
               "#b37570", "#d3d3d3")

sub_sam <- intersect(colnames(otu), NoInoc)
otu1 <- make_taxa_sums(otu, sub_sam, cutoff)
otu1$Tissue <- factor(otu1$Tissue, levels=c("Ileum", "Ceca"))
taxa_plot <- ggplot(otu1, aes(x = SampleID , y = Relative_Abundance)) + 
  geom_bar(stat="identity",aes(fill=Taxa)) +
  facet_grid(. ~ Tissue + Collection_Day, scales = "free_x") +
  theme_bw() +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "in"),
        legend.text = element_text(size=5),
        legend.position = 'bottom',
        axis.text = element_text(size=5),
        axis.text.x = element_blank(),
        axis.title = element_text(size=8)) + 
  scale_fill_manual(values= new_cols)

######################################################################
#PCoAs of each bodysite, colored by day
#Make PCOA plots

###########
otu_test <- read.table("data/Turkey_10k_Taxa.txt", 
                       sep="\t", 
                       comment="", 
                       header=T, 
                       skip=1, 
                       as.is=T, 
                       check.names=F, row=1)
tree = read_tree_greengenes('data/97_otus.tree')

beta_keeps <- intersect(rownames(mapping[NoInoc,]), colnames(otu_test))

cols3 <- c("#2C763B", "#ABB0DA", "#723181")

plot_list <- list()

for(i in 1:length(Bodysites)){
  body <- Bodysites[[i]]
  samples <- intersect(body, beta_keeps)
  
  #subset to keep samples
  m <- mapping[samples,]
  keeps <- intersect(rownames(m), colnames(otu_test))
  otu_test1 <- otu_test[,keeps]
  otu_test1 <- otu_test1[rowSums(otu_test1 > 0) > 1, ]
  m <- m[keeps,]
  #generate phyloseq objects
  OTU <- otu_table(otu_test1, taxa_are_rows=T)
  sampledata = sample_data(m)
  #manditory input for UniFrac
  physeq1 = phyloseq(OTU, sampledata, tree)
  #distances
  uni1 <- UniFrac(physeq1, weighted=FALSE, normalized=TRUE)
  
  #beta_table <-  beta_div[samples, samples]
  PCOA <- pcoa(uni1)$vectors
  var_exp <- pcoa(uni1)$values
  for(c in 1:ncol(PCOA)){
    colnames(PCOA)[c] <- paste("PC",c, sep="")
  }
  PCOA <- cbind(PCOA, rownames(PCOA))
  colnames(PCOA)[ncol(PCOA)] <- "SampleID"

  PCOA <- merge(PCOA, m, by="SampleID")
  PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
  PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
  PCOA <- as.data.frame(PCOA)

  #Run stats for diff. centroids
  beta_dist = as.dist(uni1)
  ad = adonis(beta_dist ~ m[,"Collection_Day"], data=m, permutations=999)
  p_val <- ad$aov.tab[1,6]
  r_sq <- ad$aov.tab[1,5]
  #Run Stats for diff. dispersion
  beta_out <- betadisper(beta_dist, m$Collection_Day)
  p_val_disp <- permutest(beta_out)$tab[1, 6]
  
  keeps <- intersect(rownames(otu_test1), colnames(CLR_otutable))
  x <- CLR_otutable[as.character(PCOA$SampleID),keeps]
  x <- x[,colSums(x) > 0.1]
  if(i == 2){
    correlation.table_i <- associate(x, PCOA[,2:3], method = "pearson", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
    correlation.table_i <- correlation.table_i[correlation.table_i$p.adj<0.05,]
    correlation.table_i <- correlation.table_i[correlation.table_i$Correlation<-0.6 | correlation.table_i$Correlation > 0.6,]
    correlation.table_i <- droplevels.data.frame(correlation.table_i[correlation.table_i$X2 == "PC1",])
    correlation.table_i$taxa <- taxonomy3[correlation.table_i$X1, "V2"]
    x <- data.frame(x[,correlation.table_i$X1])
    colnames(x) <- correlation.table_i$taxa
    x$SampleID <- rownames(x)
    x$Day <- PCOA$Collection_Day
    x <- melt(x)
    x$Day[x$Day == "D03"] <- 3
    x$Day[x$Day == "D06"] <- 6
    x$Day[x$Day == "D13"] <- 13
    x$variable <- as.character(x$variable)
    x$variable <- factor(x$variable, levels = unique(as.character(correlation.table_i$taxa)))
    x_i <- x
  }
  if(i == 3){
    correlation.table_c <- associate(x, PCOA[,2:3], method = "pearson", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
    correlation.table_c <- correlation.table_c[correlation.table_c$p.adj<0.05,]
    correlation.table_c <- droplevels.data.frame(correlation.table_c[correlation.table_c$X2 == "PC1",])
    correlation.table_c$taxa <- taxonomy3[correlation.table_c$X1, "V2"]
    x <- data.frame(x[,correlation.table_c$X1])
    colnames(x) <- correlation.table_c$taxa
    x$SampleID <- rownames(x)
    x$Day <- PCOA$Collection_Day
    x <- melt(x)
    x$Day[x$Day == "D03"] <- 3
    x$Day[x$Day == "D06"] <- 6
    x$Day[x$Day == "D13"] <- 13
    x$variable <- as.character(x$variable)
    x$variable <- factor(x$variable, levels = unique(as.character(correlation.table_c$taxa)))
    x_c <- x
  }
  plot1 <- ggplot(PCOA) +
    geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Collection_Day"), alpha =0.8) + 
    scale_color_manual(values= cols3) +
    annotate("text", x=0.2, y=0.2, label= paste("P=", p_val), size=2) +
    annotate("text", x=0.2, y=0.25, label= paste("R2=", round(r_sq, digits=3)), size=2) +
    labs(x=paste("PC1 (", round(var_exp$Relative_eig[1], digits=3)*100, "%)", sep=""), y= paste("PC2 (", round(var_exp$Relative_eig[2], digits=3)*100, "%)", sep=""))
  plot_list[[names(Bodysites)[i]]] <- plot1
}

######################################################################
#Diff taxa within tissue but by day / time
#Make a table with taxonomy as the rownames
taxonomy3 <- taxonomy[colnames(CLR_otutable),]
working_table <- as.matrix(t(CLR_otutable))
#rownames(working_table) <- taxonomy3$V2

#Set up tests to run
test.ixs <- list(Days[[1]], Days[[2]], Days[[3]])
names(test.ixs) <- c("Day3", "Day6", "Day13")

pvals<- c()
ALPHA <- 0.05
diff_list <- list()

for(i in 1:length(Bodysites)){
  Bodysite <- names(Bodysites)[i]
  union1 <- intersect(Bodysites[[i]], NoInoc)
    for(n in 1:(length(test.ixs)-1)){
      #for(m in (n+1):length(test.ixs)){
        test.x <- test.ixs[n]
        test.y <- test.ixs[n+1]
        set1 <- intersect(union1, test.x[[1]])
        set1 <- intersect(set1, colnames(working_table))
        set2 <- intersect(union1, test.y[[1]])
        set2 <- intersect(set2, colnames(working_table))
        if(length(set1) > 2 && length(set2) > 2){
          full_set <- c(set1, set2)
          
          keep_t <- otu_table3[,full_set]
          keep_t <- otu_table3[rowSums(keep_t) > 5,]
          test_table <- working_table[rownames(keep_t), full_set]
          
          rownames(test_table) <- taxonomy3[rownames(test_table), 2]
          
          #keep taxa and the samples you are testing
          test_table <- t(test_table)
        
          #Keep taxa that have a mean rel abundance of at least 0.01
          #test_table <- test_table[,colMeans(test_table)>0.009, drop=F]
          map_test <- mapping[rownames(test_table),]
          difftest <- test.otu.features(test_table, response=map_test$Collection_Day, sig.level = 0.10)
          #difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
          
          if(any(difftest$pvals <= ALPHA)){
            signif.ix <- which(difftest$pvals <= ALPHA)
            signif.ix <- signif.ix[order(difftest$pvals[signif.ix])]
            sink("list of diff taxa.txt", append=T)
            cat(paste(Bodysite, "\n"))
            cat(paste(names(test.ixs)[n], " vs ", names(test.ixs)[n+1], "\n"))
            print(table(names(signif.ix)))
            cat("\n")
            sink()
            diff_list[[paste(Bodysite, names(test.ixs)[n])]] <- names(signif.ix)
          } else {
          cat("Not Sign.")
        }
      }
    }
}
rownames(working_table) <- taxonomy3[rownames(working_table),2] 
keep_taxa <- unique(c(diff_list[["Ileum Day3"]], diff_list[["Ileum Day6"]]))
x <- data.frame(working_table[keep_taxa,intersect(NoInoc, Ileum)])
x$day3 <- 1
x$day6 <- 1
x[rownames(x) %in% diff_list[["Ileum Day3"]], "day3"] <- 0.04
x[rownames(x) %in% diff_list[["Ileum Day6"]], "day6"] <- 0.04

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
i_table <- data.frame(t(x))
ps <- i_table[c("day3", "day6"),]
i_table <- i_table[!rownames(i_table) %in% c("day3", "day6"),]
i_map <- mapping[rownames(i_table),]
i_table$Day <- i_map$Collection_Day
i_table <- melt(i_table)
i_t2 <-aggregate(i_table, by=list(i_table$Day, i_table$variable), 
                            FUN=mean, na.rm=TRUE)
ps <- data.frame(t(ps))
ps$Group.2 <- rownames(ps)
i_t2 <- merge(i_t2, ps, by="Group.2")
i_t2$Group.1  <- factor(i_t2$Group.1, levels=c("D03", "D06", "D13"))
taxa_order <- i_t2[i_t2$Group.1 == "D03",]
taxa_order <- taxa_order[order(-taxa_order$value),]
i_t2$Group.2 <- factor(i_t2$Group.2, levels=taxa_order$Group.2)
i_t2$day3 <- as.numeric(as.character(i_t2$day3))
i_t2$day6 <- as.numeric(as.character(i_t2$day6))
ileum_days <- ggplot(i_t2) + 
  geom_tile(aes(x = Group.1, y = Group.2, fill = value)) +
  scale_fill_gradientn("CLR RA", 
                              breaks = seq(from = -1, to = 8, by = 1), 
                              colours = c("#FEBF2C",  "#C3003A",  "#551B44"), 
                              limits = c(-1,8)) +
  geom_text(data = i_t2[i_t2$day3 < 0.05,], 
              aes(x = 1, y = Group.2, label = "*"), col = "white", size = 4) +
  geom_text(data = i_t2[i_t2$day6 < 0.05,], 
            aes(x = 3, y = Group.2, label = "*"), col = "white", size = 4)



keep_taxa <- unique(c(diff_list[["Cecum Day3"]], diff_list[["Cecum Day6"]]))
x <- data.frame(working_table[keep_taxa,intersect(NoInoc, Ceca)])
x$day3 <- 1
x$day6 <- 1
x[rownames(x) %in% diff_list[["Cecum Day3"]], "day3"] <- 0.04
x[rownames(x) %in% diff_list[["Cecum Day6"]], "day6"] <- 0.04

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
c_table <- data.frame(t(x))
ps <- c_table[c("day3", "day6"),]
c_table <- c_table[!rownames(c_table) %in% c("day3", "day6"),]

c_map <- mapping[rownames(c_table),]
c_table$Day <- c_map$Collection_Day


c_table <- melt(c_table)
c_t2 <-aggregate(c_table, by=list(c_table$Day, c_table$variable), 
                 FUN=mean, na.rm=TRUE)
ps <- data.frame(t(ps))
ps$Group.2 <- rownames(ps)
c_t2 <- merge(c_t2, ps, by="Group.2")

c_t2$Group.1  <- factor(c_t2$Group.1, levels=c("D03", "D06", "D13"))
taxa_order <- c_t2[c_t2$Group.1 == "D03",]
taxa_order <- taxa_order[order(-taxa_order$value),]
c_t2$Group.2 <- factor(c_t2$Group.2, levels=taxa_order$Group.2)
c_t2$day3 <- as.numeric(as.character(c_t2$day3))
c_t2$day6 <- as.numeric(as.character(c_t2$day6))

ceca_days <-ggplot(c_t2) + 
  geom_tile(aes(x = Group.1, y = Group.2, fill = value)) +
  scale_fill_gradientn("CLR RA", 
                       breaks = seq(from = -2, to = 7, by = 1), 
                       colours = c("#FEBF2C",  "#C3003A",  "#551B44"), 
                       limits = c(-2,7)) +
  geom_text(data = c_t2[c_t2$day3 < 0.05,], 
            aes(x = 1, y = Group.2, label = "*"), col = "white", size = 4) +
  geom_text(data = c_t2[c_t2$day6 < 0.05,], 
            aes(x = 3, y = Group.2, label = "*"), col = "white", size = 4)

######################################################################

heats <- plot_grid(ileum_days, ceca_days, ncol=2, rel_widths = c(1, 0.95))
pcoas <- plot_grid(plot_list[[2]], plot_list[[3]], nrow=2)
bottom <- plot_grid(pcoas, heats, ncol=2, rel_widths = c(0.40, 1))

pdf(paste(main_fp, "/Figure1_bottom.pdf", sep=""), height=4, width=12, useDingbats =F)
print(bottom)
dev.off()

pdf(paste(main_fp, "/Figure1_top.pdf", sep=""), height=4, width=12, useDingbats=F)
print(taxa_plot)
dev.off()


######################################################################
#Supplemental Figure 1

#alpha diversity by time and bodysite
alpha_intersect <- intersect(rownames(alpha_div), rownames(mapping))

alpha_div <- alpha_div[alpha_intersect,]
a_mapping <- mapping[alpha_intersect,]
a_mapping$shannon <- alpha_div[,"shannon"]
a_mapping$simpson <- alpha_div[,"simpson"]
a_mapping$pdwholetree <- alpha_div[,"PD_whole_tree"]
a_mapping$obsspecies <- alpha_div[,"observed_species"]

a_mapping <- a_mapping[a_mapping$Treatment2 == "No_Inoc" & a_mapping$Tissue %in% c("Ileum", "Ceca"),]

#run stats
pval_list <- list()
bodies <- Bodysites[2:3]
for(b in 1:length(bodies)){
  body <- bodies[[b]]
  for(d in 1:(length(Days)-1)){
    for(e in (d+1):length(Days)){
      day1 <- Days[[d]]
      day2 <- Days[[e]]
      samples1 <- intersect(body,day1)
      samples2 <- intersect(body, day2)
      samples1 <- intersect(samples1, rownames(a_mapping))
      samples2 <- intersect(samples2, rownames(a_mapping))
      a_test <- t.test(a_mapping[samples1, "pdwholetree"], a_mapping[samples2, "pdwholetree"])$p.value
      pval_list[[paste(names(bodies)[b], names(Days)[d], names(Days)[e], sep="_")]] <- a_test
    }
  }
}
sink(paste(main_fp, "s_Fig1_alpha_stats.txt", sep="/"))
pval_list
sink()

#run stats as linear correlation 
pval_list <- list()
bodies <- Bodysites[2:3]
a_mapping$Day <- a_mapping$Collection_Day
a_mapping$Day <- gsub("D", "", a_mapping$Day)
a_mapping$Day <- as.numeric(a_mapping$Day)
for(b in 1:length(bodies)){
  body <- bodies[[b]]
  samples1 <- intersect(body, rownames(a_mapping))
  l_test <- summary(lm(a_mapping[samples1, "pdwholetree"] ~ a_mapping[samples1, "Day"]))
  pval_list[[paste(names(bodies)[b], "Time", sep="_")]] <- l_test
}
sink(paste(main_fp, "s_Fig1_alpha_stats.txt", sep="/"), append=T)
pval_list
sink()

#Stats comparing body sites by day
pval_list <- list()
bodies <- Bodysites[2:3]
for(d in 1:length(Days)){
  day1 <- Days[[d]]
  body1 <- bodies[[1]]
  body2 <- bodies[[2]]
  samples1 <- intersect(body1,day1)
  samples2 <- intersect(body2, day1)
  samples1 <- intersect(samples1, rownames(a_mapping))
  samples2 <- intersect(samples2, rownames(a_mapping))
  a_test <- t.test(a_mapping[samples1, "pdwholetree"], a_mapping[samples2, "pdwholetree"])$p.value
  pval_list[[paste(names(bodies)[1], names(bodies)[2], names(Days)[d], sep="_")]] <- a_test
}
sink(paste(main_fp, "s_Fig1_alpha_stats_body.txt", sep="/"))
pval_list
sink()

a_mapping$Tissue <- factor(a_mapping$Tissue, levels=c("Ileum", "Ceca"))
alpha_div_plot <- ggplot(a_mapping, aes(x=Day, y=pdwholetree))+
  geom_jitter(alpha=0.8, width=.15, aes(group=Day, color=Collection_Day)) +
  geom_boxplot(fill=NA, aes(group=Day)) +
  #facet_grid(Tissue ~ ., switch="both") +
  facet_wrap(~ Tissue, strip.position="top", ncol=1) +
  scale_color_manual(values = cols3) +
  guides(color=F) +
  geom_smooth(color=cols3[2], se=F) +
  scale_x_continuous(breaks=c(3,6,13))


#beta diversity by three ways:
bray <- read.table(bray_fp, sep="\t", header=T, row=1)
uni <- read.table(uni_fp, sep="\t", header=T, row=1)
wuni <- read.table(wuni_fp, sep="\t", header=T, row=1)

beta_keeps <- intersect(rownames(mapping), rownames(bray))
mapping2 <- mapping[beta_keeps,]

beta_tables <- list(uni, wuni, bray)
beta_names <- c("unifrac", "weighted_unifrac", "bray_curtis")

bodies <- Bodysites[2:3]
plot_list <- list()
for(b in 1:length(beta_tables)){
  beta_div <- as.data.frame(beta_tables[b])
  for(i in 1:length(bodies)){
    body <- bodies[[i]]
    samples <- intersect(body, rownames(beta_div))
    samples <- intersect(samples, rownames(mapping[mapping$Treatment2 == "No_Inoc",]))
    beta_table <-  beta_div[samples, samples]
    
    #Run stats for diff. centroids
    beta_dist = as.dist(beta_table)
    m <- mapping[samples,]
    ad = adonis(beta_dist ~ m[,"Collection_Day"], data=m, permutations=999)
    p_val <- ad$aov.tab[1,6]
    r_sq <- ad$aov.tab[1,5]
    #Run Stats for diff. dispersion
    beta_out <- betadisper(beta_dist, m$Collection_Day)
    p_val_disp <- permutest(beta_out)$tab[1, 6]
    
    
    PCOA <- pcoa(beta_table)$vectors
    var_exp <- pcoa(beta_table)$values
    
    for(c in 1:ncol(PCOA)){
      colnames(PCOA)[c] <- paste("PC",c, sep="")
    }
    PCOA <- cbind(PCOA, rownames(PCOA))
    colnames(PCOA)[ncol(PCOA)] <- "SampleID"
    mapping2 <- mapping[samples,]
    PCOA <- merge(PCOA, mapping2, by="SampleID")
    PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
    PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
    PCOA <- as.data.frame(PCOA)
    
    plot1 <- ggplot(PCOA) +
      geom_point(size = 2, aes_string(x = "PC1", y = "PC2", color = "Collection_Day"), alpha=0.7) + 
      scale_color_manual(values= cols3) +
      annotate("text", x=0.2, y=-0.2, label= paste("P=", p_val), size=2) +
      annotate("text", x=0.2, y=-0.22, label= paste("R2=", round(r_sq, digits=3)), size=2) +
      labs(x=paste("PC1 (", round(var_exp$Relative_eig[1], digits=3)*100, "%)", sep=""), y= paste("PC2 (", round(var_exp$Relative_eig[2], digits=3)*100, "%)", sep=""))
    
    plot_list[[paste(names(bodies[i]), beta_names[b], sep="_")]] <- plot1
  }
}
ileum_all <- plot_grid(plot_list[["Ileum_unifrac"]],plot_list[["Ileum_weighted_unifrac"]],plot_list[["Ileum_bray_curtis"]], 
                       nrow=2, ncol=2, labels= c("unifrac", "weighted unifrac", "bray curtis"))
ceca_all <- plot_grid(plot_list[["Cecum_unifrac"]],plot_list[["Cecum_weighted_unifrac"]],plot_list[["Cecum_bray_curtis"]], 
                       nrow=2, ncol=2, labels= c("unifrac", "weighted unifrac", "bray curtis"))

together2 <- plot_grid(alpha_div_plot, NULL, nrow=2, rel_heights = c(1,0.8))
t3 <- plot_grid(alpha_div_plot, ileum_all, rel_widths = c(0.4, 1))
t4 <- plot_grid(ceca_all, NULL, rel_widths=c(1, 0.4))
t5 <- plot_grid(t3, t4, nrow=2)

plot_this <- paste(main_fp, "Supplemental_Figure1.pdf", sep='/')
pdf(plot_this, width=8, height=6, useDingbats = F)
t5
dev.off()
