#impact of abx on the ileum and cecum

#Ileum
beta <- as.matrix(vegdist(CLR_otutable, method="euclidean"))
abx_no <- union(NoInoc, BMD)
m <- mapping[mapping$Tissue == "Ileum",]

uni_fp <- paste(data_dir, "beta_diversity10k/unifrac_Turkey_10k_Taxa.txt", sep='')
uni <- read.table(uni_fp, sep="\t", head=T, row=1)

beta_keeps <- intersect(rownames(m[abx_no,]), rownames(uni))
mapping2 <- m[beta_keeps,]
beta_div <- data.frame(uni)
beta_div <- beta_div[beta_keeps, beta_keeps]
abx_col <- c("#1B9E77", "#D95F02")
mapping2$Treatment2 <- factor(mapping2$Treatment2)
# "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"
plot_list <- list()

for(i in 1:length(Days)){
  day <- Days[[i]]
  samples <- intersect(day, rownames(beta_div))
  beta_table <-  beta_div[samples, samples]
  
  PCOA <- pcoa(beta_table)$vectors
  var_exp <- pcoa(beta_table)$values
  
  #Run stats for diff. centroids
  beta_dist = as.dist(beta_table)
  ad = adonis(beta_dist ~ m[samples,"Treatment2"], data=m, permutations=999)
  p_val <- ad$aov.tab[1,6]
  r_sq <- ad$aov.tab[1,5]
  #Run Stats for diff. dispersion
  beta_out <- betadisper(beta_dist, m[samples, "Treatment2"])
  p_val_disp <- permutest(beta_out)$tab[1, 6]
  
  for(c in 1:ncol(PCOA)){
    colnames(PCOA)[c] <- paste("PC",c, sep="")
  }
  PCOA <- cbind(PCOA, rownames(PCOA))
  colnames(PCOA)[ncol(PCOA)] <- "SampleID"
  mapping2 <- mapping[samples,]
  PCOA <- merge(PCOA, mapping2, by="SampleID")
  PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
  PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
  PCOA$PC3 <- as.numeric(levels(PCOA$PC3))[PCOA$PC3]
  PCOA$PC4 <- as.numeric(levels(PCOA$PC4))[PCOA$PC4]
  PCOA <- as.data.frame(PCOA)
  PCOA$Treatment2 <- factor(PCOA$Treatment2, levels = c("No_Inoc", "BMD"))
  plot1 <- ggplot(PCOA) +
    geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Treatment2"), alpha =0.8) + 
    scale_color_manual(values= abx_col) +
    stat_ellipse(aes_string(x = "PC1", y = "PC2", color = "Treatment2"), linetype = 2) +
    annotate("text", x=0.2, y=0.2, label= paste("P=", p_val), size=2) +
    annotate("text", x=0.2, y=0.22, label= paste("R2=", round(r_sq, digits=3)), size=2) +
    labs(x=paste("PC1 (", round(var_exp$Relative_eig[1], digits=3)*100, "%)", sep=""), y= paste("PC2 (", round(var_exp$Relative_eig[2], digits=3)*100, "%)", sep=""))
  
  plot_list[[names(Days)[i]]] <- plot1
}
ileum_a_pcoa <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol=2)

######################################################################
#Add alpha diversity to map, and it will be way easier to plot
alpha_intersect <- intersect(rownames(alpha_div), rownames(mapping))

alpha_div <- alpha_div[alpha_intersect,]
a_mapping <- mapping[alpha_intersect,]
a_mapping$shannon <- alpha_div[,"shannon"]
a_mapping$simpson <- alpha_div[,"simpson"]
a_mapping$pdwholetree <- alpha_div[,"PD_whole_tree"]
a_mapping$obsspecies <- alpha_div[,"observed_species"]

#Make plot of alpha diversity
test_a <- alpha_div[beta_keeps,]
test_a$Collection_Day <- m[beta_keeps,"Collection_Day"]
test_a$Treatment <- m[beta_keeps,"Treatment2"]
test_a$Treatment <- factor(test_a$Treatment, levels = c("No_Inoc", "BMD"))
ileum_a_alpha <- ggplot(test_a) +
  geom_boxplot(aes(x=Treatment, y=shannon), outlier.colour = NA) +
  geom_jitter(aes(x=Treatment, y=shannon, color=Treatment), width=0.1, size=2, alpha=0.8) +
  facet_grid(. ~ Collection_Day) +
  scale_color_manual(values = abx_col) +
  guides(color=F)


#run stats
pval_list <- list()
bodies <- Bodysites[2:3]
for(b in 1:length(bodies)){
  body <- bodies[[b]]
  for(d in 1:length(Days)){
    day1 <- Days[[d]]
    samples <- intersect(body,day1)
    samples1 <- intersect(samples, rownames(a_mapping[a_mapping$Treatment2 %in% c("BMD"),]))
    samples2 <- intersect(samples, rownames(a_mapping[a_mapping$Treatment2 %in% c("No_Inoc"),]))
    a_test <- t.test(a_mapping[samples1, "shannon"], a_mapping[samples2, "shannon"])$p.value
    pval_list[[paste(names(bodies)[b], names(Days)[d], sep="_")]] <- a_test
  }
}
sink(paste(main_fp, "Fig2_alpha_stats.txt", sep="/"))
pval_list
sink()

######################################################################
#Find taxa that are different between the treatments
taxonomy3 <- taxonomy[colnames(CLR_otutable),]
working_table <- as.matrix(t(CLR_otutable))

#Set up tests to run
test.ixs <- list(NoInoc, BMD)
names(test.ixs) <- c("NoInoc", "BMD")

pvals<- c()
ALPHA <- 0.05
diff_list <- list()

for(i in 1:length(Days)){
  Day <- names(Days)[i]
  union1 <- intersect(Days[[i]], Ileum)
  for(n in 1:(length(test.ixs)-1)){
    for(m in (n+1):length(test.ixs)){
    test.x <- test.ixs[n]
    test.y <- test.ixs[m]
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
      map_test <- mapping[rownames(test_table),]
      difftest <- test.otu.features(test_table, response=map_test$Treatment, sig.level = 0.10)
      #difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
      
      if(any(difftest$pvals <= ALPHA)){
        signif.ix <- which(difftest$pvals <= ALPHA)
        signif.ix <- signif.ix[order(difftest$pvals[signif.ix])]
        sink("list of diff taxa ileum.txt", append=T)
        cat(paste(Day, "\n"))
        cat(paste(names(test.ixs)[n], " vs ", names(test.ixs)[m], "\n"))
        print(table(names(signif.ix)))
        cat("\n")
        sink()
        diff_list[[Day]] <- names(signif.ix)
      } 
    } 
    }
  }
}
rownames(working_table) <- taxonomy3[rownames(working_table),2] 
keep_taxa <- unique(c(diff_list[["D03"]], diff_list[["D06"]], diff_list[["D13"]]))
x <- data.frame(working_table[keep_taxa,intersect(abx_no, Ileum)])
x$day3 <- 1
x$day6 <- 1
x$day13 <- 1
x[rownames(x) %in% diff_list[["D03"]], "day3"] <- 0.04
x[rownames(x) %in% diff_list[["D06"]], "day6"] <- 0.04
x[rownames(x) %in% diff_list[["D13"]], "day13"] <- 0.04

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
ps <- i_table[c("day3", "day6", "day13"),]
i_table <- i_table[!rownames(i_table) %in% c("day3", "day6", "day13"),]
i_map <- mapping[rownames(i_table),]
i_table$Day <- i_map$Collection_Day
i_table$Treatment <- i_map$Treatment2
i_table <- melt(i_table)
i_t2 <-aggregate(i_table, by=list(i_table$Day, i_table$variable, i_table$Treatment), 
                 FUN=mean, na.rm=TRUE)
ps <- data.frame(t(ps))
ps$Group.2 <- rownames(ps)
i_t2 <- merge(i_t2, ps, by="Group.2")
i_t2$Group.1  <- factor(i_t2$Group.1, levels=c("D03", "D06", "D13"))
taxa_order <- i_t2[i_t2$Group.1 == "D03" & i_t2$Group.3 == "No_Inoc",]
taxa_order <- taxa_order[order(-taxa_order$value),]
i_t2$Group.2 <- factor(i_t2$Group.2, levels=taxa_order$Group.2)
i_t2$day3 <- as.numeric(as.character(i_t2$day3))
i_t2$day6 <- as.numeric(as.character(i_t2$day6))
i_t2$day13 <- as.numeric(as.character(i_t2$day13))
ileum_abx <- ggplot(i_t2) + 
  geom_tile(aes(x = Group.3, y = Group.2, fill = value)) +
  scale_fill_gradientn("CLR RA", 
                       breaks = seq(from = min(i_t2$value), to = max(i_t2$value), by = 1), 
                       colours = c("#FEBF2C",  "#C3003A",  "#551B44"), 
                       limits = c(min(i_t2$value),max(i_t2$value))) +
  facet_grid(.~ Group.1) +
  geom_text(data = i_t2[i_t2$day3 < 0.05 & i_t2$Group.1 == "D03",], 
            aes(x = 1 , y = Group.2, label = "*"), col = "white", size = 4) +
  geom_text(data = i_t2[i_t2$day6 < 0.05 & i_t2$Group.1 == "D06",], 
            aes(x = 1 , y = Group.2, label = "*"), col = "white", size = 4) +
  geom_text(data = i_t2[i_t2$day13 < 0.05 & i_t2$Group.1 == "D13",], 
            aes(x = 1 , y = Group.2, label = "*"), col = "white", size = 4) 

taxa_all_sig <- i_t2[i_t2$day3 ==0.04 &i_t2$day6 ==0.04 &i_t2$day13 ==0.04, ]
taxa_all_sig <- as.character(unique(taxa_all_sig$Group.2))

sink(paste(main_fp, "ileum_taxa_sig_all_times.txt", sep="/"))
taxa_all_sig
sink()

ileum_abx <- ggplot(i_t2) + 
  geom_tile(aes(x = Group.2, y = Group.3, fill = value)) +
  scale_fill_gradientn("CLR RA", 
                       breaks = seq(from = min(i_t2$value), to = max(i_t2$value), by = 1), 
                       colours = c("#FEBF2C",  "#C3003A",  "#551B44"), 
                       limits = c(min(i_t2$value),max(i_t2$value))) +
  facet_grid(Group.1 ~.) +
  geom_text(data = i_t2[i_t2$day3 < 0.05 & i_t2$Group.1 == "D03",], 
            aes(x = Group.2 , y = 1, label = "*"), col = "white", size = 4) +
  geom_text(data = i_t2[i_t2$day6 < 0.05 & i_t2$Group.1 == "D06",], 
            aes(x = Group.2 , y = 1, label = "*"), col = "white", size = 4) +
  geom_text(data = i_t2[i_t2$day13 < 0.05 & i_t2$Group.1 == "D13",], 
            aes(x = Group.2 , y = 1, label = "*"), col = "white", size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), axis.title.y = element_blank(),
        axis.title.x = element_blank())

topside <- plot_grid(NULL, ileum_a_alpha, ncol=2, rel_widths = c(0.2, 1))
side1 <- plot_grid(topside, ileum_abx, nrow=2, rel_heights = c(1,1))
top <- plot_grid(ileum_a_pcoa, side1, ncol=2, rel_widths = c(1.8,1.75))


pdf(paste(main_fp, "/Figure2_top.pdf", sep=""), height=4, width=13, useDingbats =F)
print(top)
dev.off()


######################################################################
#Ceca
beta <- as.matrix(vegdist(CLR_otutable, method="euclidean"))
abx_no <- union(NoInoc, BMD)
m <- mapping[mapping$Tissue == "Ceca",]
beta_keeps <- intersect(rownames(m[abx_no,]), rownames(uni))
mapping2 <- m[beta_keeps,]
beta_div <- data.frame(uni)
beta_div <- beta_div[beta_keeps, beta_keeps]
abx_col <- c("#1B9E77", "#D95F02")

plot_list <- list()

for(i in 1:length(Days)){
  day <- Days[[i]]
  samples <- intersect(day, rownames(beta_div))
  beta_table <-  beta_div[samples, samples]
  
  PCOA <- pcoa(beta_table)$vectors
  var_exp <- pcoa(beta_table)$values
  
  #Run stats for diff. centroids
  beta_dist = as.dist(beta_table)
  ad = adonis(beta_dist ~ m[samples,"Treatment2"], data=m, permutations=999)
  p_val <- ad$aov.tab[1,6]
  r_sq <- ad$aov.tab[1,5]
  #Run Stats for diff. dispersion
  beta_out <- betadisper(beta_dist, m[samples, "Treatment2"])
  p_val_disp <- permutest(beta_out)$tab[1, 6]
  
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
  PCOA$Treatment2 <- factor(PCOA$Treatment2, levels = c("No_Inoc", "BMD"))
  plot1 <- ggplot(PCOA) +
    geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Treatment2"), alpha =0.8) + 
    scale_color_manual(values= abx_col) +
    stat_ellipse(aes_string(x = "PC1", y = "PC2", color = "Treatment2"), linetype = 2) +
    annotate("text", x=0.2, y=0.2, label= paste("P=", p_val), size=2) +
    annotate("text", x=0.2, y=0.22, label= paste("R2=", round(r_sq, digits=3)), size=2) +
    labs(x=paste("PC1 (", round(var_exp$Relative_eig[1], digits=3)*100, "%)", sep=""), y= paste("PC2 (", round(var_exp$Relative_eig[2], digits=3)*100, "%)", sep=""))
  
  plot_list[[names(Days)[i]]] <- plot1
}
ceca_a_pcoa <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol=2)


######################################################################

#make plot of alpha diversity
test_a <- alpha_div[beta_keeps,]
test_a$Collection_Day <- m[beta_keeps,"Collection_Day"]
test_a$Treatment <- m[beta_keeps,"Treatment2"]
test_a$Treatment <- factor(test_a$Treatment, levels = c("No_Inoc", "BMD"))
ceca_a_alpha <- ggplot(test_a) +
  geom_boxplot(aes(x=Treatment, y=shannon), outlier.colour = NA) +
  geom_jitter(aes(x=Treatment, y=shannon, color=Treatment), width=0.1, size=2, alpha=0.8) +
  facet_grid(. ~ Collection_Day) +
  scale_color_manual(values = abx_col) +
  guides(color=F)


######################################################################
#Find taxa that are different between the treatments
taxonomy3 <- taxonomy[colnames(CLR_otutable),]
working_table <- as.matrix(t(CLR_otutable))

#Set up tests to run
test.ixs <- list(NoInoc, BMD)
names(test.ixs) <- c("NoInoc", "BMD")

pvals<- c()
ALPHA <- 0.05
diff_list <- list()

for(i in 1:length(Days)){
  Day <- names(Days)[i]
  union1 <- intersect(Days[[i]], Ceca)
  for(n in 1:(length(test.ixs)-1)){
    for(m in (n+1):length(test.ixs)){
      test.x <- test.ixs[n]
      test.y <- test.ixs[m]
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
        map_test <- mapping[rownames(test_table),]
        difftest <- test.otu.features(test_table, response=map_test$Treatment, sig.level = 0.10)
        #difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
        
        if(any(difftest$pvals <= ALPHA)){
          signif.ix <- which(difftest$pvals <= ALPHA)
          signif.ix <- signif.ix[order(difftest$pvals[signif.ix])]
          sink("list of diff taxa ileum.txt", append=T)
          cat(paste(Day, "\n"))
          cat(paste(names(test.ixs)[n], " vs ", names(test.ixs)[m], "\n"))
          print(table(names(signif.ix)))
          cat("\n")
          sink()
          diff_list[[Day]] <- names(signif.ix)
        } 
      } 
    }
  }
}
rownames(working_table) <- taxonomy3[rownames(working_table),2] 

keep_taxa <- unique(c(diff_list[["D03"]], diff_list[["D06"]], diff_list[["D13"]]))
x <- data.frame(working_table[keep_taxa,intersect(abx_no, Ceca), drop=F])
x_total <- x
x$day3 <- 1
x$day6 <- 1
x$day13 <- 1
x[rownames(x) %in% diff_list[["D03"]], "day3"] <- 0.04
x[rownames(x) %in% diff_list[["D06"]], "day6"] <- 0.04
x[rownames(x) %in% diff_list[["D13"]], "day13"] <- 0.04

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
ps <- c_table[c("day3", "day6", "day13"),,drop=F]
c_table <- c_table[!rownames(c_table) %in% c("day3", "day6", "day13"),,drop=F]
c_map <- mapping[rownames(c_table),]
c_table$Day <- c_map$Collection_Day
c_table$Treatment <- c_map$Treatment2
c_table <- melt(c_table)
c_t2 <-aggregate(c_table, by=list(c_table$Day, c_table$variable, c_table$Treatment), 
                 FUN=mean, na.rm=TRUE)
ps <- data.frame(t(ps))
ps$Group.2 <- rownames(ps)
c_t2 <- merge(c_t2, ps, by="Group.2")
c_t2$Group.1  <- factor(c_t2$Group.1, levels=c("D03", "D06", "D13"))
taxa_order <- c_t2[c_t2$Group.1 == "D03" & c_t2$Group.3 == "No_Inoc",]
taxa_order <- taxa_order[order(-taxa_order$value),]
c_t2$Group.2 <- factor(c_t2$Group.2, levels=taxa_order$Group.2)
c_t2$day3 <- as.numeric(as.character(c_t2$day3))
c_t2$day6 <- as.numeric(as.character(c_t2$day6))
c_t2$day13 <- as.numeric(as.character(c_t2$day13))

x$Day <- mapping[rownames(x), "Collection_Day"]
x$Treatment <- mapping[rownames(x), "Treatment2"]
colnames(x)[1] <- "CLR_Abundance"
x$Treatment <- factor(x$Treatment, levels = c("No_Inoc", "BMD"))

ceca_abx <- ggplot(c_t2, aes(x = Group.2, y = Group.3)) + 
  geom_tile(aes(x = Group.2, y = Group.3, fill = value)) +
  scale_fill_gradientn("CLR RA", 
                       breaks = seq(from = -2, to = 4, by = 1), 
                       colours = c("#FEBF2C",  "#C3003A",  "#551B44"), 
                       limits = c(-2,4)) +
  #facet_grid(. ~ Group.1) +
  facet_grid(Group.1 ~.) +
  geom_text(data = c_t2[c_t2$day13 < 0.05 & c_t2$Group.1 == "D13",], 
            aes(x = Group.2 , y = 1, label = "*"), col = "white", size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), axis.title.y = element_blank(),
        axis.title.x = element_blank())

x_total <- t(x_total)
x_total <- cbind(x_total, mapping[rownames(x_total), c("Treatment2", "Collection_Day")])
x_total2 <- melt(x_total)
unique(x_total2$variable)
x <- data.frame(t(working_table[keep_taxa,intersect(abx_no, Ceca), drop=F]))

x_total2$variable <- gsub("k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Buchnera; s__", "Buchnera 1", x_total2$variable)
x_total2$variable <- gsub("k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__Rickettsiaceae; g__Wolbachia; s__", "Wolbachia 1", x_total2$variable)
x_total2$variable <- gsub("k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__; s__", "Enterobacteriaseae 1", x_total2$variable)
x_total2$variable <- gsub("k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Arsenophonus; s__Arsenophonus endosymbiont of Dermacentor variabilis", "Arsenophonus 1", x_total2$variable)
x_total2$variable <- gsub(" k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Coprococcus; s__", "Coprococcus 1", x_total2$variable)
x_total2$Treatment2 <- factor(x_total2$Treatment2, levels = c("No_Inoc", "BMD"))
ceca_abx2 <- ggplot(x_total2, aes(x=Treatment2, y=value)) +
  geom_jitter(aes(color=Treatment2), width=0.15, alpha=0.8) +
  geom_boxplot(fill=NA, outlier.colour = NA) +
  facet_grid(variable ~ Collection_Day) +
  guides(color=F) +
  scale_color_manual(values = abx_col)

side1 <- plot_grid(NULL, ceca_a_alpha, ncol=2, rel_widths=c(0.1,1))
side2 <- plot_grid(side1, ceca_abx, nrow=2)
side <- plot_grid(ceca_a_alpha, NULL, ceca_abx, nrow=2, rel_widths = c(1,1.5))


topside <- plot_grid(NULL, ceca_a_alpha, ncol=2, rel_widths = c(0.2, 1))
side1 <- plot_grid(topside, ceca_abx, nrow=2, rel_heights = c(1,1))

bottom <- plot_grid(ceca_a_pcoa, side1, ncol=2, rel_widths = c(1.8,1.75))
#bottom <- plot_grid(ceca_a_pcoa, ceca_a_alpha, ceca_abx, ncol=3, rel_widths = c(1.8,0.65,1.1))


pdf(paste(main_fp, "/Supplemental_Figure2.pdf", sep=""), height=4, width=13, useDingbats =F)
print(bottom)
dev.off()

