#Treatment also modifies Fungi
source("bin/load.fungal.r")
#1.Make PCOA plots
f_map2 <- f_map

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

beta_table <- as.matrix(dist(t(test)))

#Run stats for diff. centroids
beta_dist = as.dist(beta_table)
ad = adonis(beta_dist ~ f_map2[colnames(test),"Treatment2"], data=f_map2, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
#Run Stats for diff. dispersion
beta_out <- betadisper(beta_dist, f_map2[colnames(test),"Treatment2"])
p_val_disp <- permutest(beta_out)$tab[1, 6]
  
PCOA <- pcoa(beta_table)$vectors
var_exp <- pcoa(beta_table)$values
for(c in 1:ncol(PCOA)){
  colnames(PCOA)[c] <- paste("PC",c, sep="")
}
PCOA <- cbind(PCOA, rownames(PCOA))
colnames(PCOA)[ncol(PCOA)] <- "SampleID"
  
PCOA <- merge(PCOA, f_map2, by="SampleID")
PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
PCOA <- as.data.frame(PCOA)
PCOA$Treatment2 <- factor(PCOA$Treatment2, levels = c("BMD", "TJPbx", "FMB11", "No_Inoc", "GroGel"))

centroids <- aggregate(cbind(PCOA$PC1,PCOA$PC2) ~ PCOA$Treatment2,PCOA,mean)
colnames(centroids) <- c("Treatment2", "PC1", "PC2")

#pairwise adonis for differences:
pair_ad_list <- list()
for(t in 1:4){
  for(u in (t+1):5){
    t1 <- unique(f_map2$Treatment2)[[t]]
    t2 <- unique(f_map2$Treatment2)[[u]]
    sams <- rownames(f_map2[f_map2$Treatment2 == t1 | f_map2$Treatment2 == t2,])
    sams <- intersect(rownames(beta_table), sams)
    sub_table <- beta_table[sams, sams]
    beta_dist = as.dist(sub_table)
    ad = adonis(beta_dist ~ f_map2[sams,"Treatment2"], data=f_map2, permutations=999)
    p_val <- ad$aov.tab[1,6]
    r_sq <- ad$aov.tab[1,5]
    pair_ad_list[[paste(t1, "_", t2, "Fungal", sep="")]] <- c("pval", p_val, "r_sq", r_sq)
  }
}

sink(paste(main_fp, "/fungal_treatment_pairwise_centroids.txt", sep=""))
pair_ad_list
sink()

fungal_6 <- ggplot(PCOA) +
    geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Treatment2"), alpha =0.8) + 
    scale_color_manual(values= total_col) +
    scale_fill_manual(values= total_col) +
    stat_ellipse(aes_string(x = "PC1", y = "PC2", color = "Treatment2"), linetype = 2) +
    annotate("text", x=0.2, y=0.25, label= paste("P=", p_val), size=2) +
    annotate("text", x=0.2, y=0.20, label= paste("R2=", round(r_sq, digits=3)), size=2) +
    coord_fixed(xlim = c(-10, 10), ylim = c(-10, 10)) +
    labs(x=paste("PC1 (", round(var_exp$Relative_eig[1], digits=3)*100, "%)", sep=""), y= paste("PC2 (", round(var_exp$Relative_eig[2], digits=3)*100, "%)", sep="")) +
    geom_point(data=centroids, aes(x=PC1, y=PC2, fill=Treatment2), shape = 23, size=3, alpha=0.8, show.legend = F)

######################################################################
#taxa summary showing the major members decreasing

#Make a taxa table
ft_table <- fungal_otu3
ftaxonomy2 <- ftaxonomy[rownames(ft_table),]
rownames(ft_table) <- ftaxonomy2$taxonomy7 # 205 x 221

#Collapse same taxonomies (from 142 to 123)
ftaxa_table <- aggregate(ft_table[,1:ncol(ft_table)], by=list(rownames(ft_table)), FUN = sum)
rownames(ftaxa_table) <- ftaxa_table[,1]
ftaxa_table <- ftaxa_table[,2:ncol(ftaxa_table)] #172 x 221

#Assign taxa colors
ftaxa_cols <- sample(cols2(length(rownames(ftaxa_table))))
names(ftaxa_cols) <- rownames(ftaxa_table)
ftaxa_cols <- c(ftaxa_cols, "#d3d3d3")
names(ftaxa_cols)[nrow(ftaxa_table)+1] <- "Other"

#convert to RA
otu <- sweep(ftaxa_table,2,colSums(ftaxa_table),`/`)

f_cols <- c( "#fad4e7", "#D95F02", "#361b9e", "#771b9e", "#666666", "#E6AB02", "#c978a3", "#A6761D", "#731445", "#7570B3",  
             "#E7298A", "#1B9E77", "#1e66a6","#92a49e", 
               "#66A61E","#a5a2ce" , "#a3c978","#e3e2ef",
               "#b37570", "#123123", "#FFFF99", "#d3d3d3")

cutoff <- 0.2

day_test <- Days_f[[2]]
day <- names(Days_f[2])

source("bin/taxa.sum.func.r")
otu1 <- make_taxa_sums2(otu, day_test, cutoff)

#otu1 <- make_taxa_sums2(ft_table, day_test, cutoff)
otu1$Treatment <- factor(otu1$Treatment, levels= c("1_NoInoc", "5_BMD", "4_TJPbxGroGel", "3_FMB11GroGel", "2_GroGelNoPbx"))
  
#Make plot
taxa_plot_f <- ggplot(otu1, aes(x = SampleID , y = Relative_Abundance)) + 
    geom_bar(stat="identity",aes(fill=Taxa)) +
    facet_wrap(facets=~Treatment, scales = "free_x", nrow=1) +
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
    guides(fill=guide_legend(ncol=9)) +
    scale_fill_manual(name= names(ftaxa_cols), values= f_cols)

test_taxa <- as.character(unique(otu1$Taxa))


######################################################################
#Plot the CLR relative abundance of taxa that are different:
####Filter####
#drop samples below 1500 counts
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

#add taxonomy 
ft_table <- day6_f
ftaxonomy2 <- ftaxonomy[colnames(ft_table),]
colnames(ft_table) <- ftaxonomy2$taxonomy7 
ft_table <- t(ft_table)
#Collapse same taxonomies 
ftaxa_table <- aggregate(ft_table[,1:ncol(ft_table)], by=list(rownames(ft_table)), FUN = mean)
rownames(ftaxa_table) <- ftaxa_table[,1]
ftaxa_table <- ftaxa_table[,2:ncol(ftaxa_table)] 


plot_these <- data.frame(t(ftaxa_table[rownames(ftaxa_table) %in% test_taxa,]))

#Set up tests to run
NoInoc_f <- intersect(row.names(f_map[f_map$Treatment == "1_NoInoc",]), rownames(plot_these))
FMB11_f <- intersect(row.names(f_map[f_map$Treatment == "3_FMB11GroGel",]), rownames(plot_these))
BMD_f <- intersect(row.names(f_map[f_map$Treatment == "5_BMD",]), rownames(plot_these))
TJPbx_f <- intersect(row.names(f_map[f_map$Treatment == "4_TJPbxGroGel",]), rownames(plot_these))
GroGel_f <- intersect(row.names(f_map[f_map$Treatment == "2_GroGelNoPbx",]), rownames(plot_these))

test.ixs <- list(NoInoc_f, FMB11_f, BMD_f, TJPbx_f, GroGel_f)
names(test.ixs) <- c("NoInoc", "FMB11", "BMD", "TJPbX", "GroGel")

pvals<- c()
ALPHA <- 0.05
diff_list <- list()

  for(n in 1:(length(test.ixs)-1)){
    for(m in (n+1):length(test.ixs)){
    test.x <- test.ixs[n]
    test.y <- test.ixs[m]
    working_table <- t(plot_these)
    set1 <- intersect(test.x[[1]], colnames(working_table))
    set2 <- intersect(test.y[[1]], colnames(working_table))
    if(length(set1) > 2 && length(set2) > 2){
      full_set <- c(set1, set2)
      test_table <- working_table[, full_set]
    
      test_table <- t(test_table)
      
      #Keep taxa that have a mean rel abundance of at least 0.01
      #test_table <- test_table[,colMeans(test_table)>0.009, drop=F]
      map_test <- f_map[rownames(test_table),]
      difftest <- test.otu.features(test_table, response=map_test$Treatment, sig.level = 0.10)
      #difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
      
      if(any(difftest$pvals <= ALPHA)){
        signif.ix <- which(difftest$pvals <= ALPHA)
        signif.ix <- signif.ix[order(difftest$pvals[signif.ix])]
        sink("list of diff fungal taxa.txt", append=T)
        cat(paste(names(test.ixs)[n], " vs ", names(test.ixs)[m], "\n"))
        print(table(names(signif.ix)))
        cat("\n")
        sink()
        diff_list[[paste(names(test.ixs)[n], names(test.ixs)[m])]] <- names(signif.ix)
      } else {
        cat("Not Signif.")
      }
    }
  }
}

keep_taxa <- unique(c(diff_list[["NoInoc BMD"]], diff_list[["TJPbX GroGel"]]))
x <- data.frame(t(plot_these[, keep_taxa]))
x_totalf <- x
x$BMD <- 1
x$TJPbx <- 1
x[rownames(x) %in% diff_list[["NoInoc BMD"]], "BMD"] <- 0.04
x[rownames(x) %in% diff_list[["TJPbX GroGel"]], "TJPbx"] <- 0.04

i_table <- data.frame(t(x))
ps <- i_table[c("BMD", "TJPbx"),]
i_table <- i_table[!rownames(i_table) %in% c("BMD", "TJPbx"),]
i_map <- f_map[rownames(i_table),]
i_table$Treatment <- i_map$Treatment
i_table <- melt(i_table)
i_t2 <-aggregate(i_table, by=list(i_table$Treatment, i_table$variable), 
                 FUN=mean, na.rm=TRUE)
ps <- data.frame(t(ps))
ps$Group.2 <- rownames(ps)
i_t2 <- merge(i_t2, ps, by="Group.2")
taxa_order <- i_t2[i_t2$Group.1 == "5_BMD",]
taxa_order <- taxa_order[order(-taxa_order$value),]
i_t2$Group.2 <- factor(i_t2$Group.2, levels=taxa_order$Group.2)
i_t2$Group.1 <- factor(i_t2$Group.1, levels = c("5_BMD", "1_NoInoc", "4_TJPbxGroGel", "3_FMB11GroGel", "2_GroGelNoPbx"))
i_t2$A_P <- as.character(i_t2$Group.1)
i_t2$A_P <- gsub("5_BMD", "Antibiotic", i_t2$A_P)
i_t2$A_P <- gsub("1_NoInoc", "Antibiotic", i_t2$A_P)
i_t2$A_P <- gsub("4_TJPbxGroGel", "Probiotic", i_t2$A_P)
i_t2$A_P <- gsub("3_FMB11GroGel", "Probiotic", i_t2$A_P)
i_t2$A_P <- gsub("2_GroGelNoPbx", "Probiotic", i_t2$A_P)
i_t2$A_P <- factor(i_t2$A_P)
i_t2$Group.2 <- gsub("\\.", " ", i_t2$Group.2)
ileum_f_taxa <- ggplot(i_t2) + 
  geom_tile(aes(x = Group.1, y = Group.2, fill = value)) +
  scale_fill_gradientn("CLR RA", 
                       breaks = seq(from = min(i_t2$value), to = max(i_t2$value), by = 0.5), 
                       colours = c("#FEBF2C",  "#C3003A",  "#551B44"), 
                       limits = c(min(i_t2$value), max(i_t2$value))) +
  #facet_grid(.~ A_P, scales="free", space="free") +
  geom_text(data = i_t2[i_t2$BMD < 0.05,], 
            aes(x = 1, y = Group.2, label = "*"), col = "white", size = 4) +
  geom_text(data = i_t2[i_t2$TJPbx < 0.05,], 
            aes(x = 3, y = Group.2, label = "*"), col = "white", size = 4) +
  labs(x="Treatment", y="Taxa")

x_totalf <- t(x_totalf)
x_totalf <- cbind(x_totalf, f_map[rownames(x_totalf), c("Treatment2", "Treatment")])
x_totalf2 <- melt(x_totalf)
unique(x_totalf2$variable)

x_totalf2$Treatment2 <- factor(x_totalf2$Treatment2, levels = c("BMD", "No_Inoc", "TJPbx", "GroGel", "FMB11"))

ileum_f_taxa2 <- ggplot(x_totalf2, aes(x=Treatment2, y=value)) +
  geom_jitter(aes(color=Treatment2), width=0.15, alpha=0.8) +
  geom_boxplot(fill=NA, outlier.colour = NA) +
  facet_grid(variable~.) +
  guides(color=F) +
  scale_color_manual(values = cols)

#No Diff Taxa when you run on any other day
######################################################################
#compile figure
bottom_f <- plot_grid(fungal_6, ileum_f_taxa, ncol=2, rel_widths=c(0.8, 1))
f_together <- plot_grid(taxa_plot_f, bottom_f, nrow=2, rel_heights = c(1, 0.7))
pdf(paste(main_fp, "/Figure_4.pdf",sep=""), height=7, width=9)
f_together
dev.off()


######################################################################
#Supplemental Figure - Alpha Diversity
alpha_div2 <- alpha_div2[rownames(f_map),]
f_map$shannon <- alpha_div2[,"shannon"]
f_map$simpson <- alpha_div2[,"simpson"]
#f_map$pdwholetree <- alpha_div2[,"PD_whole_tree"]
f_map$obsspecies <- alpha_div2[,"observed_species"]

Controls_alpha <- f_map[f_map$Treatment == "1_NoInoc",]

#run stats
pval_list <- list()
for(d in 1:(length(Days_f)-1)){
  for(e in (d+1):length(Days_f)){
      day1 <- Days_f[[d]]
      day2 <- Days_f[[e]]
      samples1 <- intersect(day1, rownames(Controls_alpha))
      samples2 <- intersect(day2, rownames(Controls_alpha))
      a_test <- t.test(Controls_alpha[samples1, "shannon"], Controls_alpha[samples2, "shannon"])$p.value
      pval_list[[paste(names(Days_f)[d], names(Days_f)[e], sep="_")]] <- a_test
  }
}
sink(paste(main_fp, "s_Fig4_alpha_stats.txt", sep="/"))
pval_list
sink()

#run stats as linear correlation 
Controls_alpha$Day <- Controls_alpha$Collection_Day
Controls_alpha$Day <- as.numeric(gsub("D", "", Controls_alpha$Day))
l_test <- summary(lm(Controls_alpha[, "shannon"] ~ Controls_alpha[, "Day"]))

sink(paste(main_fp, "s_Fig4_alpha_stats.txt", sep="/"), append=T)
l_test
sink()

alpha_div_plot_c <- ggplot(Controls_alpha, aes(x=Day, y=shannon))+
  geom_jitter(alpha=0.8, width=.15, aes(group=Day, color=Collection_Day)) +
  geom_boxplot(fill=NA, aes(group=Day)) +
  scale_color_manual(values = cols3) +
  guides(color=F) +
  geom_smooth(color=cols3[2], se=F) +
  scale_x_continuous(breaks=c(3,6,13))

f_map$Treatment2 <- factor(f_map$Treatment2, levels=c("BMD", "No_Inoc", "TJPbx", "GroGel", "FMB11"))
alpha_div_plot_treat <- ggplot(f_map, aes(x=Treatment2, y=shannon, color=Treatment2)) +
  geom_jitter(width = 0.15) +
  geom_boxplot(fill=NA, color="black", outlier.colour = NA) +
  scale_color_manual(values=c("#D95F02", "#1B9E77", "#7570B3", "#E6AB02","#E7298A")) +
  facet_grid(.~Collection_Day) +
  guides(color=F)

#test stats for treatment vs control - ignore other ones that are printed
pval_list <- list()
controlsf <- c("No_Inoc", "GroGel")
testsf <- c("BMD", "TJPbx", "FMB11")
for(t in 1:length(Days_f)){
  days1 <- Days_f[[t]]
  for(i in controlsf){
    for(k in testsf){
      sams1 <- intersect(days1, rownames(f_map[f_map$Treatment2 == i,]))
      sams2 <- intersect(days1, rownames(f_map[f_map$Treatment2 == k,]))
      a_test <- t.test(f_map[sams1, "shannon"], f_map[sams2, "shannon"])$p.value
      pval_list[[paste(names(Days_f)[t], i,k, sep="_")]] <- a_test
    }
  }
}
sink(paste(main_fp, "s_Fig4_alpha_stats.txt", sep="/"), append=T)
pval_list
sink()

#PCOA for all 3 days
fCLR_otutable <- fotu_table3
fCLR_otutable[fCLR_otutable == 0] <- 0.65 #Convert any 0 to 0.65 to allow for CLR transform

fCLR_otutable <- t(fCLR_otutable) #convert to samples as rows
fCLR_otutable <- cenLR(fCLR_otutable)$x.clr  #transform

plot_list <- list()
for(i in 1:length(Days_f)){
  day1 <- Days_f[[i]]
  day6_f <- fCLR_otutable[day1, ]
  keep_f <- fungal_otu3[rowSums(fungal_otu3 > 0) > 5, rownames(day6_f)]
  keep_f <- intersect(rownames(keep_f), colnames(day6_f))
  day6_f <- day6_f[rowSums(day6_f > 0) > 1, keep_f]
  test <- t(day6_f)
  
  beta_table <- as.matrix(dist(t(test)))
  
  #Run stats for diff. centroids
  beta_dist = as.dist(beta_table)
  ad = adonis(beta_dist ~ f_map2[colnames(test),"Treatment2"], data=f_map2, permutations=999)
  p_val <- ad$aov.tab[1,6]
  r_sq <- ad$aov.tab[1,5]
  #Run Stats for diff. dispersion
  beta_out <- betadisper(beta_dist, f_map2[colnames(test),"Treatment2"])
  p_val_disp <- permutest(beta_out)$tab[1, 6]
  
  PCOA <- pcoa(beta_table)$vectors
  var_exp <- pcoa(beta_table)$values
  for(c in 1:ncol(PCOA)){
    colnames(PCOA)[c] <- paste("PC",c, sep="")
  }
  PCOA <- cbind(PCOA, rownames(PCOA))
  colnames(PCOA)[ncol(PCOA)] <- "SampleID"
  
  PCOA <- merge(PCOA, f_map2, by="SampleID")
  PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
  PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
  PCOA <- as.data.frame(PCOA)
  PCOA$Treatment2 <- factor(PCOA$Treatment2, levels = c("BMD", "TJPbx", "FMB11", "No_Inoc", "GroGel"))
  
  fungal_pcoa1 <- ggplot(PCOA) +
    geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Treatment2"), alpha =0.75) + 
    scale_color_manual(values= total_col) +
    scale_fill_manual(values= total_col) +
    stat_ellipse(aes_string(x = "PC1", y = "PC2", color = "Treatment2"), linetype = 2) +
    annotate("text", x=0.2, y=0.25, label= paste("P=", p_val), size=2) +
    annotate("text", x=0.2, y=0.20, label= paste("R2=", round(r_sq, digits=3)), size=2) +
    labs(x=paste("PC1 (", round(var_exp$Relative_eig[1], digits=3)*100, "%)", sep=""), y= paste("PC2 (", round(var_exp$Relative_eig[2], digits=3)*100, "%)", sep=""))
  
  name <- names(Days_f[i])
  plot_list[[name]] <- fungal_pcoa1
}


######################################################################
#add taxonomy 
ft_table <- fCLR_otutable
ftaxonomy2 <- ftaxonomy[colnames(ft_table),]
colnames(ft_table) <- ftaxonomy2$taxonomy7 
ft_table <- t(ft_table)
#Collapse same taxonomies 
ftaxa_table <- aggregate(ft_table[,1:ncol(ft_table)], by=list(rownames(ft_table)), FUN = mean)
rownames(ftaxa_table) <- ftaxa_table[,1]
ftaxa_table <- ftaxa_table[,2:ncol(ftaxa_table)] 
new_ft<- t(ftaxa_table)

test_this <- melt(new_ft)

f_map$Var1 <- f_map$SampleID
test_this2 <- merge(test_this, f_map, by="Var1")
test_this2$Treatment2 <- factor(test_this2$Treatment2, levels= c("No_Inoc", "BMD", "TJPbx", "FMB11", "GroGel"))
taxa_heat <- heat(test_this2, "Var1", "Var2", fill = "value") +
  facet_grid(.~Collection_Day + Treatment2, scales="free", space = "free") +
  theme(axis.text.x=element_blank())



plot_taxaf <- plot_grid(plot_list_t[[1]], plot_list_t[[2]],plot_list_t[[3]], nrow=3)

##### Compile Figure
plot3f <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]],
                    labels=c(names(plot_list)), ncol = 2)
plot_bottom <- plot_grid(plot3f, NULL, ncol=2, rel_widths = c(3,1))

fig_s4 <- plot_grid(alpha_div_plot_c, alpha_div_plot_treat, ncol=2, rel_widths = c(1,3))
fig_s4a <- plot_grid(fig_s4, plot_bottom, ncol=1, rel_heights = c(1,3))

pdf(paste(main_fp, "/Supplemental_Figure4.pdf", sep=""), width=8, height=5, useDingbats = F)
fig_s4a
dev.off()


pdf(paste(main_fp, "/Supplemental_Figure5.pdf", sep=""), width=12, height=6, useDingbats = F)
taxa_heat
dev.off()
