#Make PCOA plots

mapping2 <- mapping
bray <- read.table(bray_fp, sep="\t", header=T, row=1)
uni <- read.table(uni_fp, sep="\t", header=T, row=1)
wuni <- read.table(wuni_fp, sep="\t", header=T, row=1)

beta_tables <- list(uni, wuni, bray)
beta_names <- c("unifrac.pdf", "weighted_unifrac.pdf", "bray_curtis.pdf")

######################################################################
#Plot day and bodysite separate###
#Save as one figure will all 9 combos
pcoa_func2 <- function(PCOA, header){
  color_by <- header
  print(header)
  plot1 <- ggplot(PCOA) +
    geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = color_by)) + 
    theme(legend.position = 'bottom') +
    scale_color_manual(values=cols2(length(unique(PCOA[,color_by])))) +
    guides(color=guide_legend(nrow=3))
}

for(b in 1:length(beta_tables)){
  plot_list <- c()
  beta_div <- as.data.frame(beta_tables[b])
  for(i in 1:length(Bodysites)){
    body <- Bodysites[[i]]
    for(k in 1:length(Days)){
      day <- Days[[k]]
      samples <- intersect(body,day)
      beta_table <-  beta_div[samples, samples]
      
      PCOA <- pcoa(beta_table)$vectors
      for(c in 1:ncol(PCOA)){
        colnames(PCOA)[c] <- paste("PC",c, sep="")
      }
      PCOA <- cbind(PCOA, rownames(PCOA))
      colnames(PCOA)[ncol(PCOA)] <- "SampleID"
      
      PCOA <- merge(PCOA, mapping2, by="SampleID")
      PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
      PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
      PCOA <- as.data.frame(PCOA)
      
      plot1 <- ggplot(PCOA) +
        geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Treatment2")) + 
        scale_color_manual(values=cols2(length(unique(PCOA[,"Treatment2"])))) 
      
      name <- paste(names(Days[k]), names(Bodysites[i]), sep="_")
      plot_list[[name]] <- plot1
    }
  }
  
  plot3by3 <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],
                        labels=c(names(plot_list)), vjust=-0.5, hjust=-0.1,ncol = 3)
  plot_this <- paste(main_fp, "beta_div/PCoa_SiteDay", beta_names[b], sep='/')
  save_plot(plot_this, plot3by3,
            ncol = 3,
            nrow = 3,
            base_aspect_ratio = 1.6)
  
}

######################################################################
#Plot Ileum Days 6 Unifrac, to see what taxon is driving the difference
#First make the PCOA table for this, then run a correlation test between PC1 and the rel. abundance of different OTUs
#If is is significant, then plot it using amount of that OTU as the color
beta_dir <- paste(main_fp, "beta_div/Unifrac/Ileum_6/", sep='/')

samples <- intersect(Ileum,Days[[2]])
beta_table <-  uni[samples, samples]
    
PCOA <- pcoa(beta_table)$vectors
for(c in 1:ncol(PCOA)){
  colnames(PCOA)[c] <- paste("PC",c, sep="")
}
PCOA <- cbind(PCOA, rownames(PCOA))
colnames(PCOA)[ncol(PCOA)] <- "SampleID"
PCOA <- as.data.frame(PCOA)    
#PCOA <- merge(PCOA, mapping2, by="SampleID")
PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]

rownames(PCOA) <- PCOA$SampleID

mapping_fortaxa <- mapping
mapping_fortaxa <- mapping_fortaxa[rownames(PCOA),]

working_taxa <- as.data.frame(t(taxa_table[,rownames(PCOA)]))
working_taxa <- working_taxa[,colSums(working_taxa)>0]
colnames(working_taxa) <- gsub("; ", "", colnames(working_taxa))
colnames(working_taxa) <- gsub("__", "_", colnames(working_taxa))
colnames(working_taxa) <- gsub("-", "_", colnames(working_taxa))
colnames(working_taxa) <- gsub("\\[", "", colnames(working_taxa))
colnames(working_taxa) <- gsub("\\]", "", colnames(working_taxa))
colnames(working_taxa) <- gsub(" ", "_",colnames(working_taxa))
working_taxa$SampleID <- rownames(working_taxa)

####
PCOA <- merge(PCOA, mapping_fortaxa, by="SampleID")
PCOA$PC1 <- as.numeric(as.character(PCOA$PC1))
PCOA$PC2 <- as.numeric(as.character(PCOA$PC2))
for(i in 4:ncol(PCOA)){
  PCOA[,i] <- as.factor(PCOA[,i])
}

for(i in 1:ncol(mapping_fortaxa)){
  header <- colnames(mapping_fortaxa)[i]
  color_by <- header
  print(header)
  name <- paste(header, "pdf", sep=".")
  fp <- paste(beta_dir, name, sep="/")
  plot1 <- ggplot(PCOA) +
    geom_point(size = 4, aes_string(x = "PC1", y = "PC2", color = color_by)) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10),
          legend.title = element_blank(),
          legend.key.size = unit(0.2, "in"),
          legend.text = element_text(size=5),
          legend.position = 'bottom',
          axis.text = element_text(size=5),
          axis.title = element_text(size=8)) +
    scale_color_manual(values=cols2(length(unique(PCOA[,color_by])))) +
    guides(color=guide_legend(nrow=3))
  
  pdf(fp, height=4,width=6)
  print(plot1)
  dev.off()
}

#keep taxa that are great than 0.01 in abundance across all samples
working_taxa <- working_taxa[,!names(working_taxa)=="SampleID"] #there are 365 OTUs
working_taxa <- working_taxa[,colSums(working_taxa) >0] #none dropped
working_taxa <- working_taxa[,colSums(working_taxa) > 0.001] # 160
working_taxa <- working_taxa[,colSums(working_taxa) > 0.01] # 63
working_taxa$SampleID <- rownames(working_taxa)

####Add Taxa quantiles to taxa map ###

ranges <- c(0, 0.0001, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
mapping_fortaxa <- merge(mapping_fortaxa, working_taxa, by="SampleID")
taxa_mod <- colnames(working_taxa)
for(i in 1:nrow(mapping_fortaxa)){
  sample_id <- rownames(mapping_fortaxa)[i]
  for(k in 1:length(taxa_mod)){
    taxon <- taxa_mod[k]
    for(m in 1:(length(ranges)-1)){
      if((mapping_fortaxa[sample_id, taxon] >= ranges[m]) && (mapping_fortaxa[sample_id, taxon] < ranges[(m+1)])){
        mapping_fortaxa[sample_id, taxon] <- ranges[m]
      } else {
        if(mapping_fortaxa[sample_id,taxon] == ranges[m]){
          mapping_fortaxa[sample_id, taxon] <- ranges[m]
        }
      }
    }
  }
}


##############################################
#If correlated with PC1 or 2, plot taxa gradient
PCOA2 <- merge(PCOA[,1:6], working_taxa, by="SampleID")

for(c in 7:ncol(PCOA2)){
  taxon <- colnames(PCOA2)[c]
  cor.out <- cor.test(PCOA2$PC1, PCOA2[,taxon])$p.val
  cor.out2 <- cor.test(PCOA2$PC2, PCOA2[,taxon])$p.val
  if(cor.out < 0.05 | cor.out2 < 0.05){
      #values <- sort(unique(PCOA2[,ncol(PCOA2)]))
      #PCOA2[,ncol(PCOA2)] <- as.factor(PCOA2[,ncol(PCOA2)])
      #PCOA2$new_column <- .bincode(PCOA2[,taxon], ranges, TRUE, TRUE)
    print(taxon)
    plot1 <- ggplot(PCOA2) +
      geom_point(size = 4, alpha=0.75, aes_string(x = "PC1", y = "PC2", color = taxon)) +
      theme_bw() +
      theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 10),
            legend.title = element_blank(),
            legend.key.size = unit(0.2, "in"),
            legend.text = element_text(size=5),
            legend.position = 'bottom',
            axis.text = element_text(size=5),
            axis.title = element_text(size=8)) +
            scale_color_gradientn(colors= c("#E6AB02", "#66A61E", "#666666"), guide="colorbar")
            #guides(color=guide_legend(nrow=3))  
  
      name <- paste(colnames(PCOA2)[c],".pdf", sep="_")
      save_plot(paste(beta_dir, name, sep=""), plot1)
  } 
}

######################################################################

pcoa_func <- function(PCOA, pcoa_dir){
  for(i in 1:ncol(mapping2)){
  header <- colnames(mapping2)[i]
  color_by <- header
  print(header)
  name <- paste(header, "pdf", sep=".")
  fp <- paste(pcoa_dir, name, sep="/")
  plot1 <- ggplot(PCOA) +
    geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = color_by)) + 
    # theme_bw() +
    # theme(panel.background = element_blank(),
    #       panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(),
    #       plot.title = element_text(size = 10),
    #       legend.title = element_blank(),
    #       legend.key.size = unit(0.2, "in"),
    #       legend.text = element_text(size=5),
    #       legend.position = 'bottom',
    #       axis.text = element_text(size=5),
    #       axis.title = element_text(size=8)) +
    theme(legend.position = 'bottom') +
    scale_color_manual(values=cols2(length(unique(PCOA[,color_by])))) +
    guides(color=guide_legend(nrow=3))
  
  #pdf(fp, height=4,width=6,useDingbats=FALSE)
  #print(plot1)
  #dev.off()
  save_plot(fp, plot1)
  }
}


####All Samples###

#Calculate coordinates
# pcoa_dirs <- c(paste(main_fp, "beta_div/BrayCurtis", sep='/'), paste(main_fp, "beta_div/Unifrac", sep='/'), paste(main_fp, "beta_div/WUnifrac", sep='/'))
# 
# PCOA <- pcoa(beta_div)$vectors
# for(i in 1:ncol(PCOA)){
#   colnames(PCOA)[i] <- paste("PC",i, sep="")
# }
# PCOA <- cbind(PCOA, rownames(PCOA))
# colnames(PCOA)[ncol(PCOA)] <- "SampleID"
# 
# PCOA <- merge(PCOA, mapping2, by="SampleID")
# PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
# PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
# PCOA <- as.data.frame(PCOA)
# 
# pcoa_func(PCOA, pcoa_dir)
# 
# 
# ##Plot each bodysite separate###
# 
# for(i in 1:length(Bodysites)){
#   body <- Bodysites[[i]]
#   beta_table <-  beta_div[body, body]
#   
#   PCOA <- pcoa(beta_table)$vectors
#   for(c in 1:ncol(PCOA)){
#     colnames(PCOA)[c] <- paste("PC",c, sep="")
#   }
#   PCOA <- cbind(PCOA, rownames(PCOA))
#   colnames(PCOA)[ncol(PCOA)] <- "SampleID"
#   
#   PCOA <- merge(PCOA, mapping2, by="SampleID")
#   PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
#   PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
#   PCOA <- as.data.frame(PCOA)
#   
#   dir <- paste("beta_div/BrayCurtis/pcoa/", names(Bodysites[i]), sep="")
#   pcoa_dir <- paste(main_fp, dir, sep='/')
#   
#   pcoa_func(PCOA, pcoa_dir)
# }
#
