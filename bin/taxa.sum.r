#Make taxa summaries

source("bin/taxa.sum.func.r")

taxa_table1 <- taxa_table
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
    taxon <- (strsplit(taxon, "; ",))[[1]][5]
  }
  if(taxon == "f__"){
    taxon <- (strsplit(otu_names[i], holder, fixed=T))[[1]][1]
    taxon <- (strsplit(taxon, "; ",))[[1]][4]
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

cutoff <- 0.05
# Creat one plot for just no-inoc controls, one plot per day
for(k in 1:length(Bodysites)){
  body_test <- Bodysites[[k]]
  body <- names(Bodysites[k])
  union1 <- intersect(body_test, NoInoc)
  days_test <- c("D03", "D06", "D13")
  plot_list <- list()
  #Make plot
  for(d in 1:length(days_test)){
    plot_title <- body
    this <- days_test[d]
    union <- intersect(union1, Days[[d]])
    otu1 <- make_taxa_sums(otu, union, cutoff)
    otu1$Treatment2 <- factor(otu1$Treatment2, levels= c("No_Inoc"))
    otu1$Collection_Day <- as.character(otu1$Collection_Day)
    taxa_plot <- ggplot(otu1, aes(x = SampleID , y = Relative_Abundance)) + 
      geom_bar(stat="identity",aes(fill=Taxa)) +
      facet_wrap(facets=~Collection_Day, scales = "free_x", nrow=1) +
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
      guides(fill=F) +
      scale_fill_manual(name= names(taxa_cols), values= taxa_cols)
  plot_list[[this]] <- taxa_plot
  }
  plot_all <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol=3)
  #assign pdf name for plot
  name <- paste(body, ".pdf", sep='')
  file_path <- paste(taxa_dir, name, sep='')
  pdf(file_path, height=4,width=11)
  print(plot_all)
  dev.off()
}



for(k in 1:length(Bodysites)){
  body_test <- Bodysites[[k]]
  body <- names(Bodysites[k])
  for(m in 1:length(Days)){
    day_test <- Days[[m]]
    day <- names(Days[m])
    union <- intersect(body_test, day_test)
    union <- intersect(union, colnames(otu))
    otu1 <- make_taxa_sums(otu, union, cutoff)
    otu1$Treatment2 <- factor(otu1$Treatment2, levels= c("No_Inoc", "GroGel", "TJPbx", "FMB11", "BMD"))
    
    #Make plot
    plot_title <- sprintf('%s, %s',body, day)
    taxa_plot <- ggplot(otu1, aes(x = SampleID , y = Relative_Abundance)) + 
      geom_bar(stat="identity",aes(fill=Taxa)) +
      facet_wrap(facets=~Treatment2, scales = "free_x", nrow=1) +
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
      guides(fill=guide_legend(ncol=5)) +
      scale_fill_manual(name= names(taxa_cols), values= taxa_cols)
    
    #assign pdf name for plot
    name <- paste(body, day, sep="-")
    name <- paste(name, ".pdf", sep='')
    file_path <- paste(taxa_dir, name, sep='')
    pdf(file_path, height=4,width=11)
    print(taxa_plot)
    dev.off()
  }
}

union <- samples_no_con
union <- intersect(union, colnames(otu))
cutoff <- 0.05
otu1 <- make_taxa_sums(otu, union, cutoff)

plot_title <- "Effect of Day, regardless of treatment"
taxa_plot <- ggplot(otu1, aes(x = Collection_Day , y = Relative_Abundance)) + 
  geom_bar(stat="identity",position="fill", aes(fill=Taxa)) +
  facet_wrap(facets=~Tissue, scales = "free_x", nrow=1) +
  theme_bw() +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "in"),
        legend.text = element_text(size=5),
        legend.position = 'bottom',
        axis.text = element_text(size=5),
        axis.title = element_text(size=8)) + 
  guides(fill=guide_legend(nrow=6)) +
  scale_fill_manual(name= names(taxa_cols), values= taxa_cols)

#assign pdf name for plot
name <- "Day_Facet_Tissue.pdf"
file_path <- paste(taxa_dir, name, sep='')
pdf(file_path, height=4,width=11)
print(taxa_plot)
dev.off()

#Compare treatment to input
cutoff <- 0.05
for(i in 1:length(Pbx)){
  treatment_plot <- Treatments[[i]]
  treatment_plot <- intersect(treatment_plot, colnames(otu))
  input_plot <- Inputs[[i]]
  input_plot <- intersect(input_plot, colnames(otu))
  otu1 <- make_taxa_sums(otu, treatment_plot, cutoff)
  otu2 <- make_taxa_sums(otu, input_plot, cutoff)
  
  plot_title <- names(Pbx[i])
  taxa_plot1 <- ggplot(otu1, aes(x = SampleID , y = Relative_Abundance)) + 
    geom_bar(stat="identity",position="fill", aes(fill=Taxa)) +
    facet_wrap(facets=~Collection_Day, scales = "free_x", nrow=1) +
    scale_fill_manual(values= taxa_cols)
  
  taxa_plot2 <- ggplot(otu2, aes(x = SampleID , y = Relative_Abundance)) + 
    geom_bar(stat="identity",position="fill", aes(fill=Taxa)) +
    scale_fill_manual(values= taxa_cols)
  
  #assign pdf name for plot
  name <- paste(plot_title, ".pdf", sep="")
  file_path <- paste(taxa_dir, name, sep='')
  
  plot_together <- plot_grid(taxa_plot2, taxa_plot1, rel_widths = c(2,5), 
                             labels=c("input", "samples"), vjust=-0.5, hjust=-0.1,ncol = 2)
  save_plot(file_path, plot_together,
            ncol = 2,
            base_aspect_ratio = 1.8)
}
