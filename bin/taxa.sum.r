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

for(k in 1:length(Bodysites)){
  body_test <- Bodysites[[k]]
  body <- names(Bodysites[k])
  for(m in 1:length(Days)){
    day_test <- Days[[m]]
    day <- names(Days[m])
    union <- intersect(body_test, day_test)
    otu1 <- make_taxa_sums(otu, union, cutoff)
    otu1$Treatment <- factor(otu1$Treatment2, levels= c("No_Inoc", "BMD", "GroGel", "TJPbx", "FMB11"))
    
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
cutoff <- 0.20
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
name <- "Day_Facet_Tissue"
file_path <- paste(taxa_dir, name, sep='')
pdf(file_path, height=4,width=11)
print(taxa_plot)
dev.off()

#Compare treatment to input
cutoff <- 0.1
for(i in 1:length(Treatments)){
  treatment_plot <- Treatments[[i]]
  input_plot <- Inputs[[i]]
  otu1 <- make_taxa_sums(otu, treatment_plot, cutoff)
  otu2 <- make_taxa_sums(otu, input_plot, cutoff)
  
  plot_title <- names(Treatments[i])
  taxa_plot1 <- ggplot(otu1, aes(x = SampleID , y = Relative_Abundance)) + 
    geom_bar(stat="identity",position="fill", aes(fill=Taxa)) +
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
    guides(fill=guide_legend(nrow=6)) +
    scale_fill_manual(name= names(taxa_cols), values= taxa_cols)
  
  taxa_plot2 <- ggplot(otu2, aes(x = SampleID , y = Relative_Abundance)) + 
    geom_bar(stat="identity",position="fill", aes(fill=Taxa)) +
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
  name <- paste(plot_title, ".pdf", sep="")
  file_path <- paste(taxa_dir, name, sep='')
  pdf(file_path, height=4,width=11)
  grid.arrange(taxa_plot2, taxa_plot1, ncol=2, nrow=1)
  dev.off()
}
