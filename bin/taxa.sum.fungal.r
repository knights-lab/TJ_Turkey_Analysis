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



############Make taxa plot for each bodysite on each day, facet by treatment
taxa_dir <- paste(main_fp, "taxa_sum_fungal/", sep='/')

cutoff <- 0.05

for(m in 1:length(Days_f)){
  day_test <- Days_f[[m]]
  day <- names(Days_f[m])
  otu1 <- make_taxa_sums2(otu, day_test, cutoff)
  otu1$Treatment <- factor(otu1$Treatment, levels= c("1_NoInoc", "2_GroGelNoPbx", "4_TJPbxGroGel", "3_FMB11GroGel", "5_BMD"))
    
  #Make plot
  plot_title <- sprintf('%s', day)
  taxa_plot <- ggplot(otu1, aes(x = SampleID , y = Relative_Abundance)) + 
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
    guides(fill=guide_legend(ncol=5)) +
    scale_fill_manual(name= names(ftaxa_cols), values= ftaxa_cols)
    
  #assign pdf name for plot
  name <- paste(day, ".pdf", sep='')
  file_path <- paste(taxa_dir, name, sep='')
  pdf(file_path, height=13,width=20)
  print(taxa_plot)
  dev.off()
}


##### Look at all samples, color by day
union <- colnames(otu)
cutoff <- 0.05
otu1 <- make_taxa_sums2(otu, union, cutoff)

plot_title <- "Effect of Day, regardless of treatment"
taxa_plot <- ggplot(otu1, aes(x = SampleID , y = Relative_Abundance)) + 
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
        axis.title = element_text(size=8)) + 
  guides(fill=guide_legend(nrow=6)) +
  scale_fill_manual(name= names(ftaxa_cols), values= ftaxa_cols)

#assign pdf name for plot
name <- "Fungal_samples_Day.pdf"
file_path <- paste(taxa_dir, name, sep='')
pdf(file_path, height=13,width=20)
print(taxa_plot)
dev.off()
