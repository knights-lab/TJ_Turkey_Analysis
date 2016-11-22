#Make PCOA plots
mapping2 <- mapping

pcoa_func <- function(PCOA, pcoa_dir){
  for(i in 1:ncol(mapping2)){
  header <- colnames(mapping2)[i]
  color_by <- header
  print(header)
  name <- paste(header, "pdf", sep=".")
  fp <- paste(pcoa_dir, name, sep="/")
  plot1 <- ggplot(PCOA) +
    geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = color_by)) + 
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
  
  pdf(fp, height=4,width=6,useDingbats=FALSE)
  print(plot1)
  dev.off()
  }
}

####All Samples###

#Calculate coordinates
pcoa_dir <- paste(main_fp, "beta_div/BrayCurtis/pcoa/", sep='/')

PCOA <- pcoa(beta_div)$vectors
for(i in 1:ncol(PCOA)){
  colnames(PCOA)[i] <- paste("PC",i, sep="")
}
PCOA <- cbind(PCOA, rownames(PCOA))
colnames(PCOA)[ncol(PCOA)] <- "SampleID"

PCOA <- merge(PCOA, mapping2, by="SampleID")
PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
PCOA <- as.data.frame(PCOA)

pcoa_func(PCOA, pcoa_dir)



##Plot each bodysite separate###

for(i in 1:length(Bodysites)){
  body <- Bodysites[[i]]
  beta_table <-  beta_div[body, body]
  
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
  
  dir <- paste("beta_div/BrayCurtis/pcoa/", names(Bodysites[i]), sep="")
  pcoa_dir <- paste(main_fp, dir, sep='/')
  
  pcoa_func(PCOA, pcoa_dir)
}

