#Ileum and cecum day 6
otu_test <- read.table("data/Turkey_10k_Taxa.txt", 
                       sep="\t", 
                       comment="", 
                       header=T, 
                       skip=1, 
                       as.is=T, 
                       check.names=F, row=1)
tree = read_tree_greengenes('data/97_otus.tree')

m <- mapping[mapping$Tissue == "Ileum" | mapping$Tissue == "Ceca",]
beta_keeps <- intersect(rownames(m[Days[[2]],]), colnames(otu_test))

total_col <- c("#D95F02", "#7570B3", "#E7298A", "#1B9E77", "#E6AB02")
# "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"
plot_list <- list()
I_C <- list(Bodysites[[2]], Bodysites[[3]])
names(I_C) <- c("Ileum", "Ceca")
for(i in 1:length(I_C)){
  site <- I_C[[i]]
  samples <- intersect(site, Days[[2]])
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
  PCOA <- pcoa(uni1)$vectors
  var_exp <- pcoa(uni1)$values
  
  #Run stats for diff. centroids
  beta_dist = uni1
  ad = adonis(beta_dist ~ m[,"Treatment2"], data=m, permutations=999)
  p_val <- ad$aov.tab[1,6]
  r_sq <- ad$aov.tab[1,5]
  #Run Stats for diff. dispersion
  beta_out <- betadisper(beta_dist, m$Treatment2)
  p_val_disp <- permutest(beta_out)$tab[1, 6]
  
  sil_list <- c()
  for(k in 2:5){
    this <- pam(beta_dist,k)
    sil_out <- summary(silhouette(this$clustering, beta_dist))
    av_widths <- sil_out$avg.width
    sil_list <- c(sil_list, av_widths)
  }
  
  
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
  PCOA$Treatment2 <- factor(PCOA$Treatment2, levels=c("BMD", "TJPbx", "FMB11", "No_Inoc", "GroGel"))
  
  centroids <- aggregate(cbind(PCOA$PC1,PCOA$PC2) ~ PCOA$Treatment2,PCOA,mean)
  colnames(centroids) <- c("Treatment2", "PC1", "PC2")
  
  #pairwise adonis for differences:
  pair_ad_list <- list()
  for(t in 1:4){
    for(u in (t+1):5){
      t1 <- unique(m$Treatment2)[[t]]
      t2 <- unique(m$Treatment2)[[u]]
      sams <- rownames(m[m$Treatment2 == t1 | m$Treatment2 == t2,])
      beta_dist = as.dist(as.matrix(uni1)[sams,sams])
      ad = adonis(beta_dist ~ m[sams,"Treatment2"], data=m, permutations=999)
      p_val <- ad$aov.tab[1,6]
      r_sq <- ad$aov.tab[1,5]
      pair_ad_list[[paste(t1, "_", t2, names(I_C)[i], sep="")]] <- c("pval", p_val, "r_sq", r_sq)
    }
  }
  
  plot1 <- ggplot(PCOA) +
    geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Treatment2"), alpha =0.8) + 
    scale_color_manual(values= total_col) +
    scale_fill_manual(values= total_col) +
    stat_ellipse(aes_string(x = "PC1", y = "PC2", color = "Treatment2"), linetype = 2) +
    annotate("text", x=0, y=0.2, label= paste("P=", p_val), size=2) +
    annotate("text", x=0, y=0.21, label= paste("R2=", round(r_sq, digits=3)), size=2) +
    labs(x=paste("PC1 (", round(var_exp$Relative_eig[1], digits=3)*100, "%)", sep=""), y= paste("PC2 (", round(var_exp$Relative_eig[2], digits=3)*100, "%)", sep="")) +
    geom_point(data=centroids, aes(x=PC1, y=PC2, fill=Treatment2), shape = 23, size=3, alpha=0.8, show.legend = F)
  
  
  
  plot_list[[names(I_C)[i]]] <- plot1
  
  #sink(paste(main_fp, "sil_widths.txt", sep='/'), append=T)
  #print(names(I_C)[i])
  #print(sil_list)
  #sink()
  
  sink(paste(main_fp, "treatment_pairwise_centroids.txt", sep="/"), append =T)
  print(pair_ad_list)
  sink()
}
day6_all_pcoa <- plot_grid(plot_list[[1]], plot_list[[2]],ncol=2)


######################################################################
#Plot all minor-member taxa to see if there is a pattern of the minor members 
#that explains the BMD/TJ overlap
x <- as.matrix(t(CLR_otutable)) #612 by 536

#Keep only OTUs in these samples
keep_t <- otu_table3
ileum_6 <- rownames(mapping)[mapping$Treatment2 %in% c("TJPbx", "BMD", "GroGel", "No_Inoc", "FMB11") & mapping$Tissue == "Ileum" & mapping$Collection_Day == "D06"]
keep_t <- keep_t[,colnames(keep_t) %in% ileum_6 ]
keep_t <- keep_t[rowSums(keep_t)<200 & rowSums(keep_t) > 5,] #334 OTUs

x <- x[rownames(keep_t),colnames(keep_t)]
rownames(x) <- taxonomy[rownames(x), 2] #334 x 536

keeps <- intersect(Ileum, Days[[2]])
x <- data.frame(t(x[,colnames(x) %in% keeps]))
these_rows <- rownames(x)

# replace rownames with unique ones
otu_names <- as.character(colnames(x))
taxa <- c()
for (i in 1:length(otu_names)){
  genus <- strsplit(otu_names[i], "..", fixed=T)[[1]][6]
  family <- strsplit(otu_names[i], "..", fixed=T)[[1]][5]
  order <- strsplit(otu_names[i], "..", fixed=T)[[1]][4]
  class <- strsplit(otu_names[i], "..", fixed=T)[[1]][3]
  phylum <- strsplit(otu_names[i], "..", fixed=T)[[1]][2]
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
colnames(x) <- taxa$taxa

#aggregate same taxa
x2 <- t(x)
names <- rownames(x2)

x3 <- data.frame(x2)
x3$taxa <- names
x4 <- aggregate(. ~ taxa, x3, FUN = mean)
x <- t(x4)
colnames(x) <- x[1,]
x <- x[2:nrow(x),]
x <- data.frame(x)
x <- data.frame(sapply(x, function(xx) as.numeric(as.character(xx))))
rownames(x) <- these_rows

y <- mapping[rownames(x),]

x$SampleID <- rownames(x)
x$Day <- mapping[rownames(x),"Collection_Day"]
x$Treatment <- mapping[rownames(x), "Treatment2"]
new_x <- melt(x, id.vars=c("SampleID", "Day", "Treatment"))
new_x <- new_x[new_x$Treatment %in% c("BMD", "TJPbx", "GroGel", "No_Inoc", "FMB11"),]
new_x2 <-aggregate(new_x, by=list(new_x$Treatment, new_x$variable), 
                 FUN=mean, na.rm=TRUE)


heat3 <- function (df, Xvar = names(df)[[1]], Yvar = names(df)[[2]],
                   fill = names(df)[[3]], LB = names(df)[[4]], star = NULL, p.adj.threshold = 1, 
                   association.threshold = 0, step = 0.2, colours = c("#FEBF2C",  "#C3003A",  "#551B44"), limits = NULL, legend.text = "", 
                   order.rows = TRUE, order.cols = TRUE, text.size = 10, filter.significant = TRUE, 
                   star.size = NULL, plot.values = FALSE) 
{
  if (is.null(limits)) {
    maxval <- max(abs(df[[fill]]))
    if (maxval <= 1) {
      limits <- c(-1, 1)
      limits <- c(min(df[[fill]]), max(df[[fill]]))
    }
    else {
      xmaxval <- ceiling(maxval)
      limits <- c(-maxval, maxval)
      limits <- range(df[[fill]])
    }
  }
  df[df[[fill]] < limits[[1]], fill] <- limits[[1]]
  df[df[[fill]] > limits[[2]], fill] <- limits[[2]]
  if (nrow(df) == 0) {
    warning("Input data frame is empty.")
    return(NULL)
  }
  if (filter.significant & !is.null(star)) {
    keep.X <- as.character(unique(df[((df[[star]] < p.adj.threshold) & 
                                        (abs(df[[fill]]) > association.threshold)), Xvar]))
    keep.Y <- as.character(unique(df[((df[[star]] < p.adj.threshold) & 
                                        (abs(df[[fill]]) > association.threshold)), Yvar]))
    df <- df[((df[[Xvar]] %in% keep.X) & (df[[Yvar]] %in% 
                                            keep.Y)), ]
  }
  theme_set(theme_bw(text.size))
  if (any(c("XXXX", "YYYY", "ffff") %in% names(df))) {
    stop("XXXX, YYYY, ffff are not allowed in df")
  }
  df[[Xvar]] <- factor(df[[Xvar]])
  df[[Yvar]] <- factor(df[[Yvar]])
  if (is.logical(order.rows) || is.logical(order.cols)) {
    rnams <- unique(as.character(df[[Xvar]]))
    cnams <- unique(as.character(df[[Yvar]]))
    mat <- matrix(0, nrow = length(rnams), ncol = length(cnams))
    rownames(mat) <- rnams
    colnames(mat) <- cnams
    for (i in 1:nrow(df)) {
      mat[as.character(df[i, Xvar]), as.character(df[i, 
                                                     Yvar])] <- df[i, fill]
    }
    mat <- t(mat)
    cind <- 1:ncol(mat)
    rind <- 1:nrow(mat)
  }
  if (is.logical(order.rows)) {
    if (order.rows) {
      if (nrow(mat) > 1 && ncol(mat) > 1) {
        rind <- hclust(as.dist(1 - cor(t(mat), use = "pairwise.complete.obs")))$order
      }
      if (nrow(mat) > 1 && ncol(mat) == 1) {
        rind <- order(mat[, 1])
      }
      order.rows <- rownames(mat)[rind]
    }
    else {
      order.rows <- rownames(mat)[rind]
    }
  }
  if (is.logical(order.cols)) {
    if (order.cols) {
      if (ncol(mat) > 1 && nrow(mat) > 1) {
        cind <- hclust(as.dist(1 - cor(mat, use = "pairwise.complete.obs")))$order
      }
      else {
        cind <- order(mat[1, ])
      }
      order.cols <- colnames(mat)[cind]
    }
    else {
      order.cols <- colnames(mat)[cind]
    }
  }
  df[[Yvar]] <- factor(df[[Yvar]], levels = order.rows)
  df[[Xvar]] <- factor(df[[Xvar]], levels = order.cols)
  XXXX <- YYYY <- ffff <- NULL
  df[["XXXX"]] <- df[[Xvar]]
  df[["YYYY"]] <- df[[Yvar]]
  df[["ffff"]] <- df[[fill]]
  p <- ggplot(df, aes(x = XXXX, y = YYYY, fill = ffff))
  p <- p + geom_tile()
  p <- p + scale_fill_gradientn(legend.text, breaks = seq(from = min(limits), 
                                                          to = max(limits), by = step), colours = colours, limits = limits)
  p <- p + xlab("") + ylab("")
  p <- p + theme(axis.text.x = element_text(angle = 90))
  p <- p + facet_grid(.~df[[LB]], scales="free", space = "free") +
    theme(axis.text.x = element_blank())
  if (!is.null(star)) {
    inds <- which((df[[star]] < p.adj.threshold) & (abs(df[[fill]]) > 
                                                      association.threshold))
    if (!is.null(star) & length(inds) > 0) {
      df.sub <- df[inds, ]
      if (is.null(star.size)) {
        star.size <- max(1, floor(text.size/2))
      }
      p <- p + geom_text(data = df.sub, aes(x = XXXX, 
                                            y = YYYY, label = "+"), col = "white", size = star.size)
    }
  }
  if (plot.values) {
    p <- p + geom_text(aes(label = round(ffff, 2)), size = 3)
  }
  p
}

new_x <- new_x[new_x$SampleID %in% beta_keeps,]
new_x$Treatment <- factor(new_x$Treatment, levels=c("No_Inoc", "BMD", "TJPbx", "FMB11", "GroGel"))

minor_ileum <- heat3(new_x, "SampleID", "variable", "value", "Treatment")
minor_ileum2 <- heat3(new_x2, "Group.1", "Group.2", "value", "Group.1")
######################################################################
#Plot together for figure 3
Fig_3 <- plot_grid(day6_all_pcoa, minor_ileum, nrow=2, rel_heights = c(1,1.7))

pdf(paste(main_fp,"/Figure3.pdf", sep=""), width=8.4, height=7.4)
Fig_3
dev.off()

######################################################################
#Supplemental Figure 3 - ceca and ileum all time points by treatment
otu_test <- read.table("data/Turkey_10k_Taxa.txt", 
                       sep="\t", 
                       comment="", 
                       header=T, 
                       skip=1, 
                       as.is=T, 
                       check.names=F, row=1)
tree = read_tree_greengenes('data/97_otus.tree')

m2 <- mapping[mapping$Tissue == "Ileum" | mapping$Tissue == "Ceca",]
beta_keeps <- intersect(rownames(m[Days[[2]],]), colnames(otu_test))

total_col <- c("#D95F02", "#7570B3", "#E7298A", "#1B9E77", "#E6AB02")
# "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"
plot_list <- list()
I_C <- list(Bodysites[[2]], Bodysites[[3]])
names(I_C) <- c("Ileum", "Ceca")


for(i in 1:length(I_C)){
  body <- I_C[[i]]
  for(k in 1:length(Days)){
    day <- Days[[k]]
    samples <- intersect(body,day)
    keeps <- intersect(rownames(m2), colnames(otu_test))
    keeps <- intersect(keeps, samples)
    otu_test1 <- otu_test[,keeps]
    otu_test1 <- otu_test1[rowSums(otu_test1 > 0) > 1, ]
    m <- m2[keeps,]
    #generate phyloseq objects
    OTU <- otu_table(otu_test1, taxa_are_rows=T)
    sampledata = sample_data(m)
    #manditory input for UniFrac
    physeq1 = phyloseq(OTU, sampledata, tree)
    #distances
    uni1 <- UniFrac(physeq1, weighted=FALSE, normalized=TRUE)
    PCOA <- pcoa(uni1)$vectors
    var_exp <- pcoa(uni1)$values
    
    #Run stats for diff. centroids
    beta_dist = uni1
    ad = adonis(beta_dist ~ m[,"Treatment2"], data=m, permutations=999)
    p_val <- ad$aov.tab[1,6]
    r_sq <- ad$aov.tab[1,5]
    #Run Stats for diff. dispersion
    beta_out <- betadisper(beta_dist, m$Treatment2)
    p_val_disp <- permutest(beta_out)$tab[1, 6]
    
    PCOA <- pcoa(uni1)$vectors
    for(c in 1:ncol(PCOA)){
      colnames(PCOA)[c] <- paste("PC",c, sep="")
    }
    PCOA <- data.frame(cbind(PCOA, rownames(PCOA)))
    colnames(PCOA)[ncol(PCOA)] <- "SampleID"
      
    PCOA <- merge(PCOA, m, by="SampleID")
    PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
    PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
    PCOA <- as.data.frame(PCOA)
    PCOA$Treatment2 <- factor(PCOA$Treatment2, levels=c("BMD", "TJPbx", "FMB11", "No_Inoc", "GroGel"))
    
    #centroids <- aggregate(cbind(PCOA$PC1,PCOA$PC2) ~ PCOA$Treatment2,PCOA,mean)
    #colnames(centroids) <- c("Treatment2", "PC1", "PC2")
    
    #pairwise adonis for differences:
    # pair_ad_list <- list()
    # for(t in 1:4){
    #   for(u in (t+1):5){
    #     t1 <- unique(m$Treatment2)[[t]]
    #     t2 <- unique(m$Treatment2)[[u]]
    #     sams <- rownames(m[m$Treatment2 == t1 | m$Treatment2 == t2,])
    #     beta_dist = as.dist(as.matrix(uni1)[sams,sams])
    #     ad = adonis(beta_dist ~ m[sams,"Treatment2"], data=m, permutations=999)
    #     p_val <- ad$aov.tab[1,6]
    #     r_sq <- ad$aov.tab[1,5]
    #     pair_ad_list[[paste(t1, "_", t2, "_", names(I_C)[i], "_", names(Days)[k], sep="")]] <- c("pval", p_val, "r_sq", r_sq)
    #   }
    # }
    
    plot1 <- ggplot(PCOA) +
      geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Treatment2"), alpha =0.7) + 
      scale_color_manual(values= total_col) +
      scale_fill_manual(values= total_col) +
      stat_ellipse(aes_string(x = "PC1", y = "PC2", color = "Treatment2"), linetype = 2) +
      annotate("text", x=0, y=0.2, label= paste("P=", p_val), size=2) +
      annotate("text", x=0, y=0.21, label= paste("R2=", round(r_sq, digits=3)), size=2) +
      labs(x=paste("PC1 (", round(var_exp$Relative_eig[1], digits=3)*100, "%)", sep=""), y= paste("PC2 (", round(var_exp$Relative_eig[2], digits=3)*100, "%)", sep="")) +
      guides(color=F)  
    #geom_point(data=centroids, aes(x=PC1, y=PC2, fill=Treatment2), shape = 23, size=3, alpha=0.8, show.legend = F)
    
    name <- paste(names(Days[k]), names(I_C[i]), sep="_")
    plot_list[[name]] <- plot1
    
    sink(paste(main_fp, "treatment_pairwise_centroids.txt", sep="/"), append =T)
    print(pair_ad_list)
    sink()
  }
}
plot3by3 <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],plot_list[[5]],plot_list[[6]],
                        labels=c(names(plot_list)), vjust=-0.5, hjust=-0.1,ncol = 3)

plot_this <- paste(main_fp, "Supplemental_Figure3.pdf", sep='/')
pdf(plot_this, width=7, height=4)
plot3by3
dev.off()
