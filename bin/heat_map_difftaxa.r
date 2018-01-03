#Heatmap of diff. taxa

#Find top taxa for Ileum, All Days for abx vs. none
test.ixs <- list(Antibiotics, NoInoc)
names(test.ixs) <- c("Antibiotic", "NoInoc")

pvals<- c()
ALPHA <- 0.1

#For each bodysite and and timepoint run tests
heatmap_otus <- c()
Bodysite <- Bodysites[2]
for(j in 1:length(Days)){
  Day <- Days[j]
  union1 <- intersect(Bodysite[[1]], Day[[1]])
  for(n in 1:(length(test.ixs)-1)){
    for(m in (n+1):length(test.ixs)){
      test.x <- test.ixs[n]
      test.y <- test.ixs[m]
      set1 <- intersect(union1, test.x[[1]])
      set1 <- intersect(set1, colnames(working_table))
      set2 <- intersect(union1, test.y[[1]])
      set2 <- intersect(set2, colnames(working_table))
      full_set <- c(set1, set2)
      #keep taxa and the samples you are testing
      test_table <- t(working_table[,full_set,drop=F])
      #Keep taxa that have a mean rel abundance of at least 0.01
      test_table <- test_table[,colMeans(test_table)>0.009, drop=F]
      map_test <- mapping[rownames(test_table),]
      difftest <- test.otu.features(test_table, response=map_test$Treatment2, sig.level = 0.10)
      #difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
      if(any(difftest$pvals <= ALPHA)){
        yasss <- names(difftest$pvals<=ALPHA)
        heatmap_otus <- c(heatmap_otus, yasss)
      }
    }
  }
}
heatmap_otus <- unique(heatmap_otus)

x <- as.matrix(t(CLR_otutable))
rownames(x) <- taxonomy3$V2
x <- x[heatmap_otus,] #keep only the top taxa (from abx vs no inoc)

Cor_out2 <- cor(t(x))
Cor_out3 <- abs(Cor_out2) > 0.95
new_taxatable <- x

for(i in 1:nrow(Cor_out3)){
  for(j in i:nrow(Cor_out3)){
    if(i==j){
      next
    }
    if(is.na(Cor_out3[j,i])){
      next
    }
    if(Cor_out3[j,i]){
      new_taxatable[j,] <- new_taxatable[j,] + new_taxatable[i,]
      new_taxatable[i,] <- rep(0, ncol(new_taxatable))
    }
  }
}

#Collapse
taxa_table.new <- new_taxatable[rowSums(new_taxatable) > 0, ]
x <- taxa_table.new

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

x <- data.frame(t(x[,colnames(x) %in% Ileum])) #subset to just ileum
y <- mapping[rownames(x),]

dat <- NULL

for(i in colnames(x)){ # OTUs
  a<-as.numeric(x[,i])
  b <- y$Treatment2
  tmp <- c(i,summary(aov(a~b))[[1]][["Pr(>F)"]][1])
  if(is.null(dat)){
    dat <- tmp
  } else {
    dat<-rbind(dat,tmp)
  }
}

dat<- data.frame(row.names = NULL, dat, stringsAsFactors = FALSE)

colnames(dat)<-c("OTUs","Pvalue")
dat$Pvalue <- as.numeric(dat$Pvalue)
dat$fdr_pval <- (p.adjust(dat$Pvalue, 'fdr'))

dat$Significance<-cut(dat$fdr_pval, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) #converts numeric p-values to factors based on significance level

# prune for complete cases with significant correlations
dat.c<-dat[complete.cases(dat),]
taxkeep<-unique(dat.c$OTUs[dat.c$fdr_pval<1])

taxassign<-taxkeep
dat.w<-dat.c[dat.c$OTUs %in% taxassign,] 

dat.m <- dat.w

# order the OTUs
# make wide transcripts v. taxa filled with correlation

order <- dat.m %>% select(OTUs, Pvalue)
order <- order %>% spread(OTUs, Pvalue)
rownames(order) <- "Pvalue"
order <- order[,2:ncol(order)]

#can't reorder unless it's a factor first
dat.m$OTUs <- as.factor(dat.m$OTUs)

# reorder the factors in the dataframe to be in the taxaorder
dat.m$OTUs <- factor(dat.m$OTUs, levels(colnames(order)))

# get the order Samples

d3 <- intersect(Ileum, Days[[1]])
d3_b <- intersect(d3, BMD)
d3_n <- intersect(d3, NoInoc)
d3_g <- intersect(d3, GroGel)
d3_f <- intersect(d3, FMB11)
d3_t <- intersect(d3, TJPbx)

d6 <- intersect(Ileum, Days[[2]])
d6_b <- intersect(d6, BMD)
d6_n <- intersect(d6, NoInoc)
d6_g <- intersect(d6, GroGel)
d6_f <- intersect(d6, FMB11)
d6_t <- intersect(d6, TJPbx)

d13 <- intersect(Ileum, Days[[3]])
d13_b <- intersect(d13, BMD)
d13_n <- intersect(d13, NoInoc)
d13_g <- intersect(d13, GroGel)
d13_f <- intersect(d13, FMB11)
d13_t <- intersect(d13, TJPbx)

sample_order <- c(d3_b, d3_n, d3_g, d3_f, d3_t, d6_b, d6_n, d6_g, d6_f, d6_t, d13_b, d13_n, d13_g, d13_f, d13_t)
sample_order2 <- rep(sample_order, ncol(x))

# set transcripts as factor
x$SampleID <- rownames(x)
x$Day <- mapping[rownames(x),"Collection_Day"]
x$Treatment <- mapping[rownames(x), "Treatment2"]
new_x <- melt(x, id.vars=c("SampleID", "Day", "Treatment"))

## reorder the factors in the dataframe to be in the taxaorder
#new_x$SampleID <- factor(new_x$SampleID, levels(new_x$SampleID)[sample_order])
#new_x$variable <- factor(new_x$variable, levels(new_x$variable)[colnames(order)])
#new_x$variable <- factor(new_x$variable, levels=(new_x$variable)[order(new_x$value)])

#plot the correlations for the collapsed levels
# out_dir <- paste(main_fp, "diff_taxa", sep="/")
# name <- "Taxa_Heatmap2.pdf"
# fp <- paste(out_dir, name, sep="/")
# pdf(fp, height=8,width=12)
# ggplot(data = new_x, aes(x=reorder(SampleID, sample_order2), y=reorder(variable, value), fill=value)) + 
#   geom_tile(color = "white") +
#   #facet_wrap(~Day, scale="free_x")
#   facet_grid(~ Day + Treatment, scales = "free_x") +
#   scale_fill_gradient2(low="#5e3c99", mid="#f7f7f7", high="#e66101", limit = c(-2,10), midpoint = 4, space = "Lab") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 8, hjust = 1), 
#         axis.text.y = element_text(size = 8),
#         panel.spacing.x=unit(0, "lines")) +
#   labs(x = "Samples", y = "OTUs") +
#   theme(strip.text.y = element_text(angle = 0, face = "italic"), 
#         axis.text.x = element_blank(),
#         strip.text.x = element_text(size = 4))+
#   coord_fixed() #+
#   #geom_text(data = dat.m, aes(label=Significance), color="black", size=2)
# dev.off()
#dev.off()
pdf(paste(main_fp, "diff_taxa/heatmap_all.pdf", sep="/"), height=8, width=12)
heat(new_x, "SampleID", "variable", "value")
dev.off()

day3 <- new_x[new_x$Day=="D03",]

day6 <- new_x[new_x$Day=="D06",]

day13 <- new_x[new_x$Day=="D13",]



##############

heat2 <- function (df, Xvar = names(df)[[1]], Yvar = names(df)[[2]],
          fill = names(df)[[3]], LB = names(df)[[4]], star = NULL, p.adj.threshold = 1, 
          association.threshold = 0, step = 0.2, colours = c("purple4", 
                                                             "#5E3C99", "white", "#f5bf99", "#e66101"), limits = NULL, legend.text = "", 
          order.rows = TRUE, order.cols = TRUE, text.size = 10, filter.significant = TRUE, 
          star.size = NULL, plot.values = FALSE) 
{
  if (is.null(limits)) {
    maxval <- max(abs(df[[fill]]))
    if (maxval <= 1) {
      limits <- c(-1, 1)
    }
    else {
      xmaxval <- ceiling(maxval)
      limits <- c(-maxval, maxval)
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
  p <- p + facet_grid(.~df[[LB]], scales="free", space = "free")
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
out_dir <- paste(main_fp, "diff_taxa", sep="/")
name <- "Taxa_Heatmap_Day3.pdf"
fp <- paste(out_dir, name, sep="/")
pdf(fp, height=8,width=12)
heat2(day3, "SampleID", "variable", "value", "Treatment") 
dev.off()

name <- "Taxa_Heatmap_Day6.pdf"
fp <- paste(out_dir, name, sep="/")
pdf(fp, height=8,width=12)
heat2(day6, "SampleID", "variable", "value", "Treatment") 
dev.off()

name <- "Taxa_Heatmap_Day13.pdf"
fp <- paste(out_dir, name, sep="/")
pdf(fp, height=8,width=12)
heat2(day13, "SampleID", "variable", "value", "Treatment") 
dev.off()
