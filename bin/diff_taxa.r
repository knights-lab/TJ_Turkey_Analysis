#### Find and plot any differentiated taxa
source("bin/diff.test.r")

test.otu.features<-function(otu, response, sig.level)
{
  pvals <- apply(otu, 2, function(feature) 
    (kruskal.test(feature~response, data.frame(feature=feature, response=response)))$p.value)
  adj.pvals <- p.adjust(pvals, "fdr")
  
  diff.features <- names(adj.pvals)[adj.pvals <= sig.level & !is.na(adj.pvals)]
  list(features=diff.features, pvals=adj.pvals)
}

#set output dir
diff_dir <- paste(main_fp, "diff_taxa/antibiotics", sep='/')

#Make a table with taxonomy as the rownames
taxonomy3 <- taxonomy[colnames(CLR_otutable),]
working_table <- as.matrix(t(CLR_otutable))
rownames(working_table) <- taxonomy3$V2

#Set up tests to run
test.ixs <- list(Antibiotics, NoInoc)
names(test.ixs) <- c("Antibiotic", "NoInoc")

pvals<- c()
ALPHA <- 0.05

#For each bodysite and and timepoint run tests
for(i in 1:length(Bodysites)){
  Bodysite <- Bodysites[i]
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
        if(length(set1) > 2 && length(set2) > 2){
          full_set <- c(set1, set2)
          #keep taxa and the samples you are testing
          test_table <- t(working_table[,full_set,drop=F])
          #Keep taxa that have a mean rel abundance of at least 0.01
          test_table <- test_table[,colMeans(test_table)>0.009, drop=F]
          map_test <- mapping[rownames(test_table),]
          difftest <- test.otu.features(test_table, response=map_test$Treatment2, sig.level = 0.10)
          #difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
          
          if(any(difftest$pvals <= ALPHA)){
            signif.ix <- which(difftest$pvals <= ALPHA)
            signif.ix <- signif.ix[order(difftest$pvals[signif.ix])]
            for(k in 1:length(signif.ix)){
            #  if(!is.null(difftest$norm.test.pvals)){
            #    norm.test <- difftest$norm.test.pvals[k]
            #  } else {
            #    norm.test <- '0'
            #  }
            #  if(norm.test < 0.05){
            #    qval <- difftest$qvalues[k]
            #  } else {
            #    
              taxon <- gsub(";", "", names(signif.ix)[k])
              qval <- difftest$pvals[signif.ix[[k]]]
              name <- paste(names(Bodysite), names(Day), names(test.x), names(test.y), taxon, sep="-")
              fp_name <- paste(diff_dir, name, sep="/")
              
              #Stats output  
              sink(paste(fp_name, "txt", sep="."))
              cat(paste('q=',qval,' taxon: ',taxon,'\n',sep=''))
              sink()
              
              #boxplots
              pdf(paste(fp_name, "pdf", sep="."),width=4,height=4)
              boxplot(test_table[,names(signif.ix)[k]] ~ map_test$Treatment2, 
                      xlab='', ylab="Relative Abundance", main=name,
                      col=cols2(length(unique(map_test$Treatment2))))
              dev.off()
            }
          } else {
            cat("not significant.")
          }
        top_50_ix <- unlist(sort.int(difftest$pvals, index.return=TRUE)[[2]])[1:50]
        top_50_taxa <- colnames(test_table)[top_50_ix]
        } else {
          cat("Less than two samples in one group, skipping this test.")
        }
      }
    }
  }
}

top_50_OTUs <- as.character(droplevels(taxonomy3[which(taxonomy3$V2 %in% top_50_taxa), "V1"]))
  
#set output dir
diff_dir <- paste(main_fp, "diff_taxa/FMB11", sep='/')

#Set up tests to run
test.ixs <- list(FMB11, GroGel)
names(test.ixs) <- c("FMB11", "GroGel")

pvals<- c()
ALPHA <- 0.05

#For each bodysite and and timepoint run tests
for(i in 1:length(Bodysites)){
  Bodysite <- Bodysites[i]
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
        if(length(set1) > 2 && length(set2) > 2){
          full_set <- c(set1, set2)
          #keep taxa and the samples you are testing
          test_table <- t(working_table[,full_set,drop=F])
          #Keep taxa that have at least one count
          test_table <- test_table[,colMeans(test_table)>0.009, drop=F]
          map_test <- mapping[rownames(test_table),]
          difftest <- test.otu.features(test_table, response=map_test$Treatment2, sig.level = 0.10)
          #difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
          
          if(any(difftest$pvals <= ALPHA)){
            signif.ix <- which(difftest$pvals <= ALPHA)
            signif.ix <- signif.ix[order(difftest$pvals[signif.ix])]
            for(k in 1:length(signif.ix)){
            #  if(!is.null(difftest$norm.test.pvals)){
            #    norm.test <- difftest$norm.test.pvals[k]
            #  } else {
            #    norm.test <- '0'
            #  }
            #  if(norm.test < 0.05){
            #    qval <- difftest$qvalues[k]
            #  } else {
            #
              taxon <- gsub(";", "", names(signif.ix)[k])
              qval <- difftest$pvals[signif.ix[[k]]]
              name <- paste(names(Bodysite), names(Day), names(test.x), names(test.y), taxon, sep="-")
              fp_name <- paste(diff_dir, name, sep="/")
              
              #Stats output  
              sink(paste(fp_name, "txt", sep="."))
              cat(paste('q=',qval,' taxon: ',taxon,'\n',sep=''))
              sink()
              
              #boxplots
              pdf(paste(fp_name, "pdf", sep="."),width=4,height=4)
              boxplot(test_table[,names(signif.ix)[k]] ~ map_test$Treatment2, 
                      xlab='', ylab="Relative Abundance", main=name,
                      col=cols2(length(unique(map_test$Treatment2))))
              dev.off()
            }
          } else {
            cat("not significant.")
            }
        } else {
          cat("Less than two samples in one group, skipping this test.")
        }
      }
    }
  }
}

#set output dir
diff_dir <- paste(main_fp, "diff_taxa/TJPbx", sep='/')

#Set up tests to run
test.ixs <- list(TJPbx, GroGel)
names(test.ixs) <- c("TJPbx", "GroGel")

pvals<- c()
ALPHA <- 0.05

#For each bodysite and and timepoint run tests
for(i in 1:length(Bodysites)){
  Bodysite <- Bodysites[i]
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
        if(length(set1) > 2 && length(set2) > 2){
          full_set <- c(set1, set2)
          #keep taxa and the samples you are testing
          test_table <- t(working_table[,full_set,drop=F])
          #Keep taxa that have at least one count
          test_table <- test_table[,colMeans(test_table)>0.009, drop=F]
          map_test <- mapping[rownames(test_table),]
          difftest <- test.otu.features(test_table, response=map_test$Treatment2, sig.level = 0.10)
          #difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
          
          if(any(difftest$pvals <= ALPHA)){
            signif.ix <- which(difftest$pvals <= ALPHA)
            signif.ix <- signif.ix[order(difftest$pvals[signif.ix])]
            for(k in 1:length(signif.ix)){
            #  if(!is.null(difftest$norm.test.pvals)){
            #    norm.test <- difftest$norm.test.pvals[k]
            #  } else {
            #    norm.test <- '0'
            #  }
            #  if(norm.test < 0.05){
            #    qval <- difftest$qvalues[k]
            #  } else {
            #    
              taxon <- gsub(";", "", names(signif.ix)[k])
              qval <- difftest$pvals[signif.ix[[k]]]
              name <- paste(names(Bodysite), names(Day), names(test.x), names(test.y), taxon, sep="-")
              fp_name <- paste(diff_dir, name, sep="/")
              
              #Stats output  
              sink(paste(fp_name, "txt", sep="."))
              cat(paste('q=',qval,' taxon: ',taxon,'\n',sep=''))
              sink()
              
              #boxplots
              pdf(paste(fp_name, "pdf", sep="."),width=4,height=4)
              boxplot(test_table[,names(signif.ix)[k]] ~ map_test$Treatment2, 
                      xlab='', ylab="Relative Abundance", main=name,
                      col=cols2(length(unique(map_test$Treatment2))))
              dev.off()
            }
          } else {
            cat("not significant.")
            }
        } else {
          cat("Less than two samples in one group, skipping this test.")
        }
      }
    }
  }
}

