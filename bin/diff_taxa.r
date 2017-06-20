#### Find and plot any differentiated taxa
source("bin/diff.test.r")

#set output dir
diff_dir <- paste(main_fp, "diff_taxa/antibiotics", sep='/')

#Make a table with taxonomy as the rownames
taxonomy3 <- taxonomy[rownames(otutable4),]
working_table <- as.matrix(otutable4)
rownames(working_table) <- taxonomy3$V2

#Set up tests to run
test.ixs <- list(Antibiotics, NoInoc)
names(test.ixs) <- c("Antibiotic", "NoInoc")

pvals<- c()
ALPHA <- 0.25

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
        set2 <- intersect(union1, test.y[[1]])
        if(length(set1) > 2 && length(set2) > 2){
          full_set <- c(set1, set2)
          #keep taxa and the samples you are testing
          test_table <- t(working_table[,full_set,drop=F])
          #Keep taxa that have at least one count
          test_table <- test_table[,colSums(test_table)>0, drop=F]
          map_test <- mapping[full_set,]
          difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
          
          if(any(difftest$qvalues <= ALPHA)){
            signif.ix <- which(difftest$qvalues <= ALPHA)
            signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
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
              qval <- difftest$qvalues[signif.ix[k]]
              name <- paste(names(Bodysite), names(Day), names(test.x), names(test.y), names(signif.ix)[k], sep="-")
              fp_name <- paste(diff_dir, name, sep="/")
              
              #Stats output  
              sink(paste(fp_name, "txt", sep="."))
              cat(paste('q=',qval,' taxon: ',names(signif.ix)[k],'\n',sep=''))
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
diff_dir <- paste(main_fp, "diff_taxa/FMB11", sep='/')

#Set up tests to run
test.ixs <- list(FMB11, GroGel)
names(test.ixs) <- c("FMB11", "GroGel")

pvals<- c()
ALPHA <- 0.25

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
        set2 <- intersect(union1, test.y[[1]])
        if(length(set1) > 2 && length(set2) > 2){
          full_set <- c(set1, set2)
          #keep taxa and the samples you are testing
          test_table <- t(working_table[,full_set,drop=F])
          #Keep taxa that have at least one count
          test_table <- test_table[,colSums(test_table)>0, drop=F]
          map_test <- mapping[full_set,]
          difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
          
          if(any(difftest$qvalues <= ALPHA)){
            signif.ix <- which(difftest$qvalues <= ALPHA)
            signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
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
              qval <- difftest$qvalues[k]
              name <- paste(names(Bodysite), names(Day), names(test.x), names(test.y), names(signif.ix)[k], sep="-")
              fp_name <- paste(diff_dir, name, sep="/")
              
              #Stats output  
              sink(paste(fp_name, "txt", sep="."))
              cat(paste('q=',qval,' taxon: ',names(signif.ix)[k],'\n',sep=''))
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
ALPHA <- 0.25

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
        set2 <- intersect(union1, test.y[[1]])
        if(length(set1) > 2 && length(set2) > 2){
          full_set <- c(set1, set2)
          #keep taxa and the samples you are testing
          test_table <- t(working_table[,full_set,drop=F])
          #Keep taxa that have at least one count
          test_table <- test_table[,colSums(test_table)>0, drop=F]
          map_test <- mapping[full_set,]
          difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
          
          if(any(difftest$qvalues <= ALPHA)){
            signif.ix <- which(difftest$qvalues <= ALPHA)
            signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
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
              qval <- difftest$qvalues[signif.ix[k]]
              name <- paste(names(Bodysite), names(Day), names(test.x), names(test.y), names(signif.ix)[k], sep="-")
              fp_name <- paste(diff_dir, name, sep="/")
              
              #Stats output  
              sink(paste(fp_name, "txt", sep="."))
              cat(paste('q=',qval,' taxon: ',names(signif.ix)[k],'\n',sep=''))
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
