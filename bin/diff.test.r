# returns p-values, q-values, 
# indices of those below alpha, and renamed features with '*' etc.
# and subset of data if requested
# add.stars adds stars to the column names that are significant at .05, .01, .001, .0001
# if filename is provided, saves result of statistical test to file
# if parametric, uses t-test/ANOVA; else uses mann-whitney U/Kruskal-wallis
"differentiation.test" <- function (x,category, alpha=0.05, parametric=TRUE,
                                    include.subset=FALSE, correct.pvalues=TRUE){
  category <- as.factor(as.character(category))
  norm.test.pvals <- NULL
  if(length(unique(category)) < 2) stop('Category only has one level')
  if(parametric){
    pvals <- apply(x,2, function(taxon){
      if(var(taxon) == 0){
        NA
      } else {
        summary(lm(taxon~category))[[4]][2,4]
      }
    })
    stats <- apply(x,2, function(taxon){
      if(var(taxon) == 0){
        NA
      } else {
        summary(lm(taxon~category))[[4]][2,1]
      }
    })
    norm.test.pvals <- apply(x, 2, function(taxon){
      fit <- lm(taxon~category);
      ks.test(rstudent(fit), pnorm, mean=mean(rstudent(fit)), sd=sd(rstudent(fit)))$p.value
    })
  } else {
    if(length(levels(category)) == 2){
      ix1 <- category == levels(category)[1]
      ix2 <- category == levels(category)[2]
      pvals <- apply(x,2,function(taxon) wilcox.test(taxon[ix1], taxon[ix2],exact=FALSE)$p.value)
      stats <- apply(x,2,function(taxon) wilcox.test(taxon[ix1], taxon[ix2],exact=FALSE)$statistic)
    } else {
      pvals <- apply(x,2,function(taxon) kruskal.test(taxon ~ category)$p.value)
      stats <- apply(x,2,function(taxon) kruskal.test(taxon ~ category)$statistic)
    }
  }
  na.ix <- is.na(pvals)
  
  if(correct.pvalues){
    adj.pvals <- rep(NA,length(pvals))
    names(adj.pvals) <- names(pvals)
    adj.pvals[!na.ix] <- p.adjust(pvals[!na.ix],'fdr')
  } else {
    adj.pvals <- pvals
  }
  keep.ix <- adj.pvals < alpha
  keep.ix[is.na(keep.ix)] <- FALSE
  
  
  
  result <- list()
  
  # add stars to column names based on significance
  if(any(keep.ix)){
    annotations <- colnames(x)
    thresholds <- c(.05, .01, .001, .0001)
    for(i in seq_along(thresholds)){
      
      star.ix <- adj.pvals[!na.ix] <= thresholds[i]
      
      if(any(star.ix)){
        for(j in which(star.ix)){
          annotations[!na.ix][j] <- paste(annotations[!na.ix][j],'*',sep='')
        }
      }
    }
    result$annotations <- annotations
    result$features <- which(keep.ix)
  } else {
    result$features <- NULL
  }
  result$qvalues <- adj.pvals
  result$pvalues <- pvals
  result$stats <- stats
  result$norm.test.pvals <- norm.test.pvals
  
  # classwise means
  result$classwise.means <- t(apply(x,2,function(xx) sapply(split(xx,category),mean)))
  colnames(result$classwise.means) <- sprintf('%s mean',colnames(result$classwise.means))
  
  if(include.subset){
    result$subset <- x[,keep.ix,drop=F]
    colnames(result$subset) <- annotations[keep.ix]
  }
  
  return(result)
}


# saves list of results from differentiation.test to file (or prints)
"write.differentiation.test.results" <- function(results, filename='differentiated.features.txt'){
  if(!is.null(filename)){
    scipen.save <- options('scipen')
    options(scipen=20)
    hits <- cbind(results$pvalues, results$qvalues)
    hits <- cbind(hits, results$classwise.means)
    colnames(hits)[1:2] <- c('pvalue','qvalue')
    hits <- hits[!is.na(hits[,1]),]
    hits <- hits[order(hits[,1]),]
    sink(filename)
    cat('Feature\t')
    write.table(hits,quote=F,sep='\t')
    sink(NULL)
    options(scipen=scipen.save)
  }
}

test.otu.features<-function(otu, response, sig.level)
{
  pvals <- apply(otu, 2, function(feature) 
    (kruskal.test(feature~response, data.frame(feature=feature, response=response)))$p.value)
  adj.pvals <- p.adjust(pvals, "fdr")
  
  diff.features <- names(adj.pvals)[adj.pvals <= sig.level & !is.na(adj.pvals)]
  list(features=diff.features, pvals=adj.pvals)
}

test.otu.features2<-function(otu, response, sig.level)
{
  pvals <- apply(otu, 2, function(feature) 
    (t.test(feature~response, data.frame(feature=feature, response=response)))$p.value)
  adj.pvals <- p.adjust(pvals, "fdr")
  
  diff.features <- names(adj.pvals)[adj.pvals <= sig.level & !is.na(adj.pvals)]
  list(features=diff.features, pvals=adj.pvals)
}

# #set output dir
# diff_dir <- paste(main_fp, "diff_taxa/antibiotics", sep='/')
# 
# #Make a table with taxonomy as the rownames
# taxonomy3 <- taxonomy[colnames(CLR_otutable),]
# working_table <- as.matrix(t(CLR_otutable))
# rownames(working_table) <- taxonomy3$V2
# 
# #Set up tests to run
# test.ixs <- list(Antibiotics, NoInoc)
# names(test.ixs) <- c("Antibiotic", "NoInoc")
# 
# pvals<- c()
# ALPHA <- 0.05

# #For each bodysite and and timepoint run tests
# for(i in 1:length(Bodysites)){
#   Bodysite <- Bodysites[i]
#   for(j in 1:length(Days)){
#     Day <- Days[j]
#     union1 <- intersect(Bodysite[[1]], Day[[1]])
#     for(n in 1:(length(test.ixs)-1)){
#       for(m in (n+1):length(test.ixs)){
#         test.x <- test.ixs[n]
#         test.y <- test.ixs[m]
#         set1 <- intersect(union1, test.x[[1]])
#         set1 <- intersect(set1, colnames(working_table))
#         set2 <- intersect(union1, test.y[[1]])
#         set2 <- intersect(set2, colnames(working_table))
#         if(length(set1) > 2 && length(set2) > 2){
#           full_set <- c(set1, set2)
#           #keep taxa and the samples you are testing
#           test_table <- t(working_table[,full_set,drop=F])
#           #Keep taxa that have a mean rel abundance of at least 0.01
#           test_table <- test_table[,colMeans(test_table)>0.009, drop=F]
#           map_test <- mapping[rownames(test_table),]
#           difftest <- test.otu.features(test_table, response=map_test$Treatment2, sig.level = 0.10)
#           #difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
#           
#           if(any(difftest$pvals <= ALPHA)){
#             signif.ix <- which(difftest$pvals <= ALPHA)
#             signif.ix <- signif.ix[order(difftest$pvals[signif.ix])]
#             for(k in 1:length(signif.ix)){
#             #  if(!is.null(difftest$norm.test.pvals)){
#             #    norm.test <- difftest$norm.test.pvals[k]
#             #  } else {
#             #    norm.test <- '0'
#             #  }
#             #  if(norm.test < 0.05){
#             #    qval <- difftest$qvalues[k]
#             #  } else {
#             #    
#               taxon <- gsub(";", "", names(signif.ix)[k])
#               qval <- difftest$pvals[signif.ix[[k]]]
#               name <- paste(names(Bodysite), names(Day), names(test.x), names(test.y), taxon, sep="-")
#               fp_name <- paste(diff_dir, name, sep="/")
#               
#               #Stats output  
#               sink(paste(fp_name, "txt", sep="."))
#               cat(paste('q=',qval,' taxon: ',taxon,'\n',sep=''))
#               sink()
#               
#               #boxplots
#               pdf(paste(fp_name, "pdf", sep="."),width=4,height=4)
#               boxplot(test_table[,names(signif.ix)[k]] ~ map_test$Treatment2, 
#                       xlab='', ylab="Relative Abundance", main=name,
#                       col=cols2(length(unique(map_test$Treatment2))))
#               dev.off()
#             }
#           } else {
#             cat("not significant.")
#           }
#         top_50_ix <- unlist(sort.int(difftest$pvals, index.return=TRUE)[[2]])[1:50]
#         top_50_taxa <- colnames(test_table)[top_50_ix]
#         } else {
#           cat("Less than two samples in one group, skipping this test.")
#         }
#       }
#     }
#   }
# }
