#Calculate Beta Div
beta_div <- as.matrix(vegdist(t(taxa_table), method = "bray"))

beta_func <- function(beta_div, control, test_set, beta_dir){
  for(i in 1:length(test_set)){
    pvals <- c()
    pbx_test <- test_set[[i]]
    pbx <- names(test_set[i])
    for(k in 1:length(Bodysites)){
      body_test <- Bodysites[[k]]
      body <- names(Bodysites[k])
      for(m in 1:length(Days)){
        day_test <- Days[[m]]
        day <- names(Days[m])
        set1 <- pbx_test[pbx_test %in% body_test]
        set2 <- control[control %in% body_test]
        set1 <- set1[set1 %in% day_test]
        set2 <- set2[set2 %in% day_test]
        full_set <- c(set1, set2)
        ktest <- c()
        set1_within <- c()
        set2_within <- c()
        between_sets <- c()
        for(c in 1:length(colnames(beta_div))){
          for(r in (c+1):length(rownames(beta_div))){
            if(rownames(beta_div)[r] %in% set1 && colnames(beta_div)[c] %in% set1){
              set1_within <- c(set1_within, beta_div[r,c])
            } else {
              if(rownames(beta_div)[r] %in% set1 && colnames(beta_div)[c] %in% set2){
                between_sets <- c(between_sets, beta_div[r,c])
              } else {
                if(rownames(beta_div)[r] %in% set2 && colnames(beta_div)[c] %in% set1){
                  between_sets <- c(between_sets, beta_div[r,c])
                } else {
                  if(rownames(beta_div)[r] %in% set2 && colnames(beta_div)[c] %in% set2){
                    set2_within <- c(set2_within, beta_div[r,c])
                  }
                }
              }
            }
          }
        }
        sets_test <- c()
        sets_test <- list(set2_within, between_sets, set1_within)
        names(sets_test) <- c(pbx, paste(pbx, "control", sep='-'), "control") 
        ktest <- kruskal.test(sets_test)
        pvals <- c(pvals,ktest$p.value)
        
        #print stats to screen
        cat(sprintf('\n%s,%s, %s:\n',pbx, body, day))
        print(ktest)
        
        #write stats to file
        sink(file_name, append =TRUE)
        cat(sprintf('\n%s, %s, %s:\n',pbx, body, day))
        print(ktest)
        sink()
        
        #assign pdf name for plot
        name1 <- paste(pbx, body, day, sep="-")
        name1 <- paste(name1, ".pdf", sep='')
        file_path <- paste(beta_dir, name1, sep='')
        pdf(file_path, height=4,width=6)
        
        #make beta div box plots
        title <- sprintf('%s, %s, %s',pbx, body, day)
        boxplot(sets_test,
                xlab='',ylab='Bray Curtis', main=title,
                col=cols2(3))
        dev.off()
      }
    }
  }
  #print fdr corrected pvals to stats file
  print(pvals)
  fdr.pvals <- p.adjust(pvals, method="fdr")
  sink(file_name, append =TRUE)
  cat("\nfdr adjusted pvals:")
  print(fdr.pvals)
  sink()
  print(fdr.pvals)
}

#### Test Treatment vs Control, account for body site and day####

#set output directory
beta_dir <- paste(main_fp, "beta_div/BrayCurtis/one/", sep='/')
#make stats file
file_name <- paste(beta_dir, "Beta_Stats.txt", sep='')
sink(file_name)
sink()

###Pbx vs control###
test_set <- Pbx
control <- GroGel

beta_func(beta_div, control, test_set, beta_dir)

#set output directory
beta_dir <- paste(main_fp, "beta_div/BrayCurtis/two/", sep='/')
#make stats file
file_name <- paste(beta_dir, "Beta_Stats.txt", sep='')
sink(file_name)
sink()

###Pbx vs control###
test_set <- list(BMD)
names(test_set) <- c("BMD")
control <- NoInoc

beta_func(beta_div, control, test_set, beta_dir)

