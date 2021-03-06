#Calculate Beta Div

beta_func2 <- function(beta_div, control, test_set, beta_dir){
  for(i in 1:length(test_set)){
    pvals <- c()
    pbx_test <- test_set[[i]]
    pbx <- names(test_set[i])
    for(m in 1:length(Days)){
      day_test <- Days[[m]]
      day <- names(Days[m])
      set1 <- pbx_test[pbx_test %in% day_test]
      set2 <- control[control %in% day_test]
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
      cat(sprintf('\n%s,%s:\n',pbx, day))
      print(ktest)
        
      #write stats to file
      sink(file_name, append =TRUE)
      cat(sprintf('\n%s, %s:\n',pbx,day))
      print(ktest)
      sink()
        
      #assign pdf name for plot
      name1 <- paste(pbx, day, sep="-")
      name1 <- paste(name1, ".pdf", sep='')
      file_path <- paste(beta_dir, name1, sep='')
      pdf(file_path, height=4,width=6)
        
      #make beta div box plots
      title <- sprintf('%s, %s',pbx, day)
      boxplot(sets_test,
                xlab='',ylab='Bray_Curtis', main=title,
                col=cols2(3))
        dev.off()
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

#### Test PBX vs Control, account for body site and day####
#set output directory
beta_dir <- paste(main_fp, "beta_div/BrayCurtis/pbx_con_fungal/", sep='/')
dir.create(beta_dir)
#make stats file
file_name <- paste(beta_dir, "Beta_Stats.txt", sep='')
sink(file_name)
sink()

##Set variables
test_set <- Pbx_f
names(test_set) <- c("FMB11", "TJPbx")
control_list <- list(GroGel_f)
names(control_list) <- c("GroGel")

#Run
beta_func2(fbray, control_list[[1]], test_set, beta_dir)

#### Test PBX vs Control, account for body site and day####
#set output directory
beta_dir <- paste(main_fp, "beta_div/BrayCurtis/abx_con_fungal/", sep='/')
dir.create(beta_dir)
#make stats file
file_name <- paste(beta_dir, "Beta_Stats.txt", sep='')
sink(file_name)
sink()

##Set variables
test_set <- list(BMD_f)
names(test_set) <- c("BMD")
control_list <- list(NoInoc_f)
names(control_list) <- c("NoInoc")
#Run
beta_func2(fbray, control_list[[1]], test_set, beta_dir)
