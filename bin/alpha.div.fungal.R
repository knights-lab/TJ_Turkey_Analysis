#### Get alpha diversity, plot and run stats

##to compare alphas##
alpha_func2 <- function(alpha_test, control_list, test_set, alpha_dir, color_by){
  control <- control_list[[1]]
  con_name <- names(control_list[1])
  for(i in 1:length(test_set)){
    pvals <- c()
    pbx_test <- test_set[[i]]
    pbx <- names(test_set[i])
    for(m in 1:length(Days_f)){
      day_test <- Days_f[[m]]
      day <- names(Days_f[m])
      subset <- pbx_test[pbx_test %in% day_test]
      subset_c <- control[control %in% day_test]
      name <- paste(day, pbx, "control", sep="-")
      print(name)
      full_set <- c(subset, subset_c)
      alpha_subset <- alpha_div2[,alpha_test, drop=FALSE]
      wtest <- wilcox.test(alpha_subset[subset,], alpha_subset[subset_c,])
      pvals <- c(pvals,wtest$p.value)
        
      sink(file_name, append =TRUE)
      cat(sprintf('\n%s, %s, vs. %s:\n', day, pbx, con_name))
      print(wtest)
      sink()
        
      #assign pdf name for plot
      name <- paste(day, pbx, con_name, sep="-")
      name <- paste(name, ".pdf", sep='')
      file_path <- paste(alpha_dir, name, sep='')
      pdf(file_path, height=4,width=6);
        
      #make alpha div box plots
      map <- f_map[full_set,]
      alpha_subset <- alpha_subset[rownames(map),,drop=TRUE]
       title <- sprintf('%s, %s vs. %s:', day, pbx, con_name)
        
      boxplot(alpha_subset ~ f_map[full_set,color_by,drop=TRUE],
                xlab='',ylab=alpha_test, main=title,
                col=cols2(2))
      dev.off()
    }
}
  #print fdr corrected pvals to stats file
  print(pvals)
  fdr.pvals <- p.adjust(pvals, method="fdr")
  sink(file_name, append =TRUE)
  cat("\nfdr adjusted pvals")
  print(fdr.pvals)
  sink()
}


#### Test PBX vs Control, account for body site and day####
color_by <- "Treatment"
#set output directory
alpha_dir <- paste(main_fp, "alpha_div/pbx_con_fungal/", sep='/')
#make stats file
file_name <- paste(alpha_dir, "Alpha_Stats.txt", sep='')
sink(file_name)
sink()
###Pbx vs control###
test_set <- Pbx_f
control_list <- list(GroGel_f)
names(control_list) <- c("GroGel")
alpha_test <- "shannon"

alpha_func2(alpha_test, control_list, test_set, alpha_dir, color_by)

#### Test ABX vs Control, account for body site and day####
#set output directory
alpha_dir <- paste(main_fp, "alpha_div/abx_con_fungal/", sep='/')
#make stats file
file_name <- paste(alpha_dir, "Alpha_Stats.txt", sep='')
sink(file_name)
sink()
###BMD vs control###
test_set <- list(BMD_f)
names(test_set) <- c("BMD")
control_list <- list(NoInoc_f)
names(control_list) <- c("NoInoc")

alpha_func2(alpha_test, control_list, test_set, alpha_dir, color_by)

#### Test Pbx vs Abx, account for body site and day####
color_by <- "Treatment"
#set output directory
alpha_dir <- paste(main_fp, "alpha_div/pbx_abx_fungal/", sep='/')
#make stats file
file_name <- paste(alpha_dir, "Alpha_Stats.txt", sep='')
sink(file_name)
sink()

test_set <- Pbx_f
control_list <- list(BMD_f)
names(control_list) <- c("BMD")

alpha_func2(alpha_test, control_list, test_set, alpha_dir, color_by)


###Test abx in general or pbx in general against general controls###
#set output directory
alpha_dir <- paste(main_fp, "alpha_div/general_abx_pbx_con_fungal/", sep='/')
#make stats file
file_name <- paste(alpha_dir, "Alpha_Stats.txt", sep='')
sink(file_name)
sink()

color_by <- "Abx"
test_set <- list(Antibiotics_f)
names(test_set) <- c("Antibiotics")
control_list <- list(Controls_f)
names(control_list) <- c("Controls")

alpha_func2(alpha_test, control_list, test_set, alpha_dir, color_by)

test_set <- list(Probiotics_f)
names(test_set) <- c("Probiotics")
color_by <- "Probiotic"

alpha_func2(alpha_test, control_list, test_set, alpha_dir, color_by)

####Plot Alpha diversity over time for each bodysite, color by treatment

#Add alphas to the map
alpha_div2 <- alpha_div2[rownames(f_map),]

f_map$shannon <- alpha_div2$shannon
f_map$observed_species <- alpha_div2$observed_species
f_map$simpson <- alpha_div2$simpson
alpha_metrics <- c("shannon", "observed_species", "simpson")

for(a in 1:length(alpha_metrics)){
  alpha_plots  <- c()
  alpha_use <- alpha_metrics[a]
  working_alpha <- melt(f_map, id.vars = c('SampleID', 'Treatment', 'Collection_Day'), measure.vars = c(alpha_use))
  working_alpha$Collection_Day[working_alpha$Collection_Day == "D03"] <- 3
  working_alpha$Collection_Day[working_alpha$Collection_Day == "D06"] <- 6
  working_alpha$Collection_Day[working_alpha$Collection_Day == "D13"] <- 13
  working_alpha$Collection_Day <- as.numeric(as.character(working_alpha$Collection_Day))
  figure <- ggplot(working_alpha, aes(x=Collection_Day, y=value, color=Treatment, group=Treatment)) +
      geom_jitter(width = 0.25) +
      geom_smooth(se=FALSE) +
      scale_color_manual(values=cols2(5))
  alpha_name <- paste("alpha_time_", alpha_use, ".pdf", sep='')
  plot_this <- paste(main_fp, "alpha_div/time_fungal", alpha_name, sep='/')
  save_plot(plot_this, figure, base_aspect_ratio = 1.8)
}


