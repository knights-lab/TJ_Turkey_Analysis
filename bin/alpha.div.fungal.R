######################################################################
#Add alpha diversity to map, and it will be way easier to plot
alpha_div2 <- alpha_div2[rownames(f_map),]
f_map$shannon <- alpha_div2[,"shannon"]
f_map$simpson <- alpha_div2[,"simpson"]
#f_map$pdwholetree <- alpha_div2[,"PD_whole_tree"]
f_map$obsspecies <- alpha_div2[,"observed_species"]

######################################################################
#One plot per day. Save as compound figure, one per metric
#x=treatment
alpha_dir <- paste(main_fp, "alpha_div/treatment_fungal/", sep='/')
#alpha_metrics <- c("shannon", "simpson", "pdwholetree", "obsspecies")
alpha_metrics <- c("shannon", "simpson", "obsspecies")

for(a in 1:length(alpha_metrics)){
  a_metric <- alpha_metrics[a]
  total_plot <- c()
  max_div <- max(f_map[,a_metric])
  min_div <- min(f_map[,a_metric])
  for(d in 1:length(Days_f)){
    plot1 <- ggplot(f_map[Days_f[[d]],]) +
      geom_boxplot(aes_string(x="Treatment2", y=a_metric, fill="Treatment2")) +
      scale_fill_manual(values=cols2(5)) +
      guides(fill=F) +
      expand_limits(y=c(min_div,max_div))
    total_plot[[names(Days_f)[d]]] <- plot1
  }
  total_plot2 <- plot_grid(total_plot[["D03"]],total_plot[["D06"]],total_plot[["D13"]], ncol=3)
  file_name <- paste(alpha_dir, a_metric, "_by_treatment.pdf", sep="")
  save_plot(file_name, total_plot2, col=3, base_aspect_ratio = 3)
}

######################################################################
#One plot per day. Save as compound figure, one per metric
#x=day
alpha_dir <- paste(main_fp, "alpha_div/time_fungal/", sep='/')

for(a in 1:length(alpha_metrics)){
  a_metric <- alpha_metrics[a]
  total_plot <- c()
  working_map <- f_map
  working_map$Collection_Day[working_map$Collection_Day == "D03"] <- 3
  working_map$Collection_Day[working_map$Collection_Day == "D06"] <- 6
  working_map$Collection_Day[working_map$Collection_Day == "D13"] <- 13
  working_map$Collection_Day <- as.numeric(working_map$Collection_Day)
  plot1 <- ggplot(working_map, aes_string(x="Collection_Day", y=a_metric,color="Treatment2", group="Treatment2")) +
    geom_jitter(width = 0.25) +
    geom_smooth(se=FALSE) +  
    scale_color_manual(values=cols2(5))
  file_name <- paste(alpha_dir, a_metric, "_by_time.pdf", sep="")
  save_plot(file_name, plot1, base_aspect_ratio = 2)
}




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
      #pdf(file_path, height=4,width=6);
        
      #make alpha div box plots
      map <- f_map[full_set,]
      alpha_subset <- alpha_subset[rownames(map),,drop=TRUE]
       title <- sprintf('%s, %s vs. %s:', day, pbx, con_name)
        
      #boxplot(alpha_subset ~ f_map[full_set,color_by,drop=TRUE],
      #          xlab='',ylab=alpha_test, main=title,
      #          col=cols2(2))
      #dev.off()
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


#### Test PBX vs Control, account for day####
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

