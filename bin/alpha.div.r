#### Get alpha diversity, plot and run stats

######################################################################
#Add alpha diversity to map, and it will be way easier to plot
alpha_div <- alpha_div[rownames(mapping),]
mapping$shannon <- alpha_div[,"shannon"]
mapping$simpson <- alpha_div[,"simpson"]
mapping$pdwholetree <- alpha_div[,"PD_whole_tree"]
mapping$obsspecies <- alpha_div[,"observed_species"]

######################################################################
#One plot per bodysite, per day. Save as compound figure, one per metric
#x=treatment
alpha_dir <- paste(main_fp, "alpha_div/treatment/", sep='/')
alpha_metrics <- c("shannon", "simpson", "pdwholetree", "obsspecies")

for(a in 1:length(alpha_metrics)){
  a_metric <- alpha_metrics[a]
  total_plot <- c()
  for(s in 1:length(Bodysites)){
    bodysite_plot <- c()
    max_div <- max(mapping[Bodysites[[s]],a_metric])
    min_div <- min(mapping[Bodysites[[s]],a_metric])
    for(d in 1:length(Days)){
      samples_keep <- intersect(Bodysites[[s]], Days[[d]])
      plot1 <- ggplot(mapping[samples_keep,]) +
        geom_boxplot(aes_string(x="Treatment2", y=a_metric, fill="Treatment2")) +
        scale_fill_manual(values=cols2(5)) +
        guides(fill=F) +
        expand_limits(y=c(min_div,max_div))
      bodysite_plot[[names(Days[d])]] <- plot1
    }
    row_plot <- plot_grid(bodysite_plot[["D03"]],bodysite_plot[["D06"]],bodysite_plot[["D13"]], ncol=3)
    total_plot[[names(Bodysites[s])]] <- row_plot
  }
  file_name <- paste(alpha_dir, a_metric, "_by_treatment.pdf", sep="")
  final_plot <- plot_grid(total_plot[["Trachea"]], total_plot[["Ileum"]], total_plot[["Cecum"]], nrow=3)
  save_plot(file_name, final_plot,nrow=3, col=3, base_aspect_ratio = 3)
}

######################################################################
#One plot per bodysite, per day. Save as compound figure, one per metric
#x=day
alpha_dir <- paste(main_fp, "alpha_div/time/", sep='/')

for(a in 1:length(alpha_metrics)){
  a_metric <- alpha_metrics[a]
  total_plot <- c()
  max_div <- max(mapping[,a_metric])
  min_div <- min(mapping[,a_metric])
  for(s in 1:length(Bodysites)){
    bodysite_subset <- mapping[Bodysites[[s]],]
    bodysite_subset$Collection_Day[bodysite_subset$Collection_Day == "D03"] <- 3
    bodysite_subset$Collection_Day[bodysite_subset$Collection_Day == "D06"] <- 6
    bodysite_subset$Collection_Day[bodysite_subset$Collection_Day == "D13"] <- 13
    bodysite_subset$Collection_Day <- as.numeric(bodysite_subset$Collection_Day)
    plot1 <- ggplot(bodysite_subset, aes_string(x="Collection_Day", y=a_metric,color="Treatment2", group="Treatment2")) +
        geom_jitter(width = 0.25) +
        geom_smooth(se=FALSE) +  
        scale_color_manual(values=cols2(5)) +
        expand_limits(y=c(min_div, max_div))
    
   total_plot[[names(Bodysites[s])]] <- plot1
  }
  file_name <- paste(alpha_dir, a_metric, "_by_time.pdf", sep="")
  final_plot <- plot_grid(total_plot[["Trachea"]], total_plot[["Ileum"]], total_plot[["Cecum"]], ncol=3)
  save_plot(file_name, final_plot, col=3, base_aspect_ratio = 4)
}


######################################################################
#Pairwise alpha comparisons function
alpha_func <- function(alpha_test, control_list, test_set, alpha_dir, color_by){
  control <- control_list[[1]]
  con_name <- names(control_list[1])
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
        subset <- pbx_test[pbx_test %in% body_test]
        subset_c <- control[control %in% body_test]
        subset <- subset[subset %in% day_test]
        subset_c <- subset_c[subset_c %in% day_test]
        name <- paste(body, day, pbx, "control", sep="-")
        print(name)
        full_set <- c(subset, subset_c)
        alpha_subset <- alpha_div[,alpha_test, drop=FALSE]
        wtest <- wilcox.test(alpha_subset[subset,], alpha_subset[subset_c,])
        pvals <- c(pvals,wtest$p.value)
        
        sink(file_name, append =TRUE)
        cat(sprintf('\n%s, %s, %s, vs. %s:\n', body, day, pbx, con_name))
        print(wtest)
        sink()
        
        #assign pdf name for plot
        name <- paste(body, day, pbx, con_name, sep="-")
        name <- paste(name, ".pdf", sep='')
        file_path <- paste(alpha_dir, name, sep='')
        pdf(file_path, height=4,width=6);
        
        #make alpha div box plots
        map <- mapping[full_set,]
        alpha_subset <- alpha_subset[rownames(map),,drop=TRUE]
        title <- sprintf('%s, %s, %s vs. %s:',body, day, pbx, con_name)
        
        boxplot(alpha_subset ~ map[,color_by,drop=TRUE],
                xlab='',ylab=alpha_test, main=title,
                col=cols2(2))
        dev.off()
      }
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
alpha_dir <- paste(main_fp, "alpha_div/pbx_con/", sep='/')
#make stats file
file_name <- paste(alpha_dir, "Alpha_Stats.txt", sep='')
sink(file_name)
sink()
###Pbx vs control###
test_set <- Pbx
control_list <- list(GroGel)
names(control_list) <- c("GroGel")
alpha_test <- "PD_whole_tree"

alpha_func(alpha_test, control_list, test_set, alpha_dir, color_by)

#### Test ABX vs Control, account for body site and day####
#set output directory
alpha_dir <- paste(main_fp, "alpha_div/abx_con/", sep='/')
#make stats file
file_name <- paste(alpha_dir, "Alpha_Stats.txt", sep='')
sink(file_name)
sink()
###BMD vs control###
test_set <- list(BMD)
names(test_set) <- c("BMD")
control_list <- list(NoInoc)
names(control_list) <- c("NoInoc")

alpha_func(alpha_test, control_list, test_set, alpha_dir, color_by)

#### Test Pbx vs Abx, account for body site and day####
color_by <- "Treatment"
#set output directory
alpha_dir <- paste(main_fp, "alpha_div/pbx_abx/", sep='/')
#make stats file
file_name <- paste(alpha_dir, "Alpha_Stats.txt", sep='')
sink(file_name)
sink()

test_set <- Pbx
control_list <- list(BMD)
names(control_list) <- c("BMD")

alpha_func(alpha_test, control_list, test_set, alpha_dir, color_by)


###Test abx in general or pbx in general against general controls###
#set output directory
alpha_dir <- paste(main_fp, "alpha_div/general_abx_pbx_con/", sep='/')
#make stats file
file_name <- paste(alpha_dir, "Alpha_Stats.txt", sep='')
sink(file_name)
sink()

color_by <- "Abx"
test_set <- list(Antibiotics)
names(test_set) <- c("Antibiotics")
control_list <- list(Controls)
names(control_list) <- c("Controls")

alpha_func(alpha_test, control_list, test_set, alpha_dir, color_by)

test_set <- list(Probiotics)
names(test_set) <- c("Probiotics")
color_by <- "Probiotic"

alpha_func(alpha_test, control_list, test_set, alpha_dir, color_by)


