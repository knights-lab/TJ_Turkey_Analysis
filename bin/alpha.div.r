#### Get alpha diversity, plot and run stats

##to compare alphas##
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

####Plot Alpha diversity over time for each bodysite, color by treatment

#Add alphas to the map
alpha_div <- alpha_div[rownames(mapping),]

mapping$shannon <- alpha_div$shannon
mapping$observed_species <- alpha_div$observed_species
mapping$chao1 <- alpha_div$chao1
mapping$simpson <- alpha_div$simpson
mapping$pd_whole_tree <- alpha_div$PD_whole_tree
alpha_metrics <- c("shannon", "observed_species", "chao1", "simpson", "pd_whole_tree")

alpha_trachea <- mapping[Bodysites[[2]],]
alpha_ileum <- mapping[Bodysites[[1]],]
alpha_cecum <- mapping[Bodysites[[3]],]
alpha_tissues <- list(alpha_trachea, alpha_ileum, alpha_cecum)
names(alpha_tissues) <- c("Trachea", "Ileum", "Cecum")

for(a in 1:length(alpha_metrics)){
  alpha_plots  <- c()
  alpha_use <- alpha_metrics[a]
  for(t in 1:length(alpha_tissues)){
    working_alpha <- melt(alpha_tissues[t], id.vars = c('SampleID', 'Treatment2', 'Collection_Day'), measure.vars = c(alpha_use))
    working_alpha$Collection_Day[working_alpha$Collection_Day == "D03"] <- 3
    working_alpha$Collection_Day[working_alpha$Collection_Day == "D06"] <- 6
    working_alpha$Collection_Day[working_alpha$Collection_Day == "D13"] <- 13
    working_alpha$Collection_Day <- as.numeric(as.character(working_alpha$Collection_Day))
    figure <- ggplot(working_alpha, aes(x=Collection_Day, y=value, color=Treatment2, group=Treatment2)) +
      geom_jitter(width = 0.25) +
      geom_smooth(se=FALSE) +
      scale_color_manual(values=cols2(5))
    name <- names(alpha_tissues[t])
    alpha_plots[[name]] <- figure
  }
  plot3 <- plot_grid(alpha_plots[[1]], alpha_plots[[2]], alpha_plots[[3]],
                     labels=c(names(alpha_plots)),hjust=-2, ncol=3)
  alpha_name <- paste("alpha_time_", alpha_use, ".pdf", sep='')
  plot_this <- paste(main_fp, "alpha_div/time", alpha_name, sep='/')
  save_plot(plot_this, plot3,
            ncol = 3,
            nrow = 1,
            base_aspect_ratio = 1.8)
}

