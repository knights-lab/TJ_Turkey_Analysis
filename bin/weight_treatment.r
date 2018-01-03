# Test weight gain differents by treatment. Account for tissue and day

weight_map <- mapping[Days[[3]],]
weight_map$Treatment2 <- factor(weight_map$Treatment2, levels=c("No_Inoc", "GroGel", "TJPbx", "FMB11", "BMD"))

colnames(weight_map)[7] <- 3
colnames(weight_map)[8] <- 6
colnames(weight_map)[10] <- 13
weight_map$SampleID <- rownames(weight_map)

weight_map2 <- melt(weight_map, id.vars = c('SampleID', 'Treatment2', "Gain313", "Gain613"), measure.vars = c('3','6','13'))
weight_map2$variable <- as.numeric(as.character(weight_map2$variable))
weight_map2$value <- as.numeric(as.character(weight_map2$value))

colnames(weight_map2)[5] <- "Day"
colnames(weight_map2)[6] <- "Weight"

figure <- ggplot(weight_map2, aes(x=Day, y=Weight, color=Treatment2, group=Treatment2)) +
          geom_jitter(width = 0.25) +
          geom_smooth(se=FALSE) +
          scale_color_manual(values=cols2(5))
fp <- paste(main_fp, "weight", sep="/")
dir.create(fp)
fp <- paste(fp, "/weight_time.pdf", sep="")
save_plot(fp, figure, base_aspect_ratio = 1.8)


#Compare the day 13 weights
weight_map3 <- weight_map2[weight_map2$Day == 13,]
weight_map3$Gain313 <- as.numeric(weight_map3$Gain313)
weight_map3$Gain613 <- as.numeric(weight_map3$Gain613)

ggplot(weight_map3, aes(x=Treatment2, y=Gain313, color=Treatment2)) +
  geom_boxplot() +
  geom_jitter(width=0.15, alpha=0.75) +
  scale_color_manual(values=cols2(5))

ggplot(weight_map3, aes(x=Treatment2, y=Gain613, color=Treatment2)) +
  geom_boxplot() +
  geom_jitter(width=0.15, alpha=0.75) +
  scale_color_manual(values=cols2(5))

figure2 <- ggplot(weight_map3, aes(x=Treatment2, y=Weight, color=Treatment2)) +
            geom_boxplot() +
            geom_jitter(width=0.15, alpha=0.75) +
            scale_color_manual(values=cols2(5)) 
fp <- paste(main_fp, "weight", sep="/")
fp <- paste(fp, "/day13.pdf", sep="")
save_plot(fp, figure2, base_aspect_ratio = 1.8)

fp <- paste(main_fp, "weight", "day13.txt", sep="/")
sink(fp)
kruskal.test(weight_map3$Treatment2 ~ weight_map3$Weight)
sink()

#Compare the day 3 weights
weight_map <- mapping[Days[[1]],]
colnames(weight_map)[7] <- 3
weight_map$SampleID <- rownames(weight_map)
weight_map$Treatment2 <- factor(weight_map$Treatment2, levels=c("No_Inoc", "GroGel", "TJPbx", "FMB11", "BMD"))


weight_map2 <- melt(weight_map, id.vars = c('SampleID', 'Treatment2'), measure.vars = c('3'))
weight_map2$variable <- as.numeric(as.character(weight_map2$variable))
weight_map2$value <- as.numeric(as.character(weight_map2$value))

colnames(weight_map2)[3] <- "Day"
colnames(weight_map2)[4] <- "Weight"

figure2 <- ggplot(weight_map2, aes(x=Treatment2, y=Weight, color=Treatment2)) +
  geom_boxplot() +
  geom_jitter(width=0.15, alpha=0.75) +
  scale_color_manual(values=cols2(5)) 
fp <- paste(main_fp, "weight", sep="/")
fp <- paste(fp, "/day3.pdf", sep="")
save_plot(fp, figure2, base_aspect_ratio = 1.8)
fp <- paste(main_fp, "weight", "day3.txt", sep="/")
stats_out <- kruskal.test(weight_map2$Treatment2 ~ weight_map2$Weight)
sink(fp)
print(stats_out)
sink()

#Compare the day 6 weights
weight_map <- mapping[Days[[2]],]
colnames(weight_map)[7] <- 3
weight_map$SampleID <- rownames(weight_map)
weight_map$Treatment2 <- factor(weight_map$Treatment2, levels=c("No_Inoc", "GroGel", "TJPbx", "FMB11", "BMD"))


weight_map2 <- melt(weight_map, id.vars = c('SampleID', 'Treatment2'), measure.vars = c('3'))
weight_map2$variable <- as.numeric(as.character(weight_map2$variable))
weight_map2$value <- as.numeric(as.character(weight_map2$value))

colnames(weight_map2)[3] <- "Day"
colnames(weight_map2)[4] <- "Weight"

figure2 <- ggplot(weight_map2, aes(x=Treatment2, y=Weight, color=Treatment2)) +
  geom_boxplot() +
  geom_jitter(width=0.15, alpha=0.75) +
  scale_color_manual(values=cols2(5)) 
fp <- paste(main_fp, "weight", sep="/")
fp <- paste(fp, "/day6.pdf", sep="")
save_plot(fp, figure2, base_aspect_ratio = 1.8)

fp <- paste(main_fp, "weight", "day6.txt", sep="/")
stats_out <- kruskal.test(weight_map2$Treatment2 ~ weight_map2$Weight)

sink(fp)
print(stats_out)
sink()


# ##Pairwise comparison
# 
# pvals<- c()
# all_treats <- list(BMD, FMB11, GroGel, NoInoc, TJPbx)
# names(all_treats) <- c("BMD", "FMB11", "GroGel", "NoInoc", "TJPbx")
# for(i in 1:(length(all_treats)-1)){
#   for(j in (i+1):length(all_treats)){
#     test1 <- all_treats[[i]]
#     test2 <- all_treats[[j]]
#     all <- c(test1, test2)
#     test_map <- weight_map2[weight_map2$SampleID %in% all,]
#     stats_out <- t.test(test_map$Weight ~ test_map$Treatment2)
#     pvals <- c(pvals, stats_out$p.value)
#   }
# }

## Is diversity correlated with weight?
## Look by day and tissue
plot_list <- c()
weight_headers <- c("D3_weight", "D6_weight", "D13_weight")


for(i in 1:length(Bodysites)){
  body <- Bodysites[[i]]
  for(k in 1:length(Days)){
    day <- Days[[k]]
    samples <- intersect(body,day)
    alpha_div1 <- a_mapping[samples,]
    alpha_div1[,weight_headers[[k]]] <- as.numeric(alpha_div1[,weight_headers[k]])
    plot1 <- ggplot(alpha_div1, aes_string(x = weight_headers[k], y = "obsspecies")) +
      geom_point(size = 3, color=cols2(8)[8])
    name <- paste(names(Days[k]), names(Bodysites[i]), sep="_")
    plot_list[[name]] <- plot1
  }
}

plot3by3 <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],
                      labels=c(names(plot_list)), vjust=-0.5, hjust=-0.1,ncol = 3)
plot_this <- paste(main_fp, "weight/alpha_weight.pdf", sep='/')
save_plot(plot_this, plot3by3,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.8)


### Look at fungal diversity and weight

plot_list <- c()
weight_headers <- c("D3_weight", "D6_weight", "D13_weight")

#Add alphas to the map
alpha_div2 <- alpha_div2[rownames(f_map),]

f_map$shannon <- alpha_div2$shannon
f_map$observed_species <- alpha_div2$observed_species
f_map$simpson <- alpha_div2$simpson

for(k in 1:length(Days_f)){
  day <- Days_f[[k]]
  alpha_div1 <- f_map[day,]
  alpha_div1[,weight_headers[[k]]] <- as.numeric(alpha_div1[,weight_headers[k]])
  plot1 <- ggplot(alpha_div1, aes_string(x = weight_headers[k], y = "shannon")) +
    geom_point(size = 3, color=cols2(8)[8])
  name <- names(Days_f[k])
  plot_list[[name]] <- plot1
}
plot3 <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]],
                      labels=c(names(plot_list)), ncol = 3)
plot_this <- paste(main_fp, "weight/fungal_alpha_weight.pdf", sep='/')
save_plot(plot_this, plot3,
          ncol = 3,
          base_aspect_ratio = 1.8)
