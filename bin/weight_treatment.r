# Test weight gain differents by treatment. Account for tissue and day

weight_map <- mapping[Days[[3]],]
colnames(weight_map)[7] <- 3
colnames(weight_map)[8] <- 6
colnames(weight_map)[10] <- 13
weight_map$SampleID <- rownames(weight_map)

weight_map2 <- melt(weight_map, id.vars = c('SampleID', 'Treatment2'), measure.vars = c('3','6','13'))
weight_map2$variable <- as.numeric(as.character(weight_map2$variable))
weight_map2$value <- as.numeric(as.character(weight_map2$value))

colnames(weight_map2)[3] <- "Day"
colnames(weight_map2)[4] <- "Weight"

figure <- ggplot(weight_map2, aes(x=Day, y=Weight, color=Treatment2, group=Treatment2)) +
          geom_jitter(width = 0.25) +
          geom_smooth(se=FALSE) +
          theme_bw() +
          theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 10), 
            legend.title = element_blank(),
            legend.key.size = unit(0.2, "in"),
            legend.text = element_text(size=5),
            legend.position = 'bottom',
            axis.text = element_text(size=5),
            axis.title = element_text(size=8)) + 
          scale_color_manual(values=cols2(5))
fp <- paste(main_fp, "weight", sep="/")
dir.create(fp)
fp <- paste(fp, "/weight_time.pdf", sep="")
pdf(fp, height=4,width=6,useDingbats=FALSE)
print(figure)
dev.off()


#Compare the day 13 weights

weight_map3 <- weight_map2[weight_map2$Day == 13,]

figure2 <- ggplot(weight_map3, aes(x=Treatment2, y=Weight, color=Treatment2)) +
          geom_boxplot() +
          geom_jitter(width=0.15, alpha=0.75) +
          theme_bw() +
          theme(panel.background = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 10), 
            legend.title = element_blank(),
            legend.key.size = unit(0.2, "in"),
            legend.text = element_text(size=5),
            legend.position = 'bottom',
            axis.text = element_text(size=5),
            axis.title = element_text(size=8)) + 
         scale_color_manual(values=cols2(5)) 
fp <- paste(main_fp, "weight", sep="/")
fp <- paste(fp, "/day13.pdf", sep="")
pdf(fp, height=4,width=6,useDingbats=FALSE)
print(figure2)
dev.off()

fp <- paste(main_fp, "weight", "day13.txt", sep="/")
sink(fp)
kruskal.test(weight_map3$Treatment2 ~ weight_map3$Weight)
sink()

#Compare the day 3 weights
weight_map <- mapping[Days[[1]],]
colnames(weight_map)[7] <- 3
weight_map$SampleID <- rownames(weight_map)

weight_map2 <- melt(weight_map, id.vars = c('SampleID', 'Treatment2'), measure.vars = c('3'))
weight_map2$variable <- as.numeric(as.character(weight_map2$variable))
weight_map2$value <- as.numeric(as.character(weight_map2$value))

colnames(weight_map2)[3] <- "Day"
colnames(weight_map2)[4] <- "Weight"

figure2 <- ggplot(weight_map2, aes(x=Treatment2, y=Weight, color=Treatment2)) +
  geom_boxplot() +
  geom_jitter(width=0.15, alpha=0.75) +
  theme_bw() +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10), 
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "in"),
        legend.text = element_text(size=5),
        legend.position = 'bottom',
        axis.text = element_text(size=5),
        axis.title = element_text(size=8)) + 
  scale_color_manual(values=cols2(5)) 
fp <- paste(main_fp, "weight", sep="/")
fp <- paste(fp, "/day3.pdf", sep="")
pdf(fp, height=4,width=6,useDingbats=FALSE)
print(figure2)
dev.off()
fp <- paste(main_fp, "weight", "day3.txt", sep="/")
stats_out <- kruskal.test(weight_map2$Treatment2 ~ weight_map2$Weight)
sink(fp)
print(stats_out)
sink()

#Compare the day 6 weights
weight_map <- mapping[Days[[2]],]
colnames(weight_map)[7] <- 3
weight_map$SampleID <- rownames(weight_map)

weight_map2 <- melt(weight_map, id.vars = c('SampleID', 'Treatment2'), measure.vars = c('3'))
weight_map2$variable <- as.numeric(as.character(weight_map2$variable))
weight_map2$value <- as.numeric(as.character(weight_map2$value))

colnames(weight_map2)[3] <- "Day"
colnames(weight_map2)[4] <- "Weight"

figure2 <- ggplot(weight_map2, aes(x=Treatment2, y=Weight, color=Treatment2)) +
  geom_boxplot() +
  geom_jitter(width=0.15, alpha=0.75) +
  theme_bw() +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10), 
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "in"),
        legend.text = element_text(size=5),
        legend.position = 'bottom',
        axis.text = element_text(size=5),
        axis.title = element_text(size=8)) + 
  scale_color_manual(values=cols2(5)) 
fp <- paste(main_fp, "weight", sep="/")
fp <- paste(fp, "/day6.pdf", sep="")
pdf(fp, height=4,width=6,useDingbats=FALSE)
print(figure2)
dev.off()

fp <- paste(main_fp, "weight", "day6.txt", sep="/")
stats_out <- kruskal.test(weight_map2$Treatment2 ~ weight_map2$Weight)

sink(fp)
print(stats_out)
sink()


##Pairwise comparison

pvals<- c()
all_treats <- list(BMD, FMB11, GroGel, NoInoc, TJPbx)
for(i in 1:(length(all_treats)-1)){
  for(j in 2:length(all_treats)){
    test1 <- all_treats[[i]]
    test2 <- all_treats[[j]]
    all <- c(test1, test2)
    test_map <- weight_map2[weight_map2$SampleID %in% all,]
    stats_out <- t.test(test_map$Weight ~ test_map$Treatment2)
    pvals <- c(pvals, stats_out$p.value)
  }
}
