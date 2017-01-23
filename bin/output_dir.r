####Make output directories####

now<- format(Sys.time(), "%H-%M-%S")
main_fp <- paste("analysis",now, sep="-")
dir.create(main_fp)

dir.create(paste(main_fp, "alpha_div", sep='/'))
dir.create(paste(main_fp, "alpha_div/one", sep="/"))
dir.create(paste(main_fp, "alpha_div/two", sep="/"))
dir.create(paste(main_fp, "alpha_div/three", sep="/"))
dir.create(paste(main_fp, "alpha_div/four", sep="/"))

dir.create(paste(main_fp, "beta_div", sep='/'))

dir.create(paste(main_fp, "beta_div/BrayCurtis", sep='/'))
dir.create(paste(main_fp, "beta_div/BrayCurtis/one", sep="/"))
dir.create(paste(main_fp, "beta_div/BrayCurtis/two", sep="/"))
dir.create(paste(main_fp, "beta_div/BrayCurtis/three", sep="/"))
dir.create(paste(main_fp, "beta_div/BrayCurtis/pcoa", sep='/'))
dir.create(paste(main_fp, "beta_div/BrayCurtis/pcoa/Ileum", sep='/'))
dir.create(paste(main_fp, "beta_div/BrayCurtis/pcoa/Ceca", sep='/'))
dir.create(paste(main_fp, "beta_div/BrayCurtis/pcoa/Trachea", sep='/'))

dir.create(paste(main_fp, "taxa_sum", sep='/'))
dir.create(paste(main_fp, "taxa_sum/one", sep='/'))
dir.create(paste(main_fp, "taxa_sum/two", sep='/'))
dir.create(paste(main_fp, "taxa_sum/three", sep='/'))

dir.create(paste(main_fp, "diff_taxa", sep='/'))
dir.create(paste(main_fp, "diff_taxa/antibiotics", sep='/'))
dir.create(paste(main_fp, "diff_taxa/FMB11", sep='/'))
dir.create(paste(main_fp, "diff_taxa/TJPbx", sep='/'))




