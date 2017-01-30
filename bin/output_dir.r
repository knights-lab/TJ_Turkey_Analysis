####Make output directories####

now<- format(Sys.time(), "%H-%M-%S")
main_fp <- paste("analysis",now, sep="-")
dir.create(main_fp)

dir.create(paste(main_fp, "alpha_div", sep='/'))
dir.create(paste(main_fp, "alpha_div/pbx_con", sep="/"))
dir.create(paste(main_fp, "alpha_div/pbx_con_fungal", sep="/"))
dir.create(paste(main_fp, "alpha_div/pbx_abx", sep="/"))
dir.create(paste(main_fp, "alpha_div/pbx_abx_fungal", sep="/"))
dir.create(paste(main_fp, "alpha_div/abx_con", sep="/"))
dir.create(paste(main_fp, "alpha_div/abx_con_fungal", sep="/"))
dir.create(paste(main_fp, "alpha_div/general_abx_pbx_con", sep="/"))
dir.create(paste(main_fp, "alpha_div/general_abx_pbx_con_fungal", sep="/"))
dir.create(paste(main_fp, "alpha_div/time", sep="/"))
dir.create(paste(main_fp, "alpha_div/time_fungal", sep="/"))

dir.create(paste(main_fp, "beta_div", sep='/'))
dir.create(paste(main_fp, "beta_div/BrayCurtis", sep='/'))
dir.create(paste(main_fp, "beta_div/BrayCurtis/pbx_con_fungal", sep='/'))
dir.create(paste(main_fp, "beta_div/BrayCurtis/abx_con_fungal", sep='/'))
dir.create(paste(main_fp, "beta_div/Unifrac", sep='/'))
dir.create(paste(main_fp, "beta_div/Unifrac/pbx_con", sep='/'))
dir.create(paste(main_fp, "beta_div/Unifrac/abx_con", sep='/'))
dir.create(paste(main_fp, "beta_div/Unifrac", sep='/'))
dir.create(paste(main_fp, "beta_div/WUnifrac", sep='/'))
dir.create(paste(main_fp, "beta_div/PCoa_SiteDay", sep='/'))
dir.create(paste(main_fp, "beta_div/PCoa_Day_fungal", sep='/'))

dir.create(paste(main_fp, "taxa_sum", sep='/'))
dir.create(paste(main_fp, "taxa_sum_fungal", sep='/'))

dir.create(paste(main_fp, "diff_taxa", sep='/'))
dir.create(paste(main_fp, "diff_taxa/antibiotics", sep='/'))
dir.create(paste(main_fp, "diff_taxa/FMB11", sep='/'))
dir.create(paste(main_fp, "diff_taxa/TJPbx", sep='/'))




