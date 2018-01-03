####Run all the commands####

source("bin/load.r")

source("bin/alpha.div.r")

source("bin/pcoa.r")

source("bin/beta.div.r")

source("bin/taxa.sum.r")

source("bin/diff_taxa.r")

source("bin/heat_map_difftaxa.r")

source("bin/transcript_PCOA.r")

source("bin/reactome.r")

source("bin/heat_map.r")

#source("bin/heat_map_pathways.r")

source("bin/load.fungal.r")

source("bin/alpha.div.fungal.R")

source("bin/pcoa.fungal.r")

source("bin/beta.div.fungal.r")

source("bin/taxa.sum.fungal.r")

source("bin/diff_taxa.fungal.r")

source("bin/weight_treatment.r")

save.image(file='Analysis_Env.RData')

