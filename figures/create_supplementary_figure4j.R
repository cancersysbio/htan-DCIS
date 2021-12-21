### plot_supplementary_figures_3i.R ##############################################################
# concordance between RNA 3 and CNA 6 clusters

### PREAMBLE ######################################################################################
library(vcd)

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

# set date
date <- Sys.Date()
### MAIN ##########################################################################################
# read in tdesign 
tdesign <- read.delim(
	file.path('..', 'NMF', '2021-02-02_TBCRC_design_RNA_CNA_icluster.tsv'),
	as.is = TRUE
	)
rdesign <- read.delim(
	file.path('..', 'NMF', '2021-02-10_RAHBT_design_RNA_CNA_icluster.tsv'),
	as.is = TRUE
	)
design <- rbind(tdesign[,c('DNA_sample','RNA_3','CNA_6')], rdesign[,c('DNA_sample','RNA_3','CNA_6')])
design <- design[!is.na(design$DNA_sample),]
colnames(design) <- c('sample','RNA','CNA')

# create mosaic plot
png(paste0(date, '_RNA3_CNA6_concordance_mosaic.png'), res = 300, width = 7, height = 7, units = 'in')
mosaic(~RNA+CNA, design, shade=TRUE)
dev.off()