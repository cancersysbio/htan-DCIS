### plot_supplementary_figure3c.R #################################################################
# create scatterplot of NMF diagnostics

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(tidyr)
library(NMF)

# set working directory
setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

# set date
date <- Sys.Date()
### MAIN ##########################################################################################
basedir <- '/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/NMF/rna_nmf'
tnmf_file <- file.path(basedir, '216samples_with_ribosomal', '2020-12-17_nmf_results_with_ribosomal.RData')
rnmf_file <- file.path(basedir, 'RAHBT_265samples', '2021-01-04_nmf_results_RAHBT.RData')
# load NMF results 
load(tnmf_file)
tres_dcis <- res_dcis
load(rnmf_file)
rres_dcis <- res_dcis

# extract summary metrics 
tnmf_summary <- summary(tres_dcis)
rnmf_summary <- summary(rres_dcis)
# create plot data
tplot_data <- gather(
	tnmf_summary[,c('rank','cophenetic','silhouette.coef')],
	key = "key",
	value = "value",
	-rank
	)
tplot_data$cohort <- 'TBCRC'
rplot_data <- gather(
	rnmf_summary[,c('rank','cophenetic','silhouette.coef')],
	key = "key",
	value = "value",
	-rank
	)
rplot_data$cohort <- 'RAHBT'
col <- rep('black', nrow(tnmf_summary))
col[2] <- 'firebrick'

plot_data <- rbind(tplot_data, rplot_data)
# update key names 
plot_data$key <- gsub('.coef','', plot_data$key)
plot_data$group <- paste(plot_data$cohort, plot_data$key)

# create scatterplot 
create.scatterplot(
	value ~ rank | group,
	data = plot_data,
	ylimits = c(0,1.05),
	yaxis.tck = 0,
	xaxis.tck = 0,
	xat = seq(2,10,2),
	xaxis.cex = 1,
	yaxis.cex = 1,
	yat = seq(0,1,0.2),
	ylab.label = 'Diagnostics',
	xlab.label = 'NMF Rank',
	col = col,
	filename = file.path('figures', paste0('Supplementary_Figure_3c_nmf_diagnostic_scatterplot.png')),
	resolution = 300
	)