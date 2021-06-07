### create_figure3g.R ##############################################################################
# create figure 3g - MIBI HER2, Ki67, ER and GLUT1 associated with RNA subtypes

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(tidyr)

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

# set date
date <- Sys.Date()
### FIND SAMPLES WITH HER2 AMP ####################################################################
find_samples_with_her2_amp <- function() {
	# set gistic dir 
	gisticdir <- '/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/GISTIC2/run_brlen98_t30_armpeel_gene_res01_conf95'
	# read in gistic peaks 
	peaks <- read.delim(
		file.path(gisticdir, 'all_lesions.conf_95.txt'),
		as.is = TRUE,
		check.names = FALSE
		)
	peaks <- peaks[which(peaks$'Residual q values after removing segments shared with higher peaks' < 0.01),]
	peaks <- peaks[which(peaks$'Amplitude Threshold' != 'Actual Copy Change Given'),]
	# extract her2 CNAs
	her2 <- unlist(peaks[peaks$Descriptor == '17q12   ',-c(1:9, ncol(peaks))])
	her2 <- data.frame(
		DNA_sample = names(her2),
		HER2 = her2
		)
	return(her2)
	}
### MAIN ##########################################################################################
# read in supplementary table 
mibi <- read.csv(
	'../data/201125_PFT_Publication.csv', 
	header = TRUE
	)
mibi <- mibi[,c('TMAD_Patient','SUBTYPE','ER_subtype_status','Status_tumor_ER_freq','HER2_subtype_status','Status_tumor_HER2intense_freq', 'Status_tumor_Ki67_freq','Status_tumor_GLUT1_freq')]
colnames(mibi) <- gsub('TMAD_Patient','Donor', colnames(mibi))

# read in rahbt design file 
design <- read.delim(
	'../NMF/2021-02-10_RAHBT_design_RNA_CNA_icluster.tsv',
	as.is = TRUE
	)
#design <- design[design$NEW_Clustering_RAHBT_021121 == 1,]

# read in gistic 
her2 <- find_samples_with_her2_amp()

# create plot data 
plot_data <- merge(
	design[,c('DNA_sample','Donor','RNA_3')],
	mibi,
	by = 'Donor'
	)
# reformat 
plot_data <- gather(
	plot_data[,-c(4,5,7)],
	value = 'Freq',
	key = 'Status',
	Status_tumor_ER_freq,
	Status_tumor_HER2intense_freq,
	Status_tumor_Ki67_freq,
	Status_tumor_GLUT1_freq
	)
plot_data$Status <- gsub('Status_tumor_|_freq|intense','',plot_data$Status)
plot_data$RNA_3 <- factor(plot_data$RNA_3)

# add her2 status 
plot_data <- merge(
	plot_data,
	her2,
	by = 'DNA_sample',
	all.x = TRUE
	)

# change point colour according to HER2 amplification
subset <- unique(plot_data[,c(1,2,6)])
col <- rep('darkgrey', nrow(subset))
col[which(is.na(subset$HER2))] <- 'grey75'
col[which(subset$HER2 == 1)] <- 'firebrick2'
col[which(subset$HER2 == 2)] <- 'firebrick4'

# create boxplot 
create.boxplot(
	Freq ~ RNA_3 | Status,
	data = plot_data,
	add.stripplot = TRUE,
	ylimit = c(0,1.05),
	yat = seq(0, 1,0.5),
	ylab.axis.padding = 2,
	points.cex = 0.8, 
	xaxis.lab = c(
		expression(bold('ER'['low'])),
		expression(bold('quiescent')),
		expression(bold('ER'['high']))
		),
	xaxis.cex = 1.2,
	xaxis.rot = 45,
	ylab.cex = 1.8,
	xlab.cex = 1.8,
	strip.cex = 1.8,
	points.col = col,
	layout = c(2,2),
	yaxis.cex = 1.2,
	#xaxis.cex = 1.2,
	top.padding = 3,
	#right.padding = 5,
	ylab.label = 'Frequency of tumor cells',
	xlab.label = 'RNA clusters',
	filename = 'figures/Figure2f_MIBI_RNA_clusters_boxplot_coloured.svg',
	#height = 5,
	#width = 8,
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 1.5,
                             fill = c('darkgrey','firebrick2','firebrick4')
                             ),
                         text = list(
                             lab = c('None', 'Low','High'),
                             cex = 1.5
                             ),
                         cex = 1,
                         title = 'ERBB2 Gain'
                         )
                     ),
                 x = 0.68,
                 y = 0.92,
                 corner = c(0,1),
                 draw = FALSE
                 )
             ),
	resolution = 300
	)

