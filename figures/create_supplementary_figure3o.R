### create_PI_EMT_boxplots.R ######################################################################
# create PI and EMT boxplots

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

# set date
date <- Sys.Date()
### MAIN ##########################################################################################
# read in PI data 
pi_df <- read.csv(
	'../data/PI_score_RAHBTn265_TBCRCn216_SHS011821.csv',
	as.is = TRUE,
	header = TRUE
	)
colnames(pi_df) <- c("Sample_ID", "PI")

# read in tdesign file 
tdesign <- read.delim(
	'../NMF/2021-02-02_TBCRC_design_RNA_CNA_icluster.tsv',
	as.is = TRUE
	)
rdesign <- read.delim(
	'../NMF/2021-02-12_RAHBT_design_RNA_CNA_icluster.tsv',
	as.is = TRUE
	)

# create plot data 
tplot_data <- merge(pi_df, tdesign[,c('Sample_ID','RNA_3')], by = 'Sample_ID')
tplot_data$cohort <- 'TBCRC'
rplot_data <- merge(pi_df, rdesign[,c('Sample_ID','RNA_3')], by = 'Sample_ID')
rplot_data$cohort <- 'RAHBT'
plot_data <- rbind(tplot_data, rplot_data)
# factor RNA groups
plot_data$RNA_3 <- factor(plot_data$RNA_3)

tstats_1vs23 <- wilcox.test(
	plot_data[plot_data$cohort == 'TBCRC' & plot_data$RNA_3 == 1, 'PI'],
	plot_data[plot_data$cohort == 'TBCRC' & plot_data$RNA_3 %in% c(2,3), 'PI']
	)
median1 <- median(plot_data[plot_data$cohort == 'TBCRC' & plot_data$RNA_3 == 1, 'PI'])
median23 <- median(plot_data[plot_data$cohort == 'TBCRC' & plot_data$RNA_3 %in% c(2,3), 'PI'])
fc_1vs23 <- median1-median23
t1vs23_pvalue_sci <- scientific.notation(tstats_1vs23$p.value, digits = 2, type = 'list');

tstats_2vs13 <- wilcox.test(
	plot_data[plot_data$cohort == 'TBCRC' & plot_data$RNA_3 == 2, 'PI'],
	plot_data[plot_data$cohort == 'TBCRC' & plot_data$RNA_3 %in% c(1,3), 'PI']
	)
median2 <- median(plot_data[plot_data$cohort == 'TBCRC' & plot_data$RNA_3 == 2, 'PI'])
median13 <- median(plot_data[plot_data$cohort == 'TBCRC' & plot_data$RNA_3 %in% c(1,3), 'PI'])
fc_2vs13 <- median2-median13
t2vs13_pvalue_sci <- scientific.notation(tstats_2vs13$p.value, digits = 2, type = 'list');

# create boxplot
create.boxplot(
	PI ~ RNA_3,
	data = plot_data[plot_data$cohort == 'TBCRC',],
	add.stripplot = TRUE,
	filename = file.path('figures', paste0(date, '_TBCRC_PI_RNA3_boxplot.png')),
	ylab.axis.padding = 2,
	strip.cex = 1.5,
	ylimits = c(0,6.5),
	yat = seq(0,6,2),
	main = 'TBCRC',
	main.cex = 1.8,
	xaxis.rot = 45,
	xaxis.lab = c(
		expression(bold('ER'['low'])),
		expression(bold('quiescent')),
		expression(bold('ER'['high']))
		),
	xaxis.cex = 1.5,
	ylab.cex = 1.8,
	xlab.cex = 1.8,
	yaxis.cex = 1.5,
	height = 6,
	width = 4,
	ylab.label = 'Proliferation Index',
	xlab.label = 'RNA Subtypes',
	key = list(
			text = list(
				lab = c(
					as.expression(substitute(
						fc[compare]*value,
						list(fc = "ES", compare = "ERlow", value = paste(':', round(fc_1vs23, digits = 1)))
						)),
					as.expression(substitute(
						pvalue[compare]*base %*% 10^exponent, 
						list(pvalue = "P-value", compare = "ERlow", base = paste(':', t1vs23_pvalue_sci[[1]]), exponent = t1vs23_pvalue_sci[[2]])
						)),
					as.expression(substitute(
						fc[compare]*value,
						list(fc = "ES", compare = "quiescent", value = paste(':', round(fc_2vs13, digits = 1)))
						)),
					as.expression(substitute(
						pvalue[compare]*base %*% 10^exponent, 
						list(pvalue = "P-value", compare = "quiescent", base = paste(':', t2vs13_pvalue_sci[[1]]), exponent = t2vs13_pvalue_sci[[2]])
						))
					),
				cex = 1.1
				),
			x = 0.01,
			y = 0.01,
			corner = c(0,0)
			),
	resolution = 300
	)

rstats_1vs23 <- wilcox.test(
	plot_data[plot_data$cohort == 'RAHBT' & plot_data$RNA_3 == 1, 'PI'],
	plot_data[plot_data$cohort == 'RAHBT' & plot_data$RNA_3 %in% c(2,3), 'PI']
	)
rmedian1 <- median(plot_data[plot_data$cohort == 'RAHBT' & plot_data$RNA_3 == 1, 'PI'])
rmedian23 <- median(plot_data[plot_data$cohort == 'RAHBT' & plot_data$RNA_3 %in% c(2,3), 'PI'])
rfc_1vs23 <- rmedian1-rmedian23
r1vs23_pvalue_sci <- scientific.notation(rstats_1vs23$p.value, digits = 2, type = 'list');

rstats_2vs13 <- wilcox.test(
	plot_data[plot_data$cohort == 'RAHBT' & plot_data$RNA_3 == 2, 'PI'],
	plot_data[plot_data$cohort == 'RAHBT' & plot_data$RNA_3 %in% c(1,3), 'PI']
	)
rmedian2 <- median(plot_data[plot_data$cohort == 'RAHBT' & plot_data$RNA_3 == 2, 'PI'])
rmedian13 <- median(plot_data[plot_data$cohort == 'RAHBT' & plot_data$RNA_3 %in% c(1,3), 'PI'])
rfc_2vs13 <- rmedian2-rmedian13
r2vs13_pvalue_sci <- scientific.notation(rstats_2vs13$p.value, digits = 2, type = 'list');

create.boxplot(
	PI ~ RNA_3,
	data = plot_data[plot_data$cohort == 'RAHBT',],
	add.stripplot = TRUE,
	filename = file.path('figures', paste0(date, '_RAHBT_PI_RNA3_boxplot.png')),
	ylab.axis.padding = 2,
	main = 'RAHBT',
	main.cex = 1.8,
	strip.cex = 1.5,
	ylimits = c(0,15),
	yat = seq(0,15,5),
	xaxis.rot = 45,
	xaxis.lab = c(
		expression(bold('ER'['low'])),
		expression(bold('quiescent')),
		expression(bold('ER'['high']))
		),
	xaxis.cex = 1.5,
	ylab.cex = 1.8,
	xlab.cex = 1.8,
	yaxis.cex = 1.5,
	ylab.label = 'Proliferation Index',
	xlab.label = 'RNA Subtypes',
	height = 6,
	width = 4.1,
	key = list(
			text = list(
				lab = c(
					as.expression(substitute(
						fc[compare]*value,
						list(fc = "ES", compare = "ERlow", value = paste(':', round(rfc_1vs23, digits = 1)))
						)),
					as.expression(substitute(
						pvalue[compare]*base %*% 10^exponent, 
						list(pvalue = "P-value", compare = "ERlow", base = paste(':', r1vs23_pvalue_sci[[1]]), exponent = r1vs23_pvalue_sci[[2]])
						)),
					as.expression(substitute(
						fc[compare]*value,
						list(fc = "ES", compare = "quiescent", value = paste(':', round(rfc_2vs13, digits = 1)))
						)),
					as.expression(substitute(
						pvalue[compare]*base %*% 10^exponent, 
						list(pvalue = "P-value", compare = "quiescent", base = paste(':', r2vs13_pvalue_sci[[1]]), exponent = r2vs13_pvalue_sci[[2]])
						))
					),
				cex = 1.1
				),
			x = 0.005,
			y = 0.99,
			corner = c(0,1)
			),
	resolution = 300
	)

### CREATE EMT PLOT ###############################################################################
# read in emt 
emt <- read.delim(
	'../data/EMT.Score_BRG091520.txt',
	as.is = TRUE,
	header = FALSE,
	col.names = c('Sample_ID','EMT')
	)

# create emt plot data 
emt_plot_data <- merge(emt, rdesign[,c('Sample_ID','RNA_3')], by = 'Sample_ID')
emt_plot_data$RNA_3 <- factor(emt_plot_data$RNA_3)

# create boxplot
create.boxplot(
	EMT ~ RNA_3,
	data = emt_plot_data,
	add.stripplot = TRUE,
	filename = file.path('figures', paste0(date, '_EMT_RNA3_boxplot.png')),
	ylab.axis.padding = 2,
	#xaxis.rot = 45,
	xaxis.lab = c(
		expression(bold('ER'['low'])),
		expression(bold('quiescent')),
		expression(bold('ER'['high']))
		),
	xaxis.cex = 1.5,
	ylab.cex = 1.8,
	xlab.cex = 1.8,
	yaxis.cex = 1.5,
	ylab.label = 'EMT Signature',
	xlab.label = 'RNA Subtypes',
	resolution = 300
	)