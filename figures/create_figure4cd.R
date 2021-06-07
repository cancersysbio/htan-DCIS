### create_figure4cd.R ############################################################################
# Plot figure 4 c and d

### PREAMBLE ######################################################################################
# load libraries 
library(BoutrosLab.plotting.general)
library(argparse)

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

date <- Sys.Date()
### OBTAIN COMMAND LINE ARGUMENTS #################################################################
parser <- ArgumentParser();

parser$add_argument('-c', '--cna', type = 'character',
	default = '/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data/mergeCohorts_cn_regions_15_1_21.RData',
	help = 'filename of purity adjusted cnas');
parser$add_argument('-d', '--design', type = 'character',
	default = "/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/NMF/cna_nmf/TBCRC_RAHBT/all/2021-01-16_TBCRC_RAHBT_design_matrix_with_cna_NMF_clusters.txt",
	help = 'design file');

args <- parser$parse_args();
### MAIN ##########################################################################################
# read in design file 
design <- read.delim(
	args$design,
	as.is = TRUE
	)

tdesign <- read.delim(
	'../NMF/2021-02-02_TBCRC_design_RNA_CNA_icluster.tsv',
	as.is = TRUE
	)
rdesign <- read.delim(
	'../NMF/2021-02-10_RAHBT_design_RNA_CNA_icluster.tsv',
	as.is = TRUE
	)
rna3 <- rbind(tdesign[,c('DNA_sample','RNA_3')], rdesign[,c('DNA_sample','RNA_3')])

design <- merge(design, rna3, by = 'DNA_sample')


# load cna data 
load(args$cna)
# only keep samples used in clustering 
cn_regions <- cn_regions[rownames(cn_regions) %in% design$DNA_sample,]

# reduce cn regions for plotting 
cn_subset <- list()
cn_names <- list()
for (i in 1:ncol(cn_regions)) {
	j <- length(cn_subset)
	if (i == 1) {
		cn_subset[[i]] <- cn_regions[,i]
		cn_names[[i]] <- colnames(cn_regions)[i]
	} else {
		cor <- cor(
			cn_regions[,i],
			cn_subset[[j]],
			method = 'spearman'
			)
		if (cor < 0.9) {
			cn_subset[[j+1]] <- cn_regions[,i]
			cn_names[[j+1]] <- colnames(cn_regions)[i]
		}
	}
}
cn_subset <- do.call(cbind, cn_subset)
colnames(cn_subset) <- unlist(cn_names)


# find sample ids for each cluster
plots <- list()
for (i in 1:6) {
	# extract sample names
	tmp <- design[design$NMF_6 == i,]
	tmp <- tmp[order(tmp$PAM50, tmp$IC11_RNA),]
	nmf <- tmp$DNA_sample

	# create ER and HER2 covariate 
	covariates <- data.frame(
		ESR1 = (tmp$ER_RNA == '+')*1,
		ERBB2 = (tmp$Her2_RNA == '+')*1
		)
	# set y axis names 
	if (i == 1) {
		yaxis.lab <- NA
		bar.lab <- c('0.0','0.5','1.0')
		ylab.label <- 'CNA'
		cls.lab <-  c('DCIS','PAM50','IC10')
		bar.label <- 'Proportion'
	} else {
		yaxis.lab <- NULL
		ylab.label <- ''
		cls.lab <- rep('', 3)
		bar.lab <- rep('',3)
		bar.label <- ''
	}

	# create PAM50 and IC subtypes covariate
	subtypes <- data.frame(
		PAM50 = tmp$PAM50,
		IC11 = tmp$IC11_RNA,
		RNA_3 = tmp$RNA_3
		)
	subtypes$PAM50_num <- as.numeric(subtypes$PAM50)+11
	subtypes$IC11_num <- as.character(subtypes$IC11)
	subtypes[subtypes$IC11_num == '4ER-','IC11_num'] <- 11
	subtypes[subtypes$IC11_num == '4ER+','IC11_num'] <- 4
	subtypes$IC11_num <- as.numeric(subtypes$IC11_num)
	subtypes$RNA_3 <- as.numeric(subtypes$RNA_3)+16

	# create PAM50 barplot 
	pam50_col <- c("#7D26CD",  "#8B0000", "#00C5CD", "#0000ff", "gray")
	names(pam50_col) <- c('Basal','Her2','LumA','LumB','Normal')
	pam50 <- as.data.frame(table(design[design$NMF_6 == i,'PAM50']))
	pam50$Prop <- pam50$Freq/length(nmf)
	plots[[paste0('bar', i)]] <- create.barplot(
		Prop ~ Var1,
		data = pam50,
		ylimits = c(0,1),
		yat = c(0,0.5,1),
		yaxis.tck = 0,
		xaxis.tck = 0,
		col = pam50_col[as.character(pam50$Var1)],
		yaxis.lab = bar.lab,
		ylab.label = bar.label,
		xaxis.lab = rep('', nrow(pam50)),
		xlab.label = ''
		)

	# create heatmap 
	plots[[paste0('cls',i)]] <- create.heatmap(
		subtypes[,-c(1,2)],
		clustering.method = 'none',
		at = seq(0.5, 19.5, 1),
		colour.scheme = c("#FF5500","#00EE76","#CD3278","#00C5CD",
				"#8B0000", "#FFFF40", "#0000CD", "#FFAA00",
				"#EE82EE", "#7D26CD", "cyan", pam50_col, 
				"orange", "chartreuse4", "darkorchid4"),
		yat = c(1,2,3),
		yaxis.lab = cls.lab,
		xaxis.tck = 0,
		xaxis.top.tck = 0,
		yaxis.tck = 0,
		print.colour.key = FALSE,
		resolution = 300
		)


	# create heatmap 
	plots[[paste0('cov',i)]] <- create.heatmap(
		covariates,
		clustering.method = 'none',
		colour.scheme = c('dodgerblue','firebrick'),
		yaxis.lab = yaxis.lab,
		xaxis.tck = 0,
		xaxis.top.tck = 0,
		yaxis.tck = 0,
		print.colour.key = FALSE,
		resolution = 300
		)
	# create rna heatmap 
	#xlab.label <- ifelse(i == 2, 'Clusters','')
	plots[[paste0('main',i)]] <- create.heatmap(
		t(cn_regions[nmf,]),
		clustering.method = 'none',
		#cluster.dimensions = 'rows',
		same.as.matrix = TRUE,
		xaxis.tck = 0,
		#plot.dendrograms = plot.dendrograms,
		xaxis.top.tck = 0,
		xlab.label = i,
		ylab.label = ylab.label,
		yaxis.tck = 0,
		print.colour.key = FALSE,
		colour.scheme = c('dodgerblue','white','firebrick')
		)
	}
plots <- plots[order(names(plots))]

# create chr covariate 
chrs <- sapply(colnames(cn_regions), function(x) strsplit(x, '\\.')[[1]][1])
chrs <- as.numeric(gsub('chr','', chrs))
# create covariate 
plots[[paste0('side',i)]] <- create.heatmap(
	data.frame(chrs, chrs),
	clustering.method = 'none',
	same.as.matrix = TRUE,
	xaxis.tck = 0,
	xaxis.top.tck = 0,
	xlab.label = rep('', 10),
	ylab.label = rep('', ncol(cn_regions)),
	yaxis.tck = 0,
	print.colour.key = FALSE,
	colour.scheme = default.colours(22, palette.type = 'chromosomes'),
	at = seq(0.5,22.5,1)
	)

legend <- legend.grob(
	list(
		legend = list(
			colours = c('dodgerblue','white','firebrick'),
			labels = c(-2,2),
			title = expression(bold('log'[2]*' CN')),
			height = 5,
			cex = 1.5,
			continuous = TRUE
			),
		legend = list(
			colours = c("#7D26CD",  "#8B0000", "#00C5CD", "#0000ff", "gray"),
			labels = c('Basal','Her2','LumA','LumB','Normal'),
			border = 'black',
			title = expression(bold("PAM50"))
			),
		legend = list(
			colours = c("orange", "chartreuse4", "darkorchid4"),
			labels = c(
				expression('ER'['low']),
				expression('quiescent'),
				expression('ER'['high'])
				),
			border = 'black',
			title = expression(bold("DCIS Subtypes"))
			),
		legend = list(
			colours = c("#FF5500","#00EE76","#CD3278","#00C5CD","cyan",
				"#8B0000", "#FFFF40", "#0000CD", "#FFAA00",
				"#EE82EE", "#7D26CD"),
			labels = c('1','2','3','4ER+','4ER-','5','6','7','8','9','10'),
			border = 'black',
			title = expression(bold("IC10"))
			),
		legend = list(
			colours = default.colours(22, palette.type = 'chromosomes'),
			labels = as.character(1:22),
			border = 'black',
			title = expression(bold('Chromosome'))
			)
		),
	layout = c(3,2),
	label.cex = 1.5,
	title.cex = 1.5
	)

create.multipanelplot(
	plots,
	file = 'figures/Figure3a_cna_heatmap.svg',
	plot.objects.heights = c(0.2, 0.12, 0.09, 0.6),
	plot.objects.width = c(0.21, 0.12, 0.13, 0.08, 0.04, 0.48,0.05),
	layout.skip = c(rep(FALSE, 6), TRUE, rep(FALSE, 6), TRUE, rep(FALSE, 6), TRUE, rep(FALSE, 7)),
	xlab.label = 'CNA Clusters',
	layout.width = 7,
	layout.height = 4,
	resolution = 300,
	right.padding = 2,
	width = 19,
	height = 10,
	legend = list(
		right = list(
			fun = legend
			)
		)
	)


### CREATE FIGURE 3b ##############################################################################
library(tidyr)
# set informative segments
segments <- c(
	'chr17.39650001-39750001',
	'chr17.59500001-59600001',
	'chr11.69300001-69650001',
	'chr20.56350001-56400001',
	'chr8.38300001-38550001'
	)

# reformat plot data 
plot_data <- as.data.frame(cn_regions[,segments])
colnames(plot_data) <- c('17q12 (IC5)','17q23.1 (IC1)','11q13.3 (IC2)','20q13.2 (IC1/IC9)','8p11.23 (IC6)')
plot_data$DNA_sample <- rownames(plot_data)
# add clusters 
plot_data <- merge(design[,c('DNA_sample','NMF_6')], plot_data, by = 'DNA_sample')

# gather plot data
plot_data <- gather(
	plot_data[,-1],
	value = 'log2CN',
	key = 'Segment',
	-NMF_6
	)
plot_data$NMF_6 <- factor(plot_data$NMF_6)

create.boxplot(
	log2CN ~ NMF_6 | Segment,
	data = plot_data,
	filename = 'figures/Figure3b_cna_boxplot.svg',
	layout = c(5,1),
	points.cex = 1,
	strip.cex = 1.65,
	xaxis.cex = 2,
	yaxis.cex = 2,
	ylab.cex = 2,
	xlab.cex = 2,
	xlab.label = 'CNA Clusters',
	ylab.label = expression(bold('log'[2]*' CN')),
	add.stripplot = TRUE,
	width = 15,
	resolution = 300
	)



