### create_figure3de.R ################################################################################
# Plot figures 3d (TBCRC) and e (RAHBT)

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(argparse)

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

date <- Sys.Date()
### OBTAIN COMMAND LINE ARGUMENTS #################################################################
parser <- ArgumentParser();

parser$add_argument('-c', '--cohort', type = 'character',
	default = 'TBCRC',
	help = 'Which cohort to analyse. Options are TBCRC or RAHBT');

args <- parser$parse_args();
### MAIN ##########################################################################################
# set input and output file names depending on specified cohort
if (args$cohort == 'TBCRC') {
	plot_file <- 'plotdata/Figure2a_rna_matrix.txt'
	design_file <- '../NMF/rna_nmf/216samples_with_ribosomal/2021_01_06_TBCRC_design_matrix_with_NMF_clusters.txt'
	out_file <- 'figures/Figure2a_TBCRC_heatmap.svg'
} else if (args$cohort == 'RAHBT') {
	plot_file <- 'plotdata/Figure2b_rna_matrix.txt'
	design_file <- '../NMF/rna_nmf/RAHBT_265samples/2021-01-04_RAHBT_design_matrix_with_NMF_clusters.txt'
	out_file <- 'figures/Figure2b_RAHBT_heatmap.svg'
} else {
	stop("Please specify a valid cohort. Options are TBCRC or RAHBT ...")
}

# read in plot data for TBCRC 
plot_data <- read.delim(
	plot_file,
	check.names = FALSE,
	as.is = TRUE
	)

# read in desgin file 
design <- read.delim(
	design_file,
	as.is = TRUE
	)
design <- design[design$Sample_ID %in% colnames(plot_data),]
design$PAM50 <- factor(design$PAM50)

# read in important genes 
genes <- read.delim(
	'plotdata/Figure2a_imp_genes.txt',
	header = FALSE,
	as.is = TRUE
	)
genes <- unlist(genes)

# find sample ids for each cluster
plots <- list()
for (i in 1:3) {
	# extract sample names
	tmp <- design[design$NMF_3 == i,]
	nmf <- tmp[order(tmp$PAM50, tmp$IC11_RNA),'Sample_ID']

	# create ER and HER2 covariate 
	covariates <- data.frame(
		ESR1 = unlist(plot_data['ESR1',nmf]),
		ERBB2 = unlist(plot_data['ERBB2',nmf])
		)
	# set y axis names 
	if (i == 1) {
		yaxis.lab <- NA
		bar.lab <- c('0.0','0.5','1.0')
		ylab.label <- paste(length(genes), 'Informative Genes')
		cls.lab <-  c('PAM50','IC10')
		bar.label <- 'Proportion'
	} else {
		yaxis.lab <- NULL
		ylab.label <- ''
		cls.lab <- rep('', 2)
		bar.lab <- rep('',3)
		bar.label <- ''
	}

	xlab.labels <- c(
		expression(bold('ER'['low'])),
		expression(bold('quiescent')),
		expression(bold('ER'['high']))
		)

	# create PAM50 and IC subtypes covariate
	subtypes <- data.frame(
		PAM50 = design[match(nmf, design$Sample_ID),'PAM50'],
		IC11 = design[match(nmf, design$Sample_ID),'IC11_RNA'] 
		)
	subtypes$PAM50_num <- as.numeric(subtypes$PAM50)+11
	subtypes$IC11_num <- as.character(subtypes$IC11)
	subtypes[subtypes$IC11_num == '4ER-','IC11_num'] <- 11
	subtypes[subtypes$IC11_num == '4ER+','IC11_num'] <- 4
	subtypes$IC11_num <- as.numeric(subtypes$IC11_num)

	# create PAM50 barplot 
	pam50_col <- c("#7D26CD",  "#8B0000", "#00C5CD", "#0000ff", "gray")
	names(pam50_col) <- c('Basal','Her2','LumA','LumB','Normal')
	pam50 <- as.data.frame(table(design[design$NMF_3 == i,'PAM50']))
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
		at = seq(0.5, 16.5, 1),
		colour.scheme = c("#FF5500","#00EE76","#CD3278","#00C5CD",
				"#8B0000", "#FFFF40", "#0000CD", "#FFAA00",
				"#EE82EE", "#7D26CD", "cyan", pam50_col),
		yat = c(1,2),
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
		colour.scheme = c('dodgerblue','white','firebrick'),
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
		plot_data[genes,nmf],
		clustering.method = 'none',
		#cluster.dimensions = 'rows',
		same.as.matrix = TRUE,
		xaxis.tck = 0,
		#plot.dendrograms = plot.dendrograms,
		xaxis.top.tck = 0,
		xlab.label = xlab.labels[i],
		ylab.label = ylab.label,
		yaxis.tck = 0,
		print.colour.key = FALSE,
		colour.scheme = c('dodgerblue','white','firebrick')
		)
	}
plots <- plots[order(names(plots))]

legend <- legend.grob(
	list(
		legend = list(
			colours = c('dodgerblue','white','firebrick'),
			labels = c(-5,5),
			at = c(0,95),
			title = 'RNA\nAbundance',
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
			colours = c("#FF5500","#00EE76","#CD3278","#00C5CD","cyan",
				"#8B0000", "#FFFF40", "#0000CD", "#FFAA00",
				"#EE82EE", "#7D26CD"),
			labels = c('1','2','3','4ER+','4ER-','5','6','7','8','9','10'),
			border = 'black',
			title = expression(bold("IC10"))
			)
		),
	label.cex = 1.5,
	title.cex = 1.8
	)

# create multipanelplot 
create.multipanelplot(
	plots,
	file = out_file,
	plot.objects.heights = c(0.2, 0.1, 0.1, 0.6),
	plot.objects.width = table(design$NMF_3)/nrow(design),
	xlab.label = 'RNA Clusters',
	layout.width = 3,
	layout.height = 4,
	resolution = 300,
	right.padding = 2,
	width = 12,
	height = 10,
	legend = list(
		right = list(
			fun = legend
			)
		)
	)