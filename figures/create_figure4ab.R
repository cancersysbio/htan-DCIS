### plot_supplementary_figures_3f.R ##############################################################
# heatmap of recurrent CNAs

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

# set date
date <- Sys.Date()
### MAIN ##########################################################################################
# set gistic dir 
gisticdir <- '/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/GISTIC2/run_brlen98_t30_armpeel_gene_res01_conf95'
# read in gistic peaks 
peaks <- read.delim(
	file.path(gisticdir, 'all_lesions.conf_95.txt'),
	as.is = TRUE,
	check.names = FALSE
	)
peaks <- peaks[which(peaks$'Residual q values after removing segments shared with higher peaks' < 0.01),]
# extract calls only
plot_data <- peaks[which(peaks$'Amplitude Threshold' == 'Actual Copy Change Given'),]

# extract out real CNV to plot 
rownames(plot_data) <- paste(
	substr(plot_data$'Unique Name', 1, 3),
	plot_data$'Descriptor',
	sep = '_'
	)
plot_data <- plot_data[,-c(1:9, ncol(plot_data))]
rownames(plot_data) <- gsub(' ', '', rownames(plot_data))

# read in cnas associated with coverage 
coverage <- read.delim(
	file.path('..', 'NMF', 'cna_nmf', '2021-03-31_recurrent_cna_associated_with_coverage.txt'),
	as.is = TRUE
	)
coverage$cna <- gsub(' ', '', coverage$cna)
coverage_sig <- coverage[coverage$fdr < 0.05,]

# create barplot plot data 
freq_data <- peaks[which(peaks$'Amplitude Threshold' != 'Actual Copy Change Given'),]
freq <- data.frame(
	Descriptor = paste(substr(freq_data[,1], 1, 3), freq_data[,2], sep = '_'),
	Count = c(
		rowSums(freq_data[,-c(1:9, ncol(freq_data))] == 1),
		rowSums(freq_data[,-c(1:9, ncol(freq_data))] == 2)
		),
	Amplitude = c(rep(1, nrow(freq_data)), rep(2, nrow(freq_data)))
	)
freq$Freq <- freq$Count/406
# remove cnas that are significantly associated with coverage 
freq$Descriptor <- gsub(' ', '', freq$Descriptor)
freq <- freq[!freq$Descriptor %in% coverage_sig$cna,]
# total count 
total <- aggregate(freq$Count, list(freq$Descriptor), sum)
colnames(total) <- c('Descriptor', 'Total')
total <- total[order(-total$Total),]
total$Index <- nrow(total):1
freq <- merge(freq, total, by = 'Descriptor')
freq <- freq[order(-freq$Index),]
# reorder plot data 
plot_data <- plot_data[unique(as.character(freq$Descriptor)),]

cov_data <- data.frame(
	grepl('Amp', rownames(plot_data))*1,
	grepl('Amp', rownames(plot_data))*1
	)

# read in tdesign 
tdesign <- read.delim(
	file.path('..', 'NMF', '2021-02-02_TBCRC_design_RNA_CNA_icluster.tsv'),
	as.is = TRUE
	)
rdesign <- read.delim(
	file.path('..', 'NMF', '2021-02-10_RAHBT_design_RNA_CNA_icluster.tsv'),
	as.is = TRUE
	)
design <- rbind(tdesign[,c('DNA_sample','RNA_3')], rdesign[,c('DNA_sample','RNA_3')])
design <- design[!is.na(design$DNA_sample),]
design <- design[order(design$RNA_3),]

plot_data <- plot_data[,design$DNA_sample]
rna_lines <- table(design$RNA_3)
rna_lines[2] <- rna_lines[2]+rna_lines[1]
rna_lines <- rna_lines[-3]+0.5

# create dataframe with RNA3
model_data <- as.data.frame(t(plot_data))
model_data$DNA_sample <- rownames(model_data)
model_data <- merge(model_data, design, by = 'DNA_sample')

# test association with RNA3 
res <- do.call(rbind, sapply(
	grep('Amp|Del', colnames(model_data), value = TRUE),
	function(i) {
		tmp <- kruskal.test(model_data[,i], model_data$RNA_3)
		fc1vs2 <- median(model_data[model_data$RNA_3 == 1,i])/median(model_data[model_data$RNA_3 == 2,i])
		fc2vs3 <- median(model_data[model_data$RNA_3 == 2,i])/median(model_data[model_data$RNA_3 == 3,i])
		fc1vs3 <- median(model_data[model_data$RNA_3 == 1,i])/median(model_data[model_data$RNA_3 == 3,i])
		##if (tmp$p.value < 0.05) {
		#	create_cna_boxplot(cna = i, dtf = model_data, pvalue = as.numeric(tmp$pvalue))
		#}
		data.frame(
			cna = i,
			fc1vs2 = fc1vs2,
			fc2vs3 = fc2vs3,
			fc1vs3 = fc1vs3,
			p = tmp$p.value
			)
		},
	simplify = FALSE
	))
res$fdr <- p.adjust(res$p, method = 'fdr')
res <- res[unique(as.character(freq$Descriptor)),]
res$index <- nrow(res):1

bar_col <- rep('black', nrow(res))
bar_col[which(res$cna == 'Amp_17q12')] <- default.colours(1)
bar_col[which(res$cna %in% c('Del_16q23.3','Del_11q25','Del_6q16.1'))] <- "darkorchid4"

# craete barplot of fdr association with RNA 
fdrplot <- create.barplot(
	index ~ -log10(fdr),
	data = res,
	col = rev(bar_col),
	plot.horizontal = TRUE,
	xat = c(2,6),
	xaxis.lab = c(
		expression(10^-2),
		expression(10^-6)
		),
	#main = "Association with\nDCIS classifier",
	#main.cex = 2,
	abline.v = 1.3,
	xaxis.cex = 2,
	xlab.cex = 1.75,
	xaxis.tck = 0,
	yaxis.tck = 0,
	ylab.label = '',
	xlab.label = 'Association with\nSubtypes (FDR)',
	yaxis.lab = rep('', nrow(res)),
	#xaxis.lab = c('0.0','0.5','1.0'),
	resolution = 300,
	)

barplot <- create.barplot(
	Index ~ Freq,
	data = freq,
	groups = freq$Amplitude,
	col = c('grey50','black'),
	stack = TRUE,
	plot.horizontal = TRUE,
	xat = c(0,0.5,1.0),
	xlimits = c(0,1),
	xaxis.tck = 0,
	yaxis.tck = 0,
	xaxis.cex = 2,
	xlab.cex = 1.75,
	ylab.label = '',
	xlab.label = 'Proportion',
	yaxis.lab = rep('', nrow(freq)),
	xaxis.lab = c('0.0','0.5','1.0'),
	resolution = 300,
	)

rna_cov <- create.heatmap(
	data.frame(design$RNA_3, design$RNA_3),
	clustering.method = 'none',
	colour.scheme = default.colours(3),
	ylab.label = '',
	xlab.label = '',
	yaxis.tck = 0,
	xaxis.tck = 0,
	print.colour.key = FALSE,
	resolution = 300
	)

heatmap <- create.heatmap(
	plot_data,
	clustering.method = 'none',
	same.as.matrix = TRUE,
	xaxis.tck = 0,
	yaxis.tck = 0,
	at = seq(-2,4,0.1),
	grid.col = TRUE,
	col.lines = rna_lines,
	col.lwd = 3,
	force.grid.col = TRUE,
	colourkey.labels.at = c(-2, 0, 4),
	colourkey.labels = c(-2, 0, 4),
	colour.scheme = c('dodgerblue','white','firebrick'),
	resolution = 300
	)

covariate <- create.heatmap(
	cov_data,
	clustering.method = 'none',
	same.as.matrix = TRUE,
	yaxis.lab = gsub('Amp_|Del_','',rownames(plot_data)),
	yat = 1:nrow(cov_data),
	yaxis.cex = 1.75,
	xaxis.tck = 0,
	yaxis.tck = 0,
	print.colour.key = FALSE,
	colour.scheme = c('dodgerblue','firebrick'),
	at = c(-0.5, 0.5, 1.5)
	) 

legend <- legend.grob(
	list(
		legend = list(
			colours = default.colours(3),
			title = "DCIS Subtypes",
			labels = c(
				expression('ER'['low']), 
				expression('quiescent'),
				expression('ER'['high'])
				),
			size = 3,
			title.cex = 1.5,
			label.cex = 1.5,
			border = 'black'
			),
		legend = list(
			colours = c('grey50','black'),
			title = "Amplitude",
			labels = c('0.3<t<0.9\n-0.3>t>-1.3', 't>0.9\nt<-1.3'),
			size = 3,
			title.cex = 1.5,
			label.cex = 1.5,
			border = 'black'
			)
		),
	layout = c(2,1),
	label.cex = 2,
	title.cex = 2,
	title.just = 'left',
	title.fontface = 'bold.italic',
	size = 2
	)


create.multipanelplot(
	list(rna_cov, covariate, heatmap, barplot, fdrplot),
	filename = 'figures/Supplementary_figure_3f_gistic_heatmap_brlen98_v2.png',
	layout.width = 4,
	layout.height = 2,
	layout.skip = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
	plot.objects.widths = c(0.15, 0.55, 0.2, 0.2),
	plot.objects.heights = c(0.08, 0.92),
	width = 16,
	height = 10,
	right.padding = 2,
	top.padding = 3,
	legend = list(
		inside = list(fun = legend, x = 0.8, y = 0.90)
		),
	resolution = 300
	)

### TEST ASSOCIATION WITH RNA #####################################################################
# read in tdesign 
tdesign <- read.delim(
	file.path('..', 'NMF', '2021-02-02_TBCRC_design_RNA_CNA_icluster.tsv'),
	as.is = TRUE
	)
rdesign <- read.delim(
	file.path('..', 'NMF', '2021-02-10_RAHBT_design_RNA_CNA_icluster.tsv'),
	as.is = TRUE
	)

# create dataframe with RNA3
model_data <- as.data.frame(t(plot_data))
model_data$DNA_sample <- rownames(model_data)
tmodel_data <- merge(model_data, tdesign[,c('DNA_sample','RNA_3')], by = 'DNA_sample')
rmodel_data <- merge(model_data, rdesign[,c('DNA_sample','RNA_3')], by = 'DNA_sample')
model_data <- rbind(tmodel_data, rmodel_data)

# test association with RNA3 
res <- do.call(rbind, sapply(
	grep('Amp|Del', colnames(model_data), value = TRUE),
	function(i) {
		tmp <- kruskal.test(model_data[,i], model_data$RNA_3)
		fc1vs2 <- median(model_data[model_data$RNA_3 == 1,i])/median(model_data[model_data$RNA_3 == 2,i])
		fc2vs3 <- median(model_data[model_data$RNA_3 == 2,i])/median(model_data[model_data$RNA_3 == 3,i])
		fc1vs3 <- median(model_data[model_data$RNA_3 == 1,i])/median(model_data[model_data$RNA_3 == 3,i])
		if (tmp$p.value < 0.05) {
			create_cna_boxplot(cna = i, dtf = model_data, pvalue = as.numeric(tmp$pvalue))
		}
		data.frame(
			cna = i,
			fc1vs2 = fc1vs2,
			fc2vs3 = fc2vs3,
			fc1vs3 = fc1vs3,
			p = tmp$p.value
			)
		},
	simplify = FALSE
	))
res$fdr <- p.adjust(res$p, method = 'fdr')

# create cna boxplot 
create_cna_boxplot <- function(cna, dtf, pvalue) {
	# reformat plot data
	plot_data <- dtf[,c(cna, 'RNA_3')] 
	colnames(plot_data) <- c('cna','rna')
	plot_data$rna <- factor(plot_data$rna)
	ylab.label <- gsub('Amp_','Amplification ', cna)
	ylab.label <- gsub('Del_','Deletion ', ylab.label)
	# create boxplot 
	create.boxplot(
		cna ~ rna,
		data = plot_data,
		filename = paste0(date, '_', cna, '_RNA3_boxplot.tiff'),
		add.stripplot = TRUE,
		resolution = 300,
		ylab.label = ylab.label,
		 legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         text = list(
                             lab = paste('P-value:', round(pvalue, digits = 2))
                             ),
                         cex = 1
                         )
                     ),
                 x = 0.01,
                 y = 0.99,
                 corner = c(0,1),
                 draw = FALSE
                 )
             ),
		xlab.label = 'RNA Clusters'
		)
}

cnas <- as.character(res[res$fdr < 0.01,'cna'])

bin_cna <- function(cna, direction) {
	if (direction == 'gain') {
		tmp <- sapply(cna, function(x) {
			ifelse(x > 0.3, )
			if (x > 0.3 & x < 0.9) {
				return(1)
			} else if (x > 0.9) {
				return(2)
			} else {
				return(0)
				}
			})
	} else if (direction == 'deletion') {
		tmp <- sapply(cna, function(x) {
			if (x < -0.3 & x > -1.3) {
				return(1)
			} else if (x < -1.3) {
				return(2)
			} else {
				return(0)
				}
			})
	}
	return(tmp)
	}

or_stats <- do.call(rbind, sapply(cnas,
	function(i) {
		dtf <- model_data[,c(i, 'RNA_3')]
		colnames(dtf) <- c('cna','rna')
		dtf$cnabin <- sapply(dtf$cna, function(x) ifelse(abs(x) > 0.3, 1, 0))
		if (i == 'Amp_17q12') {
			dtf$rna_code <- (dtf$rna == 1)*1
		} else {
			dtf$rna_code <- (dtf$rna == 3)*1
		}
		# create cont table 
		cnt <- table(dtf[,3:4])
		stats <- fisher.test(cnt)
		data.frame(
			cna = i,
			OR = stats$estimate[[1]],
			P = stats$p.value
			)
		},
	simplify = FALSE))

### TEST IF RECURRENT CNAs ASSOCIATED WITH COHORT #################################################
# test if recurrent CNAs associated with cohort
cohort_res <- do.call(rbind, sapply(
	rownames(plot_data),
	function(i) {
		tcna <- unlist(plot_data[i,colnames(plot_data) %in% tdesign$DNA_sample])
		rcna <- unlist(plot_data[i,colnames(plot_data) %in% rdesign$DNA_sample])
		stats <- wilcox.test(
			tcna,
			rcna
			)
		data.frame(
			cna = i,
			p = stats$p.value,
			es = median(tcna)-median(rcna)
			)
		},
	simplify = FALSE
	))
cohort_res$fdr <- p.adjust(cohort_res$p, method = 'fdr')

write.table(
	cohort_res,
	file = paste0(date, '_cohort_recurrent_cnas_bias.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### TEST RNA subtypes associated with PGA #########################################################
# read in pga 
pga <- read.delim(file.path('..', 'data', '2021-03-08_TBCRC_RAHBT_pga.txt'), as.is = TRUE)
colnames(pga) <- c('DNA_sample','pga')

pga <- merge(pga, design, by = 'DNA_sample')
pga$RNA_3 <- factor(pga$RNA_3)

stats <- kruskal.test(pga$pga, pga$RNA_3)

create.boxplot(
	pga ~ RNA_3,
	filename = file.path('figures', paste0(date, '_RNA3_pga_boxplot.svg')),
	data = pga,
	add.stripplot = TRUE,
	resolution = 300,
	points.cex = 1,
	xaxis.lab = c(
		expression(bold('ER'['low'])),
		expression(bold('quiescent')),
		expression(bold('ER'['high']))
		),
	ylab.label = 'PGA',
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         text = list(
                             lab = paste('P-value:', round(stats$p.value, digits = 3))
                             ),
                         cex = 1.5
                         )
                     ),
                 x = 0.01,
                 y = 0.99,
                 corner = c(0,1),
                 draw = FALSE
                 )
             ),
		xlab.label = 'RNA Clusters'
		)