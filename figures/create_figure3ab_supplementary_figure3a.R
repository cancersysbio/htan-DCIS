### call_PAM50_subtypes.R #########################################################################
# call PAM50 subtypes

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(SummarizedExperiment)
library(argparse)
library(DESeq2)
library(genefu)
library(tidyr)
library(plyr)
library(umap)

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/PAM50')


date <- Sys.Date()
### OBTAIN COMMAND LINE ARGUMENTS #################################################################
parser <- ArgumentParser();

parser$add_argument('-c', '--cohort', type = 'character',
	default = 'TBCRC',
	help = 'Which cohort to analyse. Options are TBCRC or RAHBT');

args <- parser$parse_args();
### MAIN ##########################################################################################
dir <- "/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS"
# load ssBC data 
load(file.path(dir, 'bin','ssBCsubtyping_release_20150319','ssBCsubtyping.Rdata'))
# set input and output file names depending on specified cohort
if (args$cohort == 'TBCRC') {
	rna_file <- file.path(dir, 'data','TBCRC_resequenced_merged_n227.raw_reads.csv')
	design_file <- file.path(dir, 'NMF','rna_nmf','216samples_with_ribosomal','2021_01_06_TBCRC_design_matrix_with_NMF_clusters.txt')
} else if (args$cohort == 'RAHBT') {
	rna_file <- file.path(dir, 'data','SU_WU_Read_Matrix_n933_20200505.csv')
	design_file <- file.path(dir, 'NMF','rna_nmf','RAHBT_265samples','2021-01-04_RAHBT_design_matrix_with_NMF_clusters.txt')
} else {stop("Please specify a valid cohort. Options are TBCRC or RAHBT ...")
}

# read in desgin file 
design <- read.delim(
	design_file,
	as.is = TRUE
	)

# load RNA matrix
rna <- load_rnaseq(
	filename = rna_file,
	design = design
	)
# extract mRNA abundance matrix as well as dictionary that lists gene symbol and gene names
rna_matrix <- rna$abundance
rna_dict <- rna$dictionary

# apply variance stabilizing transformation 
rna_matrix_1 <- apply_variance_stabilizing_transformation(rna_matrix)

# call PR status if not in design file
if (!any(grepl('PR_RNA', colnames(design)))) {
	pr_clust <- Mclust(data = rna_matrix_1['PGR',], modelNames = 'E',
		G = 2, na.rm = TRUE)
	design$PR_RNA <- ifelse(pr_clust$classification[match(design$Sample_ID, colnames(rna_matrix_1))] == 1, "-", "+")
}

### ssBC SUBTYPING ################################################################################
# split up samples into ER positive, ER negative and TN 
er_positive_ids <- design[design$ER_RNA == '+','Sample_ID']
er_negative_ids <- design[design$ER_RNA == '-', 'Sample_ID']
tn_ids <- design[design$ER_RNA == '-' & design$Her2_RNA == '-' & design$PR_RNA == '-','Sample_ID']
# remove triple negative samples from er negative ids 
er_negative_ids <- er_negative_ids[!er_negative_ids %in% tn_ids]

# call PAM50 per subtype
erpos_pam50 <- pam50.symbol2symbol.v2(rna_matrix_1[,er_positive_ids], s = 'ER_pos')
erneg_pam50 <- pam50.symbol2symbol.v2(rna_matrix_1[,er_negative_ids], s = 'ER_neg')
ertn_pam50 <- pam50.symbol2symbol.v2(rna_matrix_1[,tn_ids], s = 'TN')

# add to design file
pam50 <- rbind(erpos_pam50[,c('sample','call')], erneg_pam50[,c('sample','call')], ertn_pam50[,c('sample','call')])
colnames(pam50) <- c('Sample_ID','ssBC_PAM50')
design <- merge(design, pam50, by = 'Sample_ID')

### WEIGHTED MEAN SUBTYPING #######################################################################
# calculate mean ER positive and ER negative patients separately 
erpos_mean <- apply(rna_matrix_1[,er_positive_ids], 1, mean)
erneg_mean <- apply(rna_matrix_1[,c(er_negative_ids,tn_ids)], 1, mean)
# weight means 
gene_mean <- 0.6*erpos_mean+0.4*erneg_mean
# subtract gene means 
rna_matrix_2 <- rna_matrix_1-gene_mean

# predict pam50
pam50_predict <- molecular.subtyping(
	sbt.model = 'pam50',
	data = scale(t(rna_matrix_2), center = FALSE),
	annot = data.frame(Gene.Symbol = rownames(rna_matrix_2))
	)
# reformat predictions
wm <- data.frame(
	Sample_ID = names(pam50_predict$subtype),
	wm_PAM50 = pam50_predict$subtype
	)
design <- merge(design, wm, by = 'Sample_ID')

# write to file 
write.table(
	design[,grep('NMF|PCA|UMAP',colnames(design), invert = TRUE)],
	file = paste0(date, '_', args$cohort, '_PAM50.txt'),
	quote = FALSE,
	sep = '\t',
	row.names = FALSE
	)

### TRADITIONAL SUBTYPING #########################################################################
# scale 
trad_predict <- molecular.subtyping(
	sbt.model = 'pam50',
	data = scale(t(rna_matrix_1)),
	annot = data.frame(Gene.Symbol = rownames(rna_matrix_1))
	)
# reformat predictions
trad <- data.frame(
	Sample_ID = names(trad_predict$subtype),
	trad_PAM50 = trad_predict$subtype
	)
design <- merge(design, trad, by = 'Sample_ID')

### CALCULATE CORRELATION #########################################################################
# calculate distance to each centroid 
# extract pam50 centroids
pam50_centroid <- pam50$centroids
# find intersecting genes
zpam <- intersect(
	rownames(pam50_centroid),
	rownames(rna_matrix_2)
	)
rna_matrix_3 <- rna_matrix_2[zpam,]
pam_2_centroid <- pam50_centroid[zpam,]

# calculate correlation 
pam50_cor <- apply(
	rna_matrix_3,
	2,
	function(x) {
		apply(pam_2_centroid, 2, 
			function(y) cor(x, y, method='spearman'))
		})

# set pam50 covariates
sample_pam50 <- as.character(design[match(colnames(pam50_cor), design$Sample_ID),'wm_PAM50'])
sample_pam50[sample_pam50 == 'Basal'] <- "#7D26CD"
sample_pam50[sample_pam50 == 'Her2'] <- "#8B0000"
sample_pam50[sample_pam50 == 'LumA'] <- "#00C5CD"
sample_pam50[sample_pam50 == 'LumB'] <- "#0000ff"
sample_pam50[sample_pam50 == 'Normal'] <- "gray"
pam50_cov <- list(
	rect = list(
		col = 'transparent',
		fill = sample_pam50,
		lwd = 1.5
		)
	)

# plot number of each pam50 
dcis_cor <- do.call(rbind, sapply(
	c('Basal','Her2','LumA','LumB','Normal'),
	function(x) {
		data.frame(
			pam50 = x,
			cor = pam50_cor[x, design[design$wm_PAM50 == x,'Sample_ID']]
			)
		},
	simplify = FALSE
	))
boxplot <- create.boxplot(
	cor ~ pam50,
	data = dcis_cor,
	add.stripplot = TRUE,
	ylab.label = 'Correlation',
	xaxis.lab = rep('', nrow(median_cor)),
	xlab.label = '',
	yaxis.tck = 0,
	xaxis.tck = 0,
	ylimits = c(0, 1),
	yat = seq(0, 1, 0.5),
	resolution = 300
	)

# create heatmap of correlatios 
heatmap <- create.heatmap(
	pam50_cor,
	cluster.dimensions = 'rows',
	xaxis.lab = paste('    ', rownames(pam50_cor)),
	covariates = pam50_cov,
	yaxis.tck = 0,
	xaxis.tck = 0,
	ylab.label = paste(ncol(pam50_cor), 'Samples'),
	colour.scheme = c('dodgerblue','white','firebrick'),
	resolution = 300
	)

create.multipanelplot(
	list(boxplot, heatmap),
	filename = paste0(date, '_', args$cohort, '_PAM50_correlations_heatmap.tiff'),
	plot.objects.heights = c(0.25,0.75),
	height = 10,
	width = 8,
	resolution = 300
	)

### CALCULATE SILHOUETTE VALUE ####################################################################
dist_matrix <- cor(rna_matrix_3, method = 'spearman')

### CREATE PAM50 PROBABILITIES BARPLOT ############################################################
# create barplot of probabilities 
trad_prob <- as.data.frame(trad_predict$subtype.proba)
trad_prob$Sample_ID <- rownames(trad_prob)
plot_data <- gather(
	trad_prob,
	key = 'PAM50',
	value = 'probability',
	-Sample_ID
	)
plot_data <- plot_data[order(-plot_data$probability),]
plot_data <- plot_data[!duplicated(plot_data$Sample_ID),]
# order by subtype then probability 
plot_data <- plot_data[order(plot_data$PAM50, -plot_data$probability),]
plot_data$Index <- 1:nrow(plot_data)
# set col 
col <- rep(NA, nrow(plot_data))
col[plot_data$PAM50 == 'Basal'] <- "#7D26CD"
col[plot_data$PAM50 == 'Her2'] <- "#8B0000"
col[plot_data$PAM50 == 'LumA'] <- "#00C5CD"
col[plot_data$PAM50 == 'LumB'] <- "#0000ff"
col[plot_data$PAM50 == 'Normal'] <- "gray"
# create plot date 
create.barplot(
	probability ~ Index,
	data = plot_data,
	col = col,
	border.col = col,
	ylimits = c(0, 1),
	yat = seq(0, 1, 0.2),
	yaxis.tck = 0,
	xaxis.tck = 0,
	ylab.label = 'Probability',
	abline.h = 0.5,
	xaxis.lab = rep('', nrow(plot_data)),
	xlab.label = paste(nrow(plot_data), 'Samples'),
	file = paste0(date, '_', args$cohort, '_PAM50_probability_barplot.tiff'),
	resolution = 300
	)

### IBC ###########################################################################################
# load counts 
load('/N/projects/curtis/TCGA/TCGA-BRCA/GeneExp_gene_counts.RData')
load('/N/projects/curtis/TCGA/TCGA-BRCA/genefuAnnotations.RData')
load('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data/gencode.v33.annotation.RData')
# extract counts matrix
rna_ibc <- assays(data.exp)[[1]]

# apply VST 
rna_ibc_2 <- apply_variance_stabilizing_transformation(rna_ibc)

# find tissue 
tissue <- sapply(colnames(rna_ibc_2), function(x) strsplit(x, '-')[[1]][4])
rna_ibc_3 <- rna_ibc_2[,which(tissue %in% c('01A','01B','01C','06A'))]

# map ensembl ids to gene names 
gtf <- unique(gtf_df[,c('gene_id','gene_name')])
gtf$gene_id <- sapply(gtf$gene_id, function(x) strsplit(x, '\\.')[[1]][1])
anno <- gtf[match(rownames(rna_ibc), gtf$gene_id),'gene_name']

# extract pam50 genes
pam50_ids <- gtf[gtf$gene_name %in% rownames(pam50_centroid),'gene_id']
rna_ibc_4 <- rna_ibc_3[pam50_ids,]
rownames(rna_ibc_4) <- gtf[match(pam50_ids, gtf$gene_id),'gene_name']

# predict pam50 
ibc_predict <- molecular.subtyping(
	sbt.model = 'pam50',
	data = scale(t(rna_ibc_4)),
	annot = data.frame(Gene.Symbol = rownames(rna_ibc_4))
	)
ibc <- data.frame(
	ID = names(ibc_predict$subtype),
	ibc_PAM50 = ibc_predict$subtype
	)
BRCA_annotation <- merge(BRCA_annotation, ibc, by = 'ID', all.x = TRUE)

### PLOT UMAP #####################################################################################
# remove genes with zero variance 
rna_ibc_5 <- remove_genes_with_zero_variance(rna_ibc_3)
# get UMAP 
colnames(BRCA_annotation) <- gsub('ID','Sample_ID', colnames(BRCA_annotation))
BRCA_annotation <- create_umap_projections(BRCA_annotation, t(scale(t(rna_ibc_5))))

# find intersecting genes
zpam <- intersect(
	rownames(pam50_centroid),
	rownames(rna_ibc_4)
	)
rna_ibc_6 <- scale(t(rna_ibc_4[zpam,]))
pam_2_centroid <- pam50_centroid[zpam,]

# here we are building a map using a svm, mapping the umap coordinates with the expression, I tried a linear model but doesn't work so well
pam_um1 <- kernlab::ksvm(BRCA_annotation$UMAP_1~rna_ibc_6,scaled=F,kernel="vanilladot" )
pam_um2 <- kernlab::ksvm(BRCA_annotation$UMAP_2~rna_ibc_6,scaled=F,kernel="vanilladot" )
# here we map the centroid using the trained svm
pam_um1.proj <- kernlab::predict(pam_um1,t(pam_2_centroid))
pam_um2.proj <- kernlab::predict(pam_um2,t(pam_2_centroid))

# create data frame
centroids <- data.frame(
	UMAP_1=pam_um1.proj,
	UMAP_2=pam_um2.proj,
	NMF_4_2=as.factor(c(1,2,3,4,5))
	)

# create scatter plot 
plot_data <- BRCA_annotation[BRCA_annotation$UMAP_2 < 10,]
pam50_col <- c("#7D26CD",  "#8B0000", "#00C5CD", "#0000ff", "gray")
create.scatterplot(
	UMAP_2 ~ UMAP_1,
	data = plot_data,
	filename = paste0(date, '_ibc_PAM50_UMAP.tiff'),
	groups = plot_data$ibc_PAM50,
	col = pam50_col,
	add.points = TRUE,
	yaxis.tck = 0,
	xaxis.tck = 0,
	ylab.label = 'UMAP 2',
	xlab.label = 'UMAP 1',
	points.x = centroids$UMAP_1,
	points.y = centroids$UMAP_2,
	points.col = pam50_col,
	points.pch = 1,
	points.cex = 3,
	key = list(
		text = list(
			c('Basal','Her2','LumA','LumB','Normal'),
			cex = 1,
			col = 'black'
			),
		points = list(
			pch = 19,
			col = pam50_col,
			cex = 1
			),
		x = 0.01,
		y = 0.99,
		padding.text = 2
		),
	resolution = 300
	)

### PLOT CORRELATION PLOT #########################################################################
# calculate correlation 
pam50_cor <- apply(
	t(rna_ibc_6),
	2,
	function(x) {
		apply(pam_2_centroid, 2, 
			function(y) cor(x, y, method='spearman'))
		})
# plot number of each pam50 
ibc_cor <- do.call(rbind, sapply(
	c('Basal','Her2','LumA','LumB','Normal'),
	function(x) {
		data.frame(
			pam50 = x,
			cor = pam50_cor[x, BRCA_annotation[BRCA_annotation$ibc_PAM50 == x,'Sample_ID']]
			)
		},
	simplify = FALSE
	))

# merge ibc and dcis 
dcis_cor$type <- 'DCIS'
ibc_cor$type <- 'IBC'
cor_data <- rbind(dcis_cor, ibc_cor)
cor_data$group <- paste(cor_data$pam50, cor_data$type, sep = '_')

res <- do.call(rbind, sapply(
	unique(cor_data$pam50),
	function(x) {
		tmp <- cor_data[cor_data$pam50 == x,]
		stats <- wilcox.test(
			tmp[tmp$type == 'IBC', 'cor'],
			tmp[tmp$type == 'DCIS', 'cor']
			)
		fc <- median(tmp[tmp$type == 'IBC', 'cor'])/median(tmp[tmp$type == 'DCIS', 'cor'])
		median_ibc <- median(tmp[tmp$type == 'IBC', 'cor'])
		median_dcis <- median(tmp[tmp$type == 'DCIS', 'cor'])
		median_diff <- median_ibc - median_dcis
		data.frame(
			pam50 = x,
			median_ibc = median_ibc,
			median_dcis = median_dcis,
			fc = fc,
			es = median_diff,
			p = stats$p.value
			)
		},
	simplify = FALSE
	))
res$fdr <- p.adjust(res$p, method = 'fdr')
res$bonferroni <- p.adjust(res$p, method = 'bonferroni')

cov_df <- data.frame(
	rep(c(1,2), 5),
	rep(c(1,2), 5)
	)

col <- rep('black', nrow(cor_data))
col[cor_data$pam50 == 'Basal'] <- "#7D26CD"
col[cor_data$pam50 == 'Her2'] <- "#8B0000"
col[cor_data$pam50 == 'LumA'] <- "#00C5CD"
col[cor_data$pam50 == 'LumB'] <- "#0000ff"
col[cor_data$pam50 == 'Normal'] <- "grey20"

boxplot <- create.boxplot(
	cor ~ group,
	data = cor_data,
	add.stripplot = TRUE,
	#filename = paste0(date, '_DCIS_IBC_PAM50_correlation_boxplot.png'),
	ylab.label = 'Correlation',
	#xaxis.lab = rep('', nrow(median_cor)),
	width = 9,
	xaxis.cex = 2,
	yaxis.cex = 2,
	points.col = col,
	xat = seq(1.5, 11.5, 2),
	xaxis.lab = unique(cor_data$pam50),
	height = 7,
	xlab.label = 'PAM50',
	yaxis.tck = 0,
	xaxis.tck = 0,
	ylab.cex = 2.5,
	xlab.cex = 2.5,
	add.rectangle = TRUE,
	ytop.rectangle = 1,
	ybottom.rectangle = 0,
	xleft.rectangle = c(0, 4.5, 8.5),
	xright.rectangle = c(2.5, 6.5, 11),
	col.rectangle = 'grey50',
	alpha.rectangle = 0.5,
	ylimits = c(0, 1),
	yat = seq(0, 1, 0.5),
	resolution = 300
	)
cov <- create.heatmap(
	cov_df,
	clustering.method = 'none',
	colour.scheme = default.colours(2),
	resolution = 300,
	yaxis.tck = 0,
	xaxis.tck = 0,
	print.colour.key = FALSE,
	ylab.label = '',
	xlab.label = ''
	)

legend <- legend.grob(
	list(
		legend = list(
			title = 'Cancer Type',
			colours = default.colours(2),
			labels = c('DCIS','IBC'),
			border = 'black'
			),
		legend = list(
			title = 'Bonferroni',
			colours = c('white','grey50'),
			labels = c('\u2265 0.05','<0.05'),
			border = 'black'
			)
		),
	layout = c(2,1),
	title.cex = 2,
	label.cex = 2,
	size = 3
	)

create.multipanelplot(
	list(cov, boxplot),
	filename = paste0(date, '_DCIS_IBC_PAM50_correlation_boxplot.svg'),
	plot.objects.heights = c(0.07,0.93),
	legend = list(inside = list(fun = legend, x = 0.55, y = 0.93)),
	top.padding = 4,
	#width = 12,
	#height = 10,
	resolution = 300
	)

