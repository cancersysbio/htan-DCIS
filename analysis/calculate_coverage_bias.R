### calculate_coverage_bias.R #####################################################################

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# set directory
basedir <- '/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/GISTIC2'
setwd(basedir)

# set date
date <- Sys.Date()
### TEST ASSOCIATION WITH COVERAGE ################################################################
test_association_with_coverage <- function(cnas, reads) {
	# test association of each cna with coverage 
	res <- do.call(rbind, sapply(
		1:nrow(cnas),
		function(i) {
			dtf <- data.frame(
				sample = colnames(cnas),
				cna = (unlist(cnas[i,]) > 0)*1
				)
			dtf <- merge(dtf, reads, by = 'sample')
			# calculate correlation 
			stats <- wilcox.test(
				dtf[dtf$cna == 1,'mapped_reads'],
				dtf[dtf$cna == 0,'mapped_reads']
				)
			estimate <- median(dtf[dtf$cna == 0,'mapped_reads'])-median(dtf[dtf$cna == 1,'mapped_reads'])
			# check if associated with coverage 
			data.frame(
				cna = rownames(plot_data)[i],
				es = estimate,
				p = stats$p.value
				)
			},
		simplify = FALSE
		))
	res$fdr <- p.adjust(res$p, method = 'fdr')
	return(res)
	}
### TEST ASSOCIATION WITH COHORT ################################################################
test_association_with_cohort <- function(cnas, design) {
	# test association of each cna with coverage 
	res <- do.call(rbind, sapply(
		1:nrow(cnas),
		function(i) {
			dtf <- data.frame(
				DNA_sample = colnames(cnas),
				cna = (unlist(cnas[i,]) > 0)*1
				)
			dtf <- merge(dtf, design[,c('DNA_sample','Cohort')], by = 'DNA_sample')
			# calculate correlation 
			stats <- fisher.test(table(dtf[,c('Cohort','cna')]))
			# check if associated with coverage 
			data.frame(
				cna = rownames(plot_data)[i],
				es = stats$estimate[[1]],
				p = stats$p.value
				)
			},
		simplify = FALSE
		))
	res$fdr <- p.adjust(res$p, method = 'fdr')
	return(res)
	}

### GET TUMOR PURITY ##############################################################################
# get tumor purity 
get_tumor_purity <- function() {
	# set cna directory
	dir <- '/oak/stanford/groups/ccurtis2/users/azizk/data/htan-pca/QDNAseq'
	ace <- sapply(
		2:4,
		function(i) {
			read.delim(
				file.path(dir, paste0('ACE_ALL_50kb_fitpicker_', i, 'N.tsv')),
				as.is = TRUE
				)
			},
		simplify = FALSE
		)
	estimates <- do.call(cbind, lapply(ace, '[[', 2))
	purity <- data.frame(
		DNA_sample = ace[[1]][,1],
		purity = apply(estimates, 1, max)
		)
	return(purity)
}

test_association_with_purity <- function(cnas, purity) {
	# test association of each cna with coverage 
	res <- do.call(rbind, sapply(
		1:nrow(cnas),
		function(i) {
			dtf <- data.frame(
				DNA_sample = colnames(cnas),
				cna = (unlist(cnas[i,]) > 0)*1
				)
			dtf <- merge(dtf, purity, by = 'DNA_sample')
			# calculate correlation 
			stats <- wilcox.test(
				dtf[dtf$cna == 1,'purity'],
				dtf[dtf$cna == 0,'purity']
				)
			estimate <- median(dtf[dtf$cna == 0,'purity'])-median(dtf[dtf$cna == 1,'purity'])
			# check if associated with coverage 
			data.frame(
				cna = rownames(plot_data)[i],
				es = estimate,
				p = stats$p.value
				)
			},
		simplify = FALSE
		))
	res$fdr <- p.adjust(res$p, method = 'fdr')
	return(res)
	}

### MAIN ##########################################################################################
# set cna directory
cnadir <- '/oak/stanford/groups/ccurtis2/users/azizk/data/htan-pca/bam_symlink_for_AA'
# read in design 
design <- read.csv(
	'../data/RAHBT_TBCRC_combined_design_file_SHS_011321.csv',
	as.is = TRUE,
	header = TRUE
	)
design <- design[design$Clustering_TBCRC == 1 | design$Clustering_RAHBT == 1,]
design <- design[!is.na(design$DNA_sample),]

# find flagstat files
flagstat <- list.files(
	path = cnadir,
	pattern = 'flagstat',
	recursive = TRUE,
	full.names = TRUE
	)

# get number of aligned reads 
reads <- do.call(rbind, sapply(
	flagstat,
	function(filetmp) {
		sample <- strsplit(filetmp, '/')[[1]][11]
		flagstat <- read.table(filetmp, header = FALSE, as.is = TRUE, fill = TRUE)
		data.frame(
			sample = sample,
			mapped_reads = flagstat[5,1]
			)
		},
	simplify = FALSE
	))
reads$mapped_reads <- as.numeric(as.character(reads$mapped_reads))

# write to file 
write.table(
	reads, 
	file = file.path('..', paste0(date, '_DCIS_wgs_reads.txt')),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

wilcox.test(
	reads[reads$sample %in% design[design$Cohort == 'TBCRC', 'DNA_sample'],'mapped_reads'],
	reads[reads$sample %in% design[design$Cohort == 'RAHBT', 'DNA_sample'],'mapped_reads']
	)

create.densityplot(
	list(
		TBCRC = reads[reads$sample %in% design[design$Cohort == 'TBCRC', 'DNA_sample'],'mapped_reads'],
		RAHBT = reads[reads$sample %in% design[design$Cohort == 'RAHBT', 'DNA_sample'],'mapped_reads']
		),
	filename = paste0(date, '_cohort_mapped_reads_densityplot.tiff'),
	resolution = 300
	)

### CREATE DENSITYPLOT ############################################################################
# create densityplot of reads 
create.densityplot(
	list(reads = reads$mapped_reads),
	file = paste0(date, '_reads_densityplot.tiff'),
	resolution = 300
	)

### TEST ASSOCIATION WITH RECURRENT CNAS ##########################################################
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
plot_data <- peaks[which(peaks$'Amplitude Threshold' != 'Actual Copy Change Given'),]

# extract out real CNV to plot 
rownames(plot_data) <- paste(
	substr(plot_data$'Unique Name', 1, 3),
	plot_data$'Descriptor',
	sep = '_'
	)
plot_data <- plot_data[,-c(1:9, ncol(plot_data))]

# test association with cna 
res <- test_association_with_coverage(cnas = plot_data, reads = reads)

# write to file 
write.table(
	res,
	file = file.path('..', 'NMF', 'cna_nmf', paste0(date, '_recurrent_cna_associated_with_coverage.txt')),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### CHECK COVERAGE OF CNA GROUPS ##################################################################
# read in cna clusters 
clusters <- read.delim(
	file.path('..', 'NMF', 'cna_nmf','TBCRC_RAHBT','all','2021-01-16_TBCRC_RAHBT_design_matrix_with_cna_NMF_clusters.txt'),
	as.is = TRUE
	)
clusters <- clusters[,c('DNA_sample','NMF_6')]
colnames(clusters) <- c('sample','NMF_6')

# merge with reads 
clusters <- merge(clusters, reads, by = 'sample')

ks_stats <- kruskal.test(clusters$mapped_reads, clusters$NMF_6)
clusters$NMF_6 <- factor(clusters$NMF_6)
# create boxplot 
create.boxplot(
	mapped_reads ~ NMF_6,
	data = clusters,
	filename = file.path('..','NMF', 'cna_nmf','TBCRC_RAHBT', paste0(date, '_CNA6_coverage_boxplot.png')),
	add.stripplot = TRUE,
	ylab.label = 'Mapped Reads',
	xlab.label = 'CNA clusters',
	resolution = 300
	)

### CHECK COVERAGE ASSOCIATION WITH PGA ###########################################################
# read in pga 
pga <- read.delim(
	file.path('..', 'data', '2021-03-08_TBCRC_RAHBT_pga.txt'),
	as.is = TRUE
	)
pga <- merge(pga, reads, by = 'sample')

# test correlation 
cor_pga <- cor.test(pga$pga, pga$mapped_reads, method = 'spearman')

# add scatterplot 
create.scatterplot(
	pga ~ mapped_reads,
	data = pga,
	filename = file.path('..', 'manuscript','figures', paste0(date, '_pga_coverage_scatterplot.png')),
	ylab.label = 'PGA', 
	xlab.label = 'Mapped Reads',
	xat = seq(0,950e6, 200e6),
	xaxis.lab = c(
		expression(bold(0)),
		expression(bold(2 %*% 10^8)),
		expression(bold(4 %*% 10^8)),
		expression(bold(6 %*% 10^8)),
		expression(bold(8 %*% 10^8))
		),
	xaxis.cex = 1.2,
	legend = list(
		inside = list(
			fun = draw.key,
			args = list(
				key = get.corr.key(
					x = pga$pga,
					y = pga$mapped_reads,
					label.items = c('spearman'),
					alpha.background = 0,
					key.cex = 1.5
					)
				),
			x = 0.65,
			y = 0.99,
			corner = c(0,1)
			)
		),
	resolution = 300
	)

### CREATE SUPPLEMENTARY TABLE 4 ##################################################################
# test association with cna 
res_cohort <- test_association_with_cohort(cnas = plot_data, design = design)

# get tumor purity 
purity <- get_tumor_purity()
#res_purity <- test_association_with_purity(cnas = plot_data, purity = purity)

colnames(res) <- c('cna','es_coverage','p_coverage','fdr_coverage')
colnames(res_cohort) <- c('cna','es_cohort','p_cohort','fdr_cohort')
colnames(res_purity) <- c('cna','es_purity','p_purity','fdr_purity')
st4 <- merge(res[res$fdr_coverage > 0.05,], res_cohort, by = 'cna')
#st4 <- merge(st4, res_purity, by = 'cna')
