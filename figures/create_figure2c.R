### plot_figure3c.R ##############################################################################
# plot average CNA profile

### PREAMBLE ######################################################################################
library(GenVisR)

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript/figures')

# set date
date <- Sys.Date()
### CREATE MERGED SEG FILE ########################################################################
create_merged_seg_file <- function(samples, cna_dir = "/oak/stanford/groups/ccurtis2/users/azizk/data/htan-pca/QDNAseq/hg38_50kb") {
	# create merged seg file 
	seg <- do.call(rbind, sapply(
		samples,
		function(i) {
			tmp <- file.path(cna_dir, i, paste0(i, '.recal.seg'))
			if (file.exists(tmp)) {
				read.delim(
					tmp,
					as.is = TRUE
					)
				}
			},
		simplify = FALSE
		))
	# remove number of markers and rename columns 
	seg <- seg[,-5]
	colnames(seg) <- c('sample','chromosome','start','end','segmean')
	seg <- seg[,c('chromosome','start','end','segmean','sample')]
	# convert log2 to CN 
	seg$segmean <- 2*(2^seg$segmean)
	return(seg)
	}

### FIND SAMPLE WITH HIGHER PURITY ################################################################
find_sample_higher_purity <- function(sample1, sample2) {
	# read in purity 
	purity <- sapply(
		2:4,
		function(i) {
			file <- paste0("/oak/stanford/groups/ccurtis2/users/azizk/data/htan-pca/QDNAseq/ACE_ALL_50kb_fitpicker_", i, "N.tsv")
			read.delim(file, as.is = TRUE)
			},
		simplify = FALSE
		)
	# find purity for each sample 
	s1 <- unlist(lapply(
		purity,
		function(x) x[x$Sample_ID == sample1,'likely_fit']
		))
	s2 <- unlist(lapply(
		purity,
		function(x) x[x$Sample_ID == sample2,'likely_fit']
		))
	# find max purity 
	max_purity <- max(c(s1, s2))
	if (max_purity %in% s1) {
		return(sample1)
	} else if (max_purity %in% s2) {
		return(sample2)
	} else {
		stop("Error")
		}
	}

### MAIN ##########################################################################################
# read in design file 
design <- read.csv(
	"/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data/RAHBT_TBCRC_combined_design_file_SHS_011321.csv",
	header = TRUE,
	as.is = TRUE
	)

# create segment files for different groups 
dcis_design <- design[design$Clustering_TBCRC == 1 | design$Clustering_RAHBT == 1,]
dcis_design <- dcis_design[!is.na(dcis_design$DNA_sample),]
# create seg file
all_dcis <- create_merged_seg_file(dcis_design$DNA_sample)

# create DCIS plot
svg(filename = paste0(date, '_all_DCIS_CN_plot.svg'), width = 20, height = 5)
cnFreq(all_dcis, CN_low_cutoff = 2*(2^-0.3), CN_high_cutoff = 2*(2^0.3), 
	plot_title = paste0('DCIS (n=', length(unique(all_dcis$sample)), ')'), genome = 'hg38', 
	plotChr = paste0('chr', 1:22), facet_lab_size = 7)
dev.off()

# IBC 
# read in IBC samples
ibc <- read.delim(
	'/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data/all_IBC_samples_CNA.txt',
	as.is = TRUE
	)
# create ibc seg file
ibc_rec <- create_merged_seg_file(ibc$V1)
# create IBC plot
svg(filename = paste0(date, '_all_IBC_recurrence_CN_plot.svg'), width = 20, height = 5)
cnFreq(ibc_rec, CN_low_cutoff = 2*(2^-0.3), CN_high_cutoff = 2*(2^0.3), 
	plot_title = paste0('IBC (n=', length(unique(ibc_rec$sample)), ')'), genome = 'hg38', 
	plotChr = paste0('chr', 1:22), facet_lab_size = 7)
dev.off()

# normal 
# read in normal samples
normal <- read.csv(
	"/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data/WGS_design_NORMALS_SHS_102320.csv",
	as.is = TRUE
	)
# create normal seg file
normal_rec <- create_merged_seg_file(normal$Old_Sample_Id)
# create normal plot
svg(filename = paste0(date, '_all_normal_recurrence_CN_plot.svg'), width = 20, height = 5)
cnFreq(normal_rec, CN_low_cutoff = 2*(2^-0.3), CN_high_cutoff = 2*(2^0.3), 
	plot_title = paste0('Normal (n=', length(unique(normal_rec$sample)), ')'), genome = 'hg38', 
	plotChr = paste0('chr', 1:22), facet_lab_size = 7)
dev.off()