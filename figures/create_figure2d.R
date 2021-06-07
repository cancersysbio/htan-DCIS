### create_figure2d.R #############################################################################
# Plot figure 2d - densityplot of PGA

### PREAMBLE ######################################################################################
# load libraries 
library(BoutrosLab.plotting.general)


setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

date <- Sys.Date()
### CALCULATE MEASURED DNA ########################################################################
calculate_measured_dna <- function(file) {
	# read in seg file 
	seg <- read.delim(
		file,
		check.names = FALSE,
		as.is = TRUE,
		row.names = 1
		)
	# reformat seg 
	seg$chr <- sapply(rownames(seg), function(x) strsplit(x, ':')[[1]][1])
	seg$start <- as.numeric(sapply(rownames(seg), function(x) strsplit(x, ':|-')[[1]][2]))
	seg$end <- as.numeric(sapply(rownames(seg), function(x) strsplit(x, '-')[[1]][2]))
	# calculate length segment 
	seg$length <- seg$end-seg$start
	# calculate total measured length 
	measured <- sum(seg$length)
	return(measured)
	}

### CALCULATE PGA #################################################################################
calculate_pga <- function(file, measured) {
	# read in seg file 
	seg <- read.delim(
		file,
		as.is = TRUE
		)
	print(unique(seg$SAMPLE_NAME))
	pga_seg <- seg[abs(seg$LOG2_RATIO_MEAN) > 0.1,]
	# calculate length of segments 
	pga_seg$length <- pga_seg$STOP-pga_seg$START
	# calculate pga 
	pga <- sum(pga_seg$length)/measured
	# calculate pga deletions 
	del_seg <- seg[seg$LOG2_RATIO_MEAN < -0.1,]
	del_seg$length <- del_seg$STOP-del_seg$START
	del <- sum(del_seg$length)/measured
	# calculate pga gain 
	gain_seg <- seg[seg$LOG2_RATIO_MEAN > 0.1,]
	gain_seg$length <- gain_seg$STOP-gain_seg$START
	gain <- sum(gain_seg$length)/measured
	# return pga 
	data.frame(
		sample = unique(seg$SAMPLE_NAME),
		pga = pga,
		del = del,
		gain = gain
		)
	}

### MAIN ##########################################################################################
# set cna directory
cnadir <- "/oak/stanford/groups/ccurtis2/users/azizk/data/htan-pca/QDNAseq/hg38_50kb"

# read in ibc to include 
ibc <- read.delim(
	'../data/all_IBC_samples_CNA.txt',
	header = FALSE,
	as.is = TRUE
	)
segfiles <- file.path(cnadir, ibc$V1, paste0(ibc$V1, '.recal.seg'))
cnfiles <- file.path(cnadir, ibc$V1, paste0(ibc$V1, '_50_copynumber.txt'))

# calculate pga
pga <- do.call(rbind, sapply(
	ibc$V1,
	function(id) {
		# calculate length DNA measured 
		measured <- calculate_measured_dna(grep(id, cnfiles, value = TRUE))
		# calculate pga 
		calculate_pga(grep(id, segfiles, value = TRUE), measured = measured)
		},
	simplify = FALSE
	))

# write to file 
write.table(
	pga,
	file = file.path("/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data", paste0(date, '_IBC_pga.txt')),
	sep = '\t',
	row.names = FALSE
	)

### NORMAL ########################################################################################
# read in normal samples to include 
# find normal samples 
normal <- read.csv(
	"/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data/WGS_design_NORMALS_SHS_102320.csv",
	as.is = TRUE
	)
norm_segfiles <- file.path(cnadir, normal$Old_Sample_Id, paste0(normal$Old_Sample_Id, '.recal.seg'))
norm_cnfiles <- file.path(cnadir, normal$Old_Sample_Id, paste0(normal$Old_Sample_Id, '_50_copynumber.txt'))

# calculate pga
normal_pga <- do.call(rbind, sapply(
	normal$Old_Sample_Id,
	function(id) {
		# calculate length DNA measured 
		measured <- calculate_measured_dna(grep(id, norm_cnfiles, value = TRUE))
		# calculate pga 
		calculate_pga(grep(id, norm_segfiles, value = TRUE), measured = measured)
		},
	simplify = FALSE
	))

### DCIS ##########################################################################################
# read in DCIS samples to consider
dcis <- read.delim(
	'../data/DCIS_samples_in_DCIS_IBC_comparison.txt',
	header = FALSE,
	as.is = TRUE
	)

# generate dcis seg files 
dcis_pga  <- do.call(rbind, sapply(
	dcis$V1,
	function(x) {
		segfile <- file.path(cnadir,x, paste0(x, '.recal.seg'))
		cnafile <- file.path(cnadir,x, paste0(x, '_50_segmented.txt'))
		measured <- calculate_measured_dna(cnafile)
		calculate_pga(segfile, measured = measured)
		},
	simplify = FALSE
	))
dcis_pga$sample <- gsub('.recal','',dcis_pga$sample)


# read in dcis pga 
write.table(
	dcis_pga,
	file.path("/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data" paste0(date,'_TBCRC_RAHBT_pga.txt')),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

# read in design file 
design <- read.csv(
	file.path("/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data", 'RAHBT_TBCRC_combined_design_file_SHS_011321.csv'),
	as.is = TRUE,
	header = TRUE
	)
design <- design[design$Clustering_TBCRC == 1 | design$Clustering_RAHBT == 1,]
design <- design[!is.na(design$DNA_sample),]
# find cases and controls
ibc_case <- design[design$Diagnostic_Group == 'DCIS_with_IBC_recurrence','DNA_sample']
dcis_case <- design[design$Diagnostic_Group == 'DCIS_with_DCIS_recurrence','DNA_sample']
control <- design[design$Diagnostic_Group == 'DCIS_no_recurrence','DNA_sample']

# create densityplot 
create.densityplot(
	list(
		IBC = pga$pga,
		DCIS_IBCcase = dcis_pga[dcis_pga$sample %in% ibc_case,'pga'],
		DCIS_DCIScase = dcis_pga[dcis_pga$sample %in% dcis_case,'pga'],
		DCIS_control = dcis_pga[dcis_pga$sample %in% control,'pga']
		),
	filename = 'figures/Figure3c_pga_densityplot_only_TBCRC.eps',
	col = default.colours(4),
	yaxis.tck = 0,
	xaxis.tck = 0,
	#ylimits = c(0, 2),
	#yat = seq(0, 2, 0.5),
	xlab.label = 'PGA',
	legend = list(
		inside = list(
			fun = draw.key,
			args = list(
				key = list(
					points = list(
						col = default.colours(4),
						pch = 21,
						cex = 1.5,
						fill = default.colours(4)
						),
					text = list(
						lab = c(
							paste0('IBC (n=', nrow(pga), ')'),
							paste0('DCIS with IBC recurrence (n=', sum(ibc_case %in% dcis_pga$sample), ')'),
							paste0('DCIS with DCIS recurrence (n=',sum(dcis_case %in% dcis_pga$sample), ')'),
							paste0('DCIS no recurrence (n=', sum(control %in% dcis_pga$sample), ')')
							)
						),
					cex = 1
					)
				),
			x = 0.01,
			y = 0.99,
			draw = FALSE
			)
		),
	resolution = 300
	)


