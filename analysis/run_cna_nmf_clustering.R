### TBCRC NMF #####################################################################################
# Run NMF to cluster TBCRC DCIS samples 

### PREAMBLE ######################################################################################
# load libraries
library(BoutrosLab.plotting.general)
library(ComplexHeatmap)
library(BiocParallel)
library(RColorBrewer)
library(survminer)
library(argparse)
library(bedr)
library(NMF)
library(vcd)

# set date
date <- Sys.Date()
### OBTAIN COMMAND LINE ARGUMENTS #################################################################
parser <- ArgumentParser();

parser$add_argument('-c', '--cna', type = 'character', 
	default = '~/DCIS/mergeCohorts_cn_regions_15_1_21.RData',
	help = 'filename of purity-adjusted cnas');
parser$add_argument('-d', '--design', type = 'character',
	default = '~/DCIS/RAHBT_TBCRC_combined_design_file_SHS_011321.csv',
	help = 'filename of design file. File should be csv file.');

args <- parser$parse_args();
### UPDATE IC11 ###################################################################################
update_IC11 <- function(design) {
	# update IC11 
	design$IC11_RNA <- as.character(design$IC10_RNA)
	design$IC11_RNA[design$IC10_RNA=="4"] <- ifelse(design$ER_RNA[design$IC10_RNA=="4"]=="+","4ER+","4ER-")
	return(design)
}

### ADD NMF CLUSTER TO DESIGN FILE ################################################################
add_nmf <- function(design, cna_nmf) {
	for (i in 2:15) {
		k <- predict(cna_nmf$fit[[as.character(i)]],"consensus")
		design$new <- k[match(design$DNA_sample,names(k))]
		colnames(design) <- gsub('new', paste0('NMF_', i), colnames(design))
	}
	return(design)
}

### CREATE HEATMAP ANNOTATIONS ####################################################################
create_heatmap_annotations <- function(design, cluster_2_palette, nmf10_col) {
	# create heatmap annotations
	heatmap_anno <- HeatmapAnnotation(
		ER=design$ER_RNA,
		Her2=design$Her2_RNA,
		TBCRC = design$Clustering_TBCRC,
		PAM50=design$PAM50,
		IC11 =design$IC11_RNA,
		NMF4 = design$NMF_4,
		#RNA3 = design$RNA_3,
		#CNA6 = design$CNA_6,
		# NMF3=design$NMF_3,
		# NMF6=design$NMF_6,
		#iClust4 = design$icluster_K4,
		#iClust7 = design$icluster_K7,
		col=list(
			ER=c("+"=cluster_2_palette[1],"-"=cluster_2_palette[2]),
			Her2=c("+"=cluster_2_palette[1],"-"=cluster_2_palette[2]),
			TBCRC = c("1"=cluster_2_palette[1], "0"=cluster_2_palette[2]),
			PAM50=c("Basal"= "#7D26CD", "Her2"= "#8B0000",  "LumB"="#0000ff",   "LumA"="#00C5CD", "Normal"="gray"),
			IC11=c("1"="#FF5500","2"= "#00EE76","3"= "#CD3278","4ER-"="cyan","4ER+"="#00C5CD",
				"5"="#8B0000", "6"="#FFFF40", "7"="#0000CD", "8"="#FFAA00",
				"9"= "#EE82EE","10"="#7D26CD"),
			#RNA3 = nmf10_col,
			#CNA6 = nmf10_col,
			#iClust4 = nmf10_col,
			#iClust7 = nmf10_col
			NMF4 = nmf10_col
			# NMF6 = nmf10_col
			),
		annotation_legend_param = list(ncol = 2),
		border = TRUE 
		)
	return(heatmap_anno)
}

### CREATE NMF CLUSTER HEATMAP ####################################################################
create_nmf_clusters_heatmap <- function(cn_matrix, cna_nmf, design, cluster_2_palette, nmf10_col, nmf) {
	#mostImpFeat <- extractFeatures(cna_nmf$fit[[as.character(nmf)]],30)
	# create plot data 
	#plot_data <- cn_matrix[unique(unlist(mostImpFeat)),]
	plot_data <- cn_matrix
	chr <- sapply(rownames(plot_data), function(x) gsub('chr', '', strsplit(x, '\\.')[[1]][1]))
	start <- sapply(rownames(plot_data), function(x) strsplit(x, '\\.|-')[[1]][2])
	plot_data <- plot_data[order(as.numeric(chr), as.numeric(start)),]
	# set chr colours
	chr_colours <- default.colours(24, palette.type = 'chromosomes')
	names(chr_colours) <- c(1:22,'X','Y')
	# create heatmap 
	hm <- Heatmap(
		plot_data,
		#cn_matrix,
		#cn_matrix[unique(unlist(mostImpFeat)),],
		cluster_rows = FALSE,
		column_split = design[,paste0('NMF_', nmf)],
		#column_split = design$icluster_K4,
		show_row_names = FALSE,
		row_names_gp = gpar(fontsize = 10),
		show_column_names = FALSE,
		show_heatmap_legend = TRUE,
		top_annotation = create_heatmap_annotations(design, cluster_2_palette, nmf10_col),
		right_annotation = rowAnnotation(
			CHR = chr[order(as.numeric(chr))],
			col = list(
				CHR = chr_colours
				)
			),
		use_raster = TRUE
		#clustering_distance_columns = "pearson"
		)
	tiff(paste0(date, '_nmf_', nmf, '_clusters_heatmap.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
	#tiff(paste0(date, '_icluster4_clusters_heatmap.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
	draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left")
	dev.off()
	}

### EXTRACT IMPORTANT FEATURES ####################################################################
extract_important_features <- function(cn_matrix, cna_nmf, nmf) {
	# extract 30 most important features for each factor
	mostImpFeat <- extractFeatures(cna_nmf$fit[[as.character(nmf)]],30)
	# get segment names 
	segnames <- lapply(mostImpFeat, function(x) {
		rownames(cn_matrix)[x]
		})
	return(segnames)
}

### CREATE BARPLOT OF CNA COPIES FOR INFORMATIVE SEGMENTS #########################################
create_barplot_informative_segments <- function(design, cn_regions, segment, nmf) {
	# create plot data 
	cnrf <- data.frame(
		DNA_sample = rownames(cn_regions2),
		cna = cn_regions2[,segment]
		)
	plot_data <- merge(
		cnrf,
		design[,c('DNA_sample',paste0('NMF_', nmf))],
		by = 'DNA_sample'
		)
	colnames(plot_data) <- c('sample','cna','nmf')
	plot_data$nmf <- factor(plot_data$nmf)

	# create boxplot
	create.boxplot(
		cna ~ nmf,
		data = plot_data,
		main = gsub('.', ':', segment, fixed = TRUE),
		main.cex = 2,
		add.stripplot = TRUE,
		ylab.label = 'Copies',
		xlab.label = 'Cluster',
		filename = paste0(date, '_', segment, '_cna_nmf', nmf, '_boxplot.tiff'),
		resolution = 300
		)
}


### MAIN ##########################################################################################
# load design matrix
design <- read.csv(
	args$design,
	header = TRUE,
	stringsAsFactors = FALSE
	)
# only keep files to be used for clustering
design <- design[design$Clustering_TBCRC == 1 | design$Clustering_RAHBT == 1,]
# only keep samples with DNA 
design <- design[!is.na(design$DNA_sample),]

# add IC11 if not present 
if (!any(grepl('IC11_RNA', colnames(design)))) {
	design <- update_IC11(design)
}

# rerun with only cluster 6 samples
# design <- design[design$NMF_6 == 6,]

# load CNA 
load(args$cna)
# only keep samples to be used in clustering 
cn_regions <- cn_regions[rownames(cn_regions) %in% design$DNA_sample,]

# convert log2 cna to copies 
cn_regions2 <- apply(cn_regions,2,function(x) 2*(2^x))
### RUN NMF #######################################################################################
# run NMF from 2 to 10 factors, 30 times
# takes time and memory
cna_nmf <- nmf(t(cn_regions2), 2:15, nrun=30,.opt="vp20")
save(cna_nmf, file = paste0(date, '_nmf_results_cna_cluster6_only.RData'))

# add NMF clusters to design matrix
design <- add_nmf(design, cna_nmf)
# write to file 
write.table(
	design,
	file = paste0(date, '_TBCRC_RAHBT_design_matrix_with_cna_NMF_clusters_cluster6_only.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

design <- design[!is.na(design$NMF_2),]

### PLOT DIAGNOSTIC PLOTS #########################################################################
# here we evealate the clusters quality to pick the best number of clusters
tiff(filename = paste0(date, '_nmf_plots_cna.tiff'), compression = 'lzw',  res = 300, width = 7, height = 7, units = 'in')
plot(cna_nmf) 
dev.off()

# plot consensus heatmap 
for (i in c(5,6)) {
	# plot consensus heatmap
	tiff(filename = paste0(date, '_nmf_consensus_', i, '.tiff'), compression = 'lzw')
	consensusmap(cna_nmf$consensus[[as.character(i)]])
	dev.off()
	# plot silhouette plots
	tiff(filename = paste0(date, '_nmf_silhouette_plots_', i, '.tiff'), compression = 'lzw')
	plot(silhouette(cna_nmf$fit[[as.character(i)]]))
	dev.off()
	}

### PLOT HEATMAP ##################################################################################
# set cluster colors
cluster_2_palette <- colorRampPalette(brewer.pal(3, "Set1"))(3)
cluster_10_palette <- c("orange", "chartreuse4", "darkorchid4", "gold", "dodgerblue",
 	"firebrick3", "yellowgreen", "darkorange1", "slateblue4" , "seagreen3")
nmf10_col 		<- cluster_10_palette
names(nmf10_col) <- 1:10

# create heatmap and extract most important features
mostImpFeat <- list()
for (i in c(7:9)) {
	create_nmf_clusters_heatmap(
		cn_matrix = t(cn_regions[design$DNA_sample,]),
		cna_nmf = cna_nmf,
		design = design,
		cluster_2_palette = cluster_2_palette,
		nmf10_col = nmf10_col,
		nmf = i
		)
	mostImpFeat[[as.character(i)]] <- extract_important_features(
		cn_matrix = t(cn_regions2),
		cna_nmf = cna_nmf,
		nmf = i
		)
}

mostImpCyto <- list()
for (i in 1:4) {
	# find cytobands 
	cytobands <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/NMF/cna_nmf/cytoBand.txt.gz',
		as.is = TRUE,
		header = FALSE
		)
	cytoindex <- paste0(cytobands$V1, ':', cytobands$V2, '-', cytobands$V3)

	query <- mostImpFeat[[1]][[i]]
	query <- gsub('.', ':', query, fixed = TRUE)

	querye <- data.frame(
		chr = sapply(query, function(x) strsplit(x, ':')[[1]][1]),
		start = sapply(query, function(x) strsplit(x, ':|-')[[1]][2]),
		end = sapply(query, function(x) strsplit(x, '-')[[1]][2])
		)


	query <- query[order(querye$chr, as.numeric(as.character(querye$start)), as.numeric(as.character(querye$end)))]
	cytos <- in.region(cytoindex, query, check.valid = FALSE)

	tmp <- cytobands[cytos,]
	tmp$cluster <- paste0('cluster', i)
	mostImpCyto[[i]] <- tmp
 	}
 mostImpCyto <- do.call(rbind, mostImpCyto)

write.table(
	mostImpCyto,
	file = paste0(date, '_most_important_cytobands_C5.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE,
	col.names = FALSE
	)

### CHECK IF BIASES TOWARDS ONE COHORT ############################################################
# check correlation with cohort
tiff(paste0(date, '_nmf_6_cohort_mosaic.tiff'), compression = 'lzw', res = 300, width = 9, height = 7, units = 'in')
mosaic(~NMF_6+Clustering_TBCRC, design, shade=TRUE,
	spacing = spacing_increase(start = 0.75, rate = 1.5),
	rot_labels=c(45,90,45,0) 
	)
dev.off()

tiff(paste0(date, '_icluster4_diagnostic_group_mosaic.tiff'), compression = 'lzw', res = 300, width = 9, height = 7, units = 'in')
mosaic(~icluster_K4+Diagnostic_Group, design, shade=TRUE,
	spacing = spacing_increase(start = 0.75, rate = 1.5),
	rot_labels=c(45,90,45,0) 
	)
dev.off()

### TEST ENRICHMENT WITH INVASIVE SUBTYPES ########################################################
# test enrichment with PAM50
tiff(paste0(date, '_nmf_5_PAM50_mosaic.tiff'), compression = 'lzw', res = 300, width = 9, height = 7, units = 'in')
mosaic(~NMF_5+PAM50, design, shade=TRUE,
	spacing = spacing_increase(start = 0.75, rate = 1.5),
	rot_labels=c(45,90,45,0) 
	)
dev.off()

# test enrichment with IC11
tiff(paste0(date, '_nmf_5_IC11_mosaic.tiff'), compression = 'lzw', res = 300, width = 9, height = 7, units = 'in')
mosaic(~NMF_5+IC11_RNA, design, shade=TRUE,
	spacing = spacing_increase(start = 0.75, rate = 1.5),
	rot_labels=c(45,90,45,0) 
	)
dev.off()


tiff(paste0(date, '_nmf_all_cna_vs_gistic_mosaic.tiff'), compression = 'lzw', res = 300, width = 9, height = 7, units = 'in')
mosaic(~gistic_6+all_cna_6, both, shade=TRUE,
	spacing = spacing_increase(start = 0.75, rate = 1.5),
	rot_labels=c(45,90,45,0) 
	)
dev.off()

### OUTCOME ASSOCIATIONS ##########################################################################
# extract individual recurrence
ibc_recurrence <- design[design$Diagnostic_Group %in% c('DCIS_with_IBC_recurrence','DCIS_no_recurrence'),]
dcis_recurrence <- design[design$Diagnostic_Group %in% c('DCIS_with_DCIS_recurrence','DCIS_no_recurrence'),]

# create survival plots, apply coxph model and anova
# test NMF 6
tiff(paste0(date, '_survival_NMF_6.tiff'), compression = 'lzw', res = 300, width = 7, height = 8, units = 'in')
ggsurvplot(survfit(Surv(Mo_FU,Sample_Type!="DCIS_only")~NMF_6,design),design,risk.table = TRUE,pval = T,palette =cluster_10_palette )
dev.off()

tiff(paste0(date, '_survival_DCIS_recurrence_NMF_6.tiff'), compression = 'lzw', res = 300, width = 7, height = 8, units = 'in')
ggsurvplot(survfit(Surv(Mo_FU,Sample_Type!="DCIS_only")~NMF_6,dcis_recurrence),dcis_recurrence,risk.table = TRUE,pval = T,palette =cluster_10_palette )
dev.off()

tiff(paste0(date, '_survival_IBC_recurrence_NMF_6.tiff'), compression = 'lzw', res = 300, width = 7, height = 8, units = 'in')
ggsurvplot(survfit(Surv(Mo_FU,Sample_Type!="DCIS_only")~NMF_6,ibc_recurrence),ibc_recurrence,risk.table = TRUE,pval = T,palette =cluster_10_palette )
dev.off()

### CREATE BARPLOT OF CNA COPIES FOR INFORMATIVE SEGMENTS #########################################
# plot segmnent including HER2
create_barplot_informative_segments(
	design = design,
	cn_regions = cn_regions2,
	segment = 'chr17.39650001-39750001', 
	nmf = 6
	)

create_barplot_informative_segments(
	design = design,
	cn_regions = cn_regions2,
	segment = 'chr11.69300001-69650001', 
	nmf = 6
	)
create_barplot_informative_segments(
	design = design,
	cn_regions = cn_regions2,
	segment = 'chr11.69650001-69700001', 
	nmf = 6
	)


for (i in 1:6) {
	for (j in 1:5) {
		segment <- mostImpFeat[[1]][[i]][j]
		create_barplot_informative_segments(
			design = design,
			cn_regions = cn_regions2,
			segment = segment, 
			nmf = 6
			)
	}
}