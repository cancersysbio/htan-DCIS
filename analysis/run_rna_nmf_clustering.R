### TBCRC NMF #####################################################################################
# Run NMF to cluster TBCRC DCIS samples 

### PREAMBLE ######################################################################################
# load libraries
library(BoutrosLab.plotting.general)
library(ComplexHeatmap)
library(BiocParallel)
library(RColorBrewer)
library(gprofiler2)
library(survminer)
library(survival)
library(argparse)
library(viridis)
library(DESeq2)
library(fgsea)
library(tidyr)
library(umap)
library(NMF)
library(rms)
library(vcd)

# set date
date <- Sys.Date()
### OBTAIN COMMAND LINE ARGUMENTS #################################################################
parser <- ArgumentParser();

parser$add_argument('-r', '--rna', type = 'character', 
	default = '/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data/TBCRC_resequenced_merged_n227.raw_reads.csv',
	help = 'filename of rna abundance matrix. File should be csv file.');
parser$add_argument('-d', '--design', type = 'character',
	default = '/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data/RAHBT_TBCRC_combined_design_file_SHS_021121.csv',
	help = 'filename of design file. File should be csv file.');
parser$add_argument('-g', '--gencode', type = 'character',
	default = '/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data/gencode.v33.annotation.RData',
	help = 'Gencode RData object');
parser$add_argument('-c', '--cohort', type = 'character',
	default = 'TBCRC',
	help = 'Which cohort to analyse. Options are TBCRC or RAHBT');

args <- parser$parse_args();
### LOAD RNASEQ ###################################################################################
load_rnaseq <- function(filename, design) {
	#load original reads (more data coming next week)
	rnaseq_c_raw <- read.csv(
		filename,
		header = TRUE,
		stringsAsFactors = FALSE,
		check.names = FALSE
		)
	# remove replicated columns 
	rnaseq_c_raw <- rnaseq_c_raw[,!duplicated(colnames(rnaseq_c_raw))]
	# reformat
	rnaseq_c_m <- as.matrix(rnaseq_c_raw[,-c(1,2)])
	rownames(rnaseq_c_m) <- rnaseq_c_raw$gene_name

	# dictionary symbol/ensmebl
	rna_seq_dic <- rnaseq_c_raw[,c(1:2)]

	# sort reads
	rna_matrix <- rnaseq_c_m[,match(design$Sample_ID,colnames(rnaseq_c_m))]
	# create output 
	output <- list()
	output[['abundance']] <- rna_matrix
	output[['dictionary']] <- rna_seq_dic
	return(output)
	}

### APPLY VARIANCE STABILIZING TRANSFORMATION #####################################################
# calculates variance stabilizing transformation which normalizes count data 
# ensures rna counts are approximately homoskedastic 
apply_variance_stabilizing_transformation <- function(rna_matrix) {
	# convert rna matrix to DESeq data set
	dds <- DESeqDataSetFromMatrix(
		countData = rna_matrix, 
		data.frame(cond=rep(1,ncol(rna_matrix))),
		design=~1
		)
	# estimate library size
	dds <- estimateSizeFactors(dds)
	# extracts counts normalized by size factors
	#rna_matrix_norm <- counts(dds,normalized=T)
	# apply variance stablizing transformation to counts normalized by size factors
	vsd <-varianceStabilizingTransformation(dds)
	# extract matrix of transformed values 
	rna_matrix_ep <- assay(vsd)
	return(rna_matrix_ep)
	}

normalize_by_size_factors <- function(rna_matrix) {
	# convert rna matrix to DESeq data set
	dds <- DESeqDataSetFromMatrix(
		countData = rna_matrix, 
		data.frame(cond=rep(1,ncol(rna_matrix))),
		design=~1
		)
	# estimate library size
	dds <- estimateSizeFactors(dds)
	# extracts counts normalized by size factors
	rna_matrix_norm <- counts(dds,normalized=T)
	return(rna_matrix_norm)
	}

### REMOVE NONCODING AND RIBOSOMAL PROTEINS #######################################################
remove_noncoding_and_ribosomoal_proteins <- function(rna_matrix, gencode, ribosomal = FALSE) {
	# load gencode
	load(gencode)
	# only keep protein coding genes
	gtf_df   		<- gtf_df[which(gtf_df$type=="gene"),]
	coding_genes 	<- unique(gtf_df$gene_name[which(gtf_df$gene_type=="protein_coding")])
	# find ribosomal genes
	rp_genes <- rownames(rna_matrix)[grep("^RP",rownames(rna_matrix))]
	if (ribosomal) {
		# remove ribosomal genes from coding genes
		genes_to_keep <- setdiff(coding_genes,rp_genes)
	} else {
		# if ribosomal is FALSE, keep ribosomal genes
		genes_to_keep <- coding_genes
		}
	# only keep coding, non-ribosomal genes
	rna_matrix_filt <- rna_matrix[genes_to_keep,]
	# return matrix and list of coding genes 
	output <- list()
	output$rna_matrix <- rna_matrix_filt
	output$coding_genes <- coding_genes
	return(output)
	}

### REMOVE GENES WITH ZERO VARIANCE ###############################################################
remove_genes_with_zero_variance <- function(rna_matrix) {
	# find genes with non-zero variance
	rvar_ep <- apply(rna_matrix, 1, mad, na.rm=TRUE) 
	id_rvar_ep <- rownames(rna_matrix)[which(rvar_ep!=0)]
	# only keep genes with no zero variance
	rna_matrix_nonzero <- rna_matrix[id_rvar_ep,]
	return(rna_matrix_nonzero)
	}

### ADD NMF CLUSTER TO DESIGN FILE ################################################################
add_nmf <- function(design, k_3, k_4, k_5, k_6, k_7, k_8, k_9, k_10) {
	# add NMF clusters to design matrix and return
	design$NMF_3 <- k_3[match(design$Sample_ID,names(k_3))]
	design$NMF_4 <- k_4[match(design$Sample_ID,names(k_4))]
	design$NMF_5 <- k_5[match(design$Sample_ID,names(k_5))]
	design$NMF_6 <- k_6[match(design$Sample_ID,names(k_6))]
	design$NMF_7 <- k_7[match(design$Sample_ID,names(k_7))]
	design$NMF_8 <- k_8[match(design$Sample_ID,names(k_8))]
	design$NMF_9 <- k_9[match(design$Sample_ID,names(k_9))]
	design$NMF_10 <- k_10[match(design$Sample_ID,names(k_10))]
	return(design)
}

### CREATE HEATMAP ANNOTATIONS ####################################################################
create_heatmap_annotations <- function(design, cluster_2_palette, nmf10_col) {
	# create heatmap annotations
	heatmap_anno <- HeatmapAnnotation(
		ER=design$ER_RNA,
		Her2=design$Her2_RNA,
		PAM50=design$PAM50,
		IC11 =design$IC11_RNA,
		#RNA3=design$RNA_3,
		#RNA6=design$RNA_6,
		#NMF3=design$NMF_3,
		#NMF4=design$NMF_4,
		#NMF8=design$NMF_8,
		#CNA3 = design$CNA_3,
		#CNA6 = design$CNA_6,
		#iClust4 = design$icluster_K4,
		#iClust7 = design$icluster_K7,
		#iClust8 = design$icluster_K8,
		col=list(
			ER=c("+"=cluster_2_palette[1],"-"=cluster_2_palette[2]),
			Her2=c("+"=cluster_2_palette[1],"-"=cluster_2_palette[2]),
			PAM50=c("Basal"= "#7D26CD", "Her2"= "#8B0000",  "LumB"="#0000ff",   "LumA"="#00C5CD", "Normal"="gray"),
			IC11=c("1"="#FF5500","2"= "#00EE76","3"= "#CD3278","4ER-"="cyan","4ER+"="#00C5CD",
							"5"="#8B0000", "6"="#FFFF40", "7"="#0000CD", "8"="#FFAA00",
							"9"= "#EE82EE","10"="#7D26CD")
			# RNA3 = nmf10_col,
			# RNA6 = nmf10_col,
			# CNA3 = nmf10_col,
			# CNA6 = nmf10_col
			#iClust4 = nmf10_col,
			#iClust7 = nmf10_col,
			#iClust8 = nmf10_col
			# NMF3 = nmf10_col,
			# NMF4 = nmf10_col,
			# NMF8 = nmf10_col
			),
		annotation_legend_param = list(ncol = 2),
		border = TRUE 
		)
	return(heatmap_anno)
}

### CREATE NMF CLUSTER HEATMAP ####################################################################
create_nmf_clusters_heatmap <- function(rna_matrix, res_dcis, design, cluster_2_palette, nmf10_col, nmf) {
	# extract 30 most important features for each factor
	mostImpFeat <- extractFeatures(res_dcis$fit[[as.character(nmf)]],30)
	# create heatmap 
	hm <- Heatmap(
		#rna_matrix,
		rna_matrix[unique(unlist(mostImpFeat)),],
		#column_split = design$icluster_K4,
		column_split = design[,paste0('NMF_', nmf)],
		show_row_names = FALSE,
		show_column_names = FALSE,
		show_heatmap_legend = TRUE,
		top_annotation = create_heatmap_annotations(design, cluster_2_palette, nmf10_col),
		use_raster = TRUE,
		clustering_distance_columns = "pearson"
		)
	tiff(paste0(date, '_nmf_', nmf, '_clusters_heatmap.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
	draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left")
	dev.off()
	# return most important features 
	return(mostImpFeat)
	}

### CREATE PCA AND UMAP PROJECTIONS ###############################################################
create_umap_projections <- function(design, rna_matrix) {
	# set up upmap parameters
	umap.conf 				<- umap.defaults
	umap.conf$spread		<- 1
	umap.conf$knn_repeats	<- 10
	umap.conf$verbose		<- TRUE
	umap.conf$random_state	<- 1980
	# run umap
	ep.umap <- umap(
		t(rna_matrix),
		config = umap.conf
		)
	# add umap coordinates to design file
	design$UMAP_1 <- ep.umap$layout[match(design$Sample_ID,colnames(rna_matrix)),1]
	design$UMAP_2 <- ep.umap$layout[match(design$Sample_ID,colnames(rna_matrix)),2]
	return(design)	
	}

create_pca_projections <- function(design, rna_matrix) {
	# run pca
	pca1 		<- prcomp(t(rna_matrix), rank.=5, scale=FALSE)
	pca1.ind 	<- factoextra::get_pca_ind(pca1)$coord
	#factoextra::fviz_screeplot(pca1, addlabels = TRUE, ylim = c(0, 50))
	# add pca components to design matrix
	design$PCA_1 <- pca1.ind[match(design$Sample_ID,colnames(rna_matrix)),1] 
	design$PCA_2 <- pca1.ind[match(design$Sample_ID,colnames(rna_matrix)),2]
	design$PCA_3 <- pca1.ind[match(design$Sample_ID,colnames(rna_matrix)),3]
	return(design)
}

### EXTRACT SIG AND WRITE TO FILE #################################################################
extract_sig_and_write_to_file <- function(results, filename) {
	# only keep sig results defined as Padj < 0.05
	results_sig <- results[which(results$padj < 0.05),]
	# order by effect size 
	results_sig <- results_sig[order(abs(results_sig$log2FoldChange)),]
	# write to file 
	write.table(
		results_sig,
		file = filename,
		sep = '\t',
		quote = FALSE
		)
	}

### RUN GSEA ######################################################################################
run_gsea <- function(pathway, ranks, filename, plot = TRUE) {
	# run fgsea
	fgseaRes <- fgsea(pathway, ranks, minSize=15, maxSize = 500, nperm=1000)
	# extract top upregulated pathways
	topPathwaysUp <- fgseaRes[ES > 0 & padj < 0.05]
	if (nrow(topPathwaysUp) > 10) {
		topPathwaysUp <- topPathwaysUp[head(order(pval), n=10),]
	} 
	# extract top downregulated pathways
	topPathwaysDown <- fgseaRes[ES < 0 & padj < 0.05]
	if (nrow(topPathwaysDown) > 10) {
		topPathwaysDown <- topPathwaysDown[head(order(pval), n=10),]
	}
	# combine upregulated and downregulated
	topPathways <- c(topPathwaysUp[,pathway], rev(topPathwaysDown[,pathway]))
	# plot GSEA results
	if (plot) {
		# plots results
		tiff(filename, compression = 'lzw', res = 300, width = 12, height = 7, units = 'in')
		plotGseaTable(pathway[topPathways], ranks, fgseaRes, gseaParam=0.5)
		dev.off()
	}
	# return top pathways 
	output <- rbind(topPathwaysUp, topPathwaysDown)
	return(output)
	}

### CREATE RNA ABUNDANCE BOXPLOT ##################################################################
create_cluster_boxplot <- function(gene, rna_matrix, design, nmf, ylimits = NULL, yat = TRUE, keyx = 0.01, keyy = 0.99) {
	# reformat rna 
	plot_data <- data.frame(
		Sample_ID = colnames(rna_matrix),
		exp = unlist(rna_matrix[gene,])
		)
	# add cluster 
	plot_data <- merge(
		plot_data, 
		design[,c('Sample_ID', paste0('NMF_', nmf))], 
		by = 'Sample_ID'
		)
	colnames(plot_data) <- c('sample','exp','cluster')
	# factor clusters 
	plot_data$cluster <- factor(plot_data$cluster)
	# calculate kruskal wallis test 
	stats <- kruskal.test(exp ~ cluster, data = plot_data)
	pvalue_sci <- scientific.notation(stats$p.value, digits = 2, type = 'list');
	# create boxplot 
	create.boxplot(
		exp ~ cluster,
		data = plot_data,
		add.stripplot = TRUE,
		filename = paste0(date, '_', gene, '_C', nmf, '_boxplot.tiff'),
		ylab.label = gene,
		xlab.label = 'Cluster',
		resolution = 300,
		ylimits = ylimits,
		yat = yat,
		key = list(
			text = list(
				lab = c(
					as.expression(substitute(
						base %*% 10^exponent, 
						list(base = paste0('P-value: ', pvalue_sci[[1]]), exponent = pvalue_sci[[2]])
						))
					),
				cex = 1.5
				),
			x = keyx,
			y = keyy,
			corner = c(0,1)
			)
		)
}

### CREATE GSEA DOTMAP ############################################################################
create_gsea_dotmap <- function(pathway, gsea, clusters, remove_from_label = NULL, width = 12, height = 11) {
	# extract pathway of interest 
	gseasubset <- gsea[grep(pathway, names(gsea))]
	# add cluster index 
	for (i in 1:length(gseasubset)) {
		cluster <- gsub(paste0(pathway, '_'), '', names(gseasubset)[i])
		gseasubset[[i]]$cluster <- cluster
	}
	plot_data <- do.call(rbind, gseasubset)
	# reformat plot data 
	plot_data_es <- spread(
		plot_data[,c('pathway','ES','cluster')],
		key = cluster,
		value = ES
		)
	plot_data_fdr <- spread(
		plot_data[,c('pathway','padj','cluster')],
		key = cluster,
		value = padj
		)

	# set spot size and colour functions
	spot.size.function <- function(x) {abs(x)*3}
	spot.colour.function <- function(x) {
	    colours <- rep('white', length(x));
	    colours[x > 0] <- 'firebrick3';
	    colours[x < 0] <- 'dodgerblue3';
	    colours[x == 0] <- 'transparent';
	    return(colours);
	    }

	# remove prefix from label 
	if (is.null(remove_from_label)) {
		yaxis.lab <- plot_data_es$pathway
	} else {
		yaxis.lab <- gsub(remove_from_label, '', plot_data_es$pathway)
	}

	# create key
	key_sizes <- seq(1, -1, -0.5);
	create.dotmap(
		x = plot_data_es[,-1],
		bg.data = -log10(plot_data_fdr[,-1]),
		filename = paste0(date, '_C', clusters, '_gsea_', pathway, '_dotmap.tiff'),
		xaxis.rot = 90,
		spot.colour.function = spot.colour.function,
		spot.size.function = spot.size.function,
		yaxis.lab = yaxis.lab,
		xaxis.lab = paste('    ', colnames(plot_data_es)[-1]),
		colour.scheme = c('white','black'),
		na.spot.size = 0,
		bg.alpha = 1,
		colourkey = TRUE,
		colourkey.labels = c(
			expression(1),
			expression(10^-1),
			expression(10^-2),
			expression(10^-3)
			),
		at = c(0, seq(1, 4, 0.1)),
		colourkey.labels.at = 0:3,
		colourkey.cex = 1.3,
		key = list(
			space = 'right',
			points = list(
				cex = spot.size.function(key_sizes),
				col = 'black',
				fill = spot.colour.function(key_sizes),
				pch = 21
				),
			text = list(
				lab = c('1.0','0.5','0.0','-0.5','-1.0'),
				cex = 1.2,
				adj = 1.0,
				fontface = 'bold'
				),
			title = expression(bold("Effect Size:")),
			cex.title = 1.2,
			background = 'white',
			padding = 8
			),
		width = width,
		height = height,
		bottom.padding = 3,
		left.padding = 10,
		right.padding = 5,
		resolution = 300
		)
}

### GPROFILER #####################################################################################
# create gprofiler dotmap 
create_gprofiler_dotmap <- function(database, direction, clusters, height = 10, width = 10) {
	# read in results 
	result_files <- list.files(pattern = 'gprofiler')
	# only keep direction specified 
	result_files <- grep('.txt', result_files, value = TRUE)
	result_files <- grep(direction, result_files, value = TRUE)
	result_files <- grep(paste0('C', clusters), result_files, value = TRUE)
	# read in results 
	results <- do.call(rbind, sapply(
		result_files,
		function(x) {
			cluster <- gsub('C','', strsplit(x, '_|\\.')[[1]][6])
			tmp <- read.delim(x, as.is = TRUE)
			# only keep data base of interest
			tmp <- tmp[tmp$source == database,]
			# only keep terms with intersection > 3
			tmp <- tmp[tmp$intersection_size > 3,]
			if (nrow(tmp) > 0) {
				tmp$cluster <- cluster
				return(tmp)
			} else {
				print("No sig pathways")
			}
			},
		simplify = FALSE
		))
	if (is.null(results)) {stop("No significant pathways with interaction term > 3")}
	# find top 10 pathways per cluster
	toppath <- list() 
	for (clust in unique(results$cluster)) {
		tmp <- results[results$cluster == clust,]
		tmp <- tmp[order(tmp$p_value),]
		if (nrow(tmp) > 10) {
			tmp <- tmp[1:10,]
		}
		toppath[[clust]] <- as.character(tmp$term_name)
	}
	toppath <- unlist(toppath)

	# set colour according to direction 
	if (direction == 'downregulated') {
		colour.scheme <- c('white','dodgerblue3')
	} else if (direction == 'upregulated') {
		colour.scheme <- c('white','firebrick3')
	} else {
		stop("please specify valid direction. Options are downregulated or upregulated ...")
	}

	if (length(unique(results$cluster)) > 1) {
		# create plot data 
		plot_data <- spread(
			results[results$term_name %in% toppath,c('term_name','p_value','cluster')],
			key = cluster,
			value = p_value
			)
		#create heatmap 
		create.heatmap(
			t(-log10(plot_data[,-1])),
			clustering.method = 'none',
			filename = paste0(date, '_gprofiler_C', clusters, '_', database, '_', direction, '_heatmap.tiff'),
			colour.scheme = colour.scheme,
			fill.colour = 'white',
			height = height,
			width = width,
			yaxis.lab = plot_data$term_name,
			yat = 1:nrow(plot_data),
			yaxis.cex = 1,
			xaxis.lab = paste0('   Cluster', colnames(plot_data)[-1]),
			xaxis.rot = 90,
			colourkey.labels = c(
				expression(1),
				expression(10^-2),
				expression(10^-4),
				expression(10^-6)
				),
			at = c(0, seq(1, 6, 0.1)),
			colourkey.labels.at = c(0, 2, 4, 6),
			colourkey.cex = 2,
			left.padding = 7,
			right.padding = 5,
			resolution = 300
			)
	} else {
		plot_data <- results[results$term_name %in% toppath,c('term_name','p_value','p_value')]
		#create heatmap 
		create.heatmap(
			t(-log10(plot_data[,-1])),
			clustering.method = 'none',
			filename = paste0(date, '_gprofiler_C', clusters, '_', database, '_', direction, '_heatmap.tiff'),
			colour.scheme = colour.scheme,
			fill.colour = 'white',
			height = height,
			width = width,
			xat = 1.5,
			xaxis.lab = paste0('   Cluster', unique(results$cluster)),
			yaxis.lab = plot_data$term_name,
			yat = 1:nrow(plot_data),
			yaxis.cex = 1,
			xaxis.rot = 90,
			colourkey.labels = c(
				expression(1),
				expression(10^-2),
				expression(10^-4),
				expression(10^-6)
				),
			at = c(0, seq(1, 6, 0.1)),
			colourkey.labels.at = c(0, 2, 4, 6),
			colourkey.cex = 2,
			left.padding = 7,
			right.padding = 5,
			resolution = 300
			)
	}

}

### RUN GPROFILER #################################################################################
run_gprofiler <- function(res, nmf, clusters) {
	# extract sig results 
	res <- res[which(res$padj < 0.05),]
	# subset results down to up and down regulated genes 
	cup 	<- res[res$log2FoldChange >= 2,]
	cdown 	<- res[res$log2FoldChange <= -2,] 

	# create list to return results
	results <- list()
	if (nrow(cup) > 0) {
		# order by effect size 
		cup <- cup[order(-cup$log2FoldChange),]
		# run upregulated gprofiler
		gostres_up <- gost(
			query = rownames(cup),
			organism = 'hsapiens',
			ordered_query = TRUE,
			exclude_iea = TRUE,
			correction_method = 'fdr',
			sources = c("GO:BP", "REAC")
			)
		if (!is.null(gostres_up)) {
			# write results to file
			write.table(
				gostres_up$result[,1:13],
				file = paste0(date, '_gprofiler_C', clusters, '_upregulated_genes_C',nmf,'.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE
				)
			# add to return object 
			results[['upregulated']] <- gostres_up$result
		}
	}

	if (nrow(cdown) > 0) {
		# order by effect size
		cdown <- cdown[order(cdown$log2FoldChange),]

		# run downregulated gprofiler 
		gostres_down <- gost(
			query = rownames(cdown),
			organism = 'hsapiens',
			ordered_query = TRUE,
			exclude_iea = TRUE,
			correction_method = 'fdr',
			sources = c("GO:BP", "REAC")
			)

		if (!is.null(gostres_down)) {
			# write genes to file
			write.table(
				gostres_down$result[,1:13],
				file = paste0(date, '_gprofiler_C', clusters, '_downregulated_genes_C', nmf,'.txt'),
				sep = '\t',
				row.names = FALSE,
				quote = FALSE
				)
			# add to return object 
			results[['downregulated']] <- gostres_down$result
			}
		}

	# return results 
	return(results)
}

### UPDATE IC11 ###################################################################################
update_IC11 <- function(design) {
	# update IC11 
	design$IC11_RNA <- as.character(design$IC10_RNA)
	design$IC11_RNA[design$IC10_RNA=="4"] <- ifelse(design$ER_RNA[design$IC10_RNA=="4"]=="+","4ER+","4ER-")
	return(design)
}

### MAIN ##########################################################################################
# load design matrix
design <- read.csv(
	args$design,
	header = TRUE,
	stringsAsFactors = FALSE
	)
# only keep samples in specified cohort
if (args$cohort == 'TBCRC') {
	# only keep samples in TBCRC 
	design <- design[which(design$Clustering_TBCRC == 1),]
} else if (args$cohort == 'RAHBT') {
	# only keep samples in TBCRC 
	design <- design[which(design$NEW_Clustering_RAHBT_021121 == 1),]
} else {
	stop("Please specify valid cohort. Options are TBCRC or RAHBT")
}

# add IC11 if not present 
if (!any(grepl('IC11_RNA', colnames(design)))) {
	design <- update_IC11(design)
}


#load original reads (more data coming next week)
rna <- load_rnaseq(
	filename = args$rna,
	design = design
	)
# extract mRNA abundance matrix as well as dictionary that lists gene symbol and gene names
rna_matrix <- rna$abundance
rna_dict <- rna$dictionary

# apply variance stabilizing transformation 
rna_matrix_1 <- apply_variance_stabilizing_transformation(rna_matrix)

# remove noncoding and ribosomal proteins (if ribsomal is FALSE, do not remove ribosomal genes)
rna_matrix_out <- remove_noncoding_and_ribosomoal_proteins(rna_matrix_1, gencode = args$gencode, ribosomal = FALSE)
rna_matrix_2 <- rna_matrix_out$rna_matrix
coding_genes <- rna_matrix_out$coding_genes

# remove variances with zero variance
rna_matrix_3 <- remove_genes_with_zero_variance(rna_matrix_2)

# scale by row means for plotting 
rna_matrix_scaled <- rna_matrix_3-rowMeans(rna_matrix_3)
# create heatmap of normalized, scaled and filtered rna 
tiff(filename = paste0(date, '_rna_abundance_heatmap.tiff'), compression = 'lzw')
ComplexHeatmap::Heatmap(
	rna_matrix_scaled,
	show_row_names = FALSE,
	show_column_names = FALSE,
	use_raster = TRUE
	)
dev.off()

### RUN NMF AND CREATE DIAGNOSTIC PLOTS ###########################################################
# remove outliers
#rna_matrix_4 <- rna_matrix_3[,!colnames(rna_matrix_3) %in% c("TBCRC_76","TBCRC_164")]
# standirize data (not center, matrix has to be positive)
rna_matrix_4 <- t(scale(t(rna_matrix_3),center = FALSE))

er_positive <- design[design$ER_RNA == '+',]
er_negative <- design[design$ER_RNA == '-',]
# run NMF from 2 to 10 factors, 30 times
# takes time and memory
res_dcis <- NMF::nmf(rna_matrix_4[,er_negative$Sample_ID], rank = 2:14, nrun=30,.options="vp20")
save(res_dcis, file = paste0(date, '_nmf_results_', args$cohort, '_ER_negative.RData'))
# here we evealate the clusters quality to pick the best number of clusters
tiff(filename = paste0(date, '_nmf_plots_', args$cohort, '.tiff'), compression = 'lzw',  res = 300, width = 7, height = 7, units = 'in')
plot(res_dcis) 
dev.off()
# same here, we evaluate how robist are the clisters
consensusmap(res_dcis)
tiff(filename = paste0(date, '_nmf_consensus_3.tiff'), compression = 'lzw')
consensusmap(res_dcis$consensus$`3`)
dev.off()
tiff(filename = paste0(date, '_nmf_consensus_4.tiff'), compression = 'lzw')
consensusmap(res_dcis$consensus$`4`)
dev.off()
consensusmap(res_dcis$consensus$`5`)
tiff(filename = paste0(date, '_nmf_consensus_6.tiff'), compression = 'lzw')
consensusmap(res_dcis$consensus$`6`)
dev.off()
tiff(filename = paste0(date, '_nmf_consensus_7.tiff'), compression = 'lzw')
consensusmap(res_dcis$consensus$`7`)
dev.off()
consensusmap(res_dcis$consensus$`8`)
tiff(filename = paste0(date, '_nmf_consensus_9.tiff'), compression = 'lzw')
consensusmap(res_dcis$consensus$`9`)
dev.off()
tiff(filename = paste0(date, '_nmf_consensus_10.tiff'), compression = 'lzw')
consensusmap(res_dcis$consensus$`10`)
dev.off()

# here we get the cluster labels
k_3 <- predict(res_dcis$fit$`3`,"consensus")
k_4 <- predict(res_dcis$fit$`4`,"consensus")
k_5 <- predict(res_dcis$fit$`5`,"consensus")
k_6 <- predict(res_dcis$fit$`6`,"consensus")
k_7 <- predict(res_dcis$fit$`7`,"consensus")
k_8 <- predict(res_dcis$fit$`8`,"consensus")
k_9 <- predict(res_dcis$fit$`9`,"consensus")
k_10 <- predict(res_dcis$fit$`10`,"consensus")

#plot the silhuette plot (avoid solutions with negative values)
tiff(filename = paste0(date, '_nmf_silhouette_plots_3.tiff'), compression = 'lzw')
plot(silhouette(res_dcis$fit$`3`))
dev.off()

plot(silhouette(res_dcis$fit$`5`))
tiff(filename = paste0(date, '_nmf_silhouette_plots_6.tiff'), compression = 'lzw')
plot(silhouette(res_dcis$fit$`6`))
dev.off()
tiff(filename = paste0(date, '_nmf_silhouette_plots_8.tiff'), compression = 'lzw')
plot(silhouette(res_dcis$fit$`8`))
dev.off()

tiff(filename = paste0(date, '_nmf_silhouette_plots_10.tiff'), compression = 'lzw')
plot(silhouette(res_dcis$fit$`10`))
dev.off()

# add NMF clusters to design matrix
design <- add_nmf(design, k_3, k_4, k_5, k_6, k_7, k_8, k_9, k_10)
# remove two samples not included 
#design <- design[!is.na(design$NMF_4),]
# write to file 
write.table(
	design,
	file = paste0(date, '_', args$cohort, '_design_matrix_with_NMF_clusters.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### OUTCOME ASSOCIATIONS ##########################################################################
# set colour palettes
cluster_6_palette <- colorRampPalette(brewer.pal(6, "Set1"))(6)
cluster_7_palette <- colorRampPalette(brewer.pal(7, "Set1"))(7)
cluster_2_palette <- colorRampPalette(brewer.pal(3, "Set1"))(3)
cluster_10_palette <- c("orange", "chartreuse4", "darkorchid4", "gold", "dodgerblue",
 	"firebrick3", "yellowgreen", "darkorange1", "slateblue4" , "seagreen3")

ibc_recurrence <- design[design$Diagnostic_Group %in% c('DCIS_with_IBC_recurrence','DCIS_no_recurrence'),]
dcis_recurrence <- design[design$Diagnostic_Group %in% c('DCIS_with_DCIS_recurrence','DCIS_no_recurrence'),]

# create survival plots, apply coxph model and anova
# test NMF 6
tiff(paste0(date, '_survival_NMF_8.tiff'), compression = 'lzw', res = 300, width = 7, height = 8, units = 'in')
ggsurvplot(survfit(Surv(Mo_FU,Sample_Type!="DCIS_only")~NMF_8,design),design,risk.table = TRUE,pval = T,palette =cluster_10_palette )
dev.off()

tiff(paste0(date, '_survival_DCIS_recurrence_NMF_8.tiff'), compression = 'lzw', res = 300, width = 7, height = 8, units = 'in')
ggsurvplot(survfit(Surv(Mo_FU,Sample_Type!="DCIS_only")~NMF_8,dcis_recurrence),dcis_recurrence,risk.table = TRUE,pval = T,palette =cluster_10_palette )
dev.off()

tiff(paste0(date, '_survival_IBC_recurrence_NMF_8.tiff'), compression = 'lzw', res = 300, width = 7, height = 8, units = 'in')
ggsurvplot(survfit(Surv(Mo_FU,Sample_Type!="DCIS_only")~NMF_8,ibc_recurrence),ibc_recurrence,risk.table = TRUE,pval = T,palette =cluster_10_palette )
dev.off()

# test NMF 5
ggsurvplot(survfit(Surv(FU,Sample_Type!="DCIS_only")~NMF_4,design),design,risk.table = TRUE,pval = T,palette =cluster_6_palette )
coxph(Surv(FU,Sample_Type!="DCIS_only")~NMF_5,design)
anova(coxph(Surv(FU,Sample_Type!="DCIS_only")~NMF_5,design))
survminer::ggforest(coxph(Surv(FU,Sample_Type!="DCIS_only")~NMF_5,design))

# test NMF 6 after releveling to consider 3 as the baseline
design$NMF_6 = relevel(design$NMF_6,ref="3")
coxph(Surv(FU,Sample_Type!="DCIS_only")~NMF_6,design)
survminer::ggforest(coxph(Surv(FU,Sample_Type!="DCIS_only")~NMF_6,design))


### PLOT HEATMAPS #################################################################################
# set colour palettes
nmf4_col 		<- rainbow(4)
names(nmf4_col) <- 1:4
nmf6_col 		<- cluster_6_palette
names(nmf6_col)	<- 1:6
nmf7_col 		<- cluster_7_palette
names(nmf7_col) <- 1:7
nmf10_col 		<- cluster_10_palette
names(nmf10_col) <- 1:10

# create heatmap 1
hm1 <- Heatmap(res_dcis@consensus,
            clustering_distance_rows = "pearson",
            clustering_distance_columns = "pearson",
            clustering_method_columns = "ward.D2",
            clustering_method_rows = "ward.D2"
            )
hm1

# ### ADD CNA CLUSTERS ##############################################################################
# design_clusters <- read.delim(
# 	'/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/NMF/2021-01-19_RAHBT_design_RNA_CNA.tsv',
# 	as.is = TRUE
# 	)
# rownames(design_clusters) <- design_clusters$Sample_ID
# design_clusters <- design_clusters[colnames(rna_matrix_scaled),]


# create heatmap and extract most important features
mostImpFeat <- list()
for (i in c(3,8)) {
	mostImpFeat[[i]] <- create_nmf_clusters_heatmap(
		rna_matrix = rna_matrix_scaled[,design$Sample_ID],
		res_dcis = res_dcis,
		design = design,
		cluster_2_palette = cluster_2_palette,
		nmf10_col = nmf10_col,
		nmf = i
		)
}

# get most important gene names
mostImpGenes_6 <- rownames(rna_matrix_scaled)[unique(unlist(mostImpFeat[[6]]))]

### CREATE PCA AND UMAP PROJECTIONS ###############################################################
# add umap projections to design matrix 
design <- create_umap_projections(design, rna_matrix_3)
# add pca projections to design matrix 
design <- create_pca_projections(design, rna_matrix_3)
# create umap and pca plots coloured by NMF clusters
tiff(paste0(date, '_nmf_3_UMAP.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
ggscatter(design,x="UMAP_1",y="UMAP_2",col = "NMF_3",size=5,palette = cluster_6_palette)
dev.off()
tiff(paste0(date, '_nmf_3_PCA.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
ggscatter(design,x="PCA_1",y="PCA_2",col = "NMF_3",size=5,palette = cluster_6_palette)
dev.off()
# create umap and pca plots coloured by PAM50 clusters
tiff(paste0(date, '_PAM50_UMAP.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
ggscatter(design,x="UMAP_1",y="UMAP_2",col = "PAM50",size=5,palette = cluster_6_palette)
dev.off()
tiff(paste0(date, '_PAM50_PCA.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
ggscatter(design,x="PCA_1",y="PCA_2",col = "PAM50",size=5,palette = cluster_6_palette)
dev.off()

# create umap and pca plots coloured by PAM50 clusters
tiff(paste0(date, '_ER_status_UMAP.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
ggscatter(design,x="UMAP_1",y="UMAP_2",col = "ER_RNA",size=5,palette = cluster_6_palette)
dev.off()
tiff(paste0(date, '_ER_status_PCA.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
ggscatter(design,x="PCA_1",y="PCA_2",col = "ER_RNA",size=5,palette = cluster_6_palette)
dev.off()

# create umap and pca plots coloured by PAM50 clusters
tiff(paste0(date, '_HER2_status_UMAP.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
ggscatter(design,x="UMAP_1",y="UMAP_2",col = "Her2_RNA",size=5,palette = cluster_6_palette)
dev.off()
tiff(paste0(date, '_HER2_status_PCA.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
ggscatter(design,x="PCA_1",y="PCA_2",col = "Her2_RNA",size=5,palette = cluster_6_palette)
dev.off()

IC11_palette <- default.colours(11)
# create umap and pca plots coloured by IC11 clusters
tiff(paste0(date, '_IC11_UMAP.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
ggscatter(design,x="UMAP_1",y="UMAP_2",col = "IC11_RNA",size=5,palette = IC11_palette)
dev.off()
tiff(paste0(date, '_IC11_PCA.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
ggscatter(design,x="PCA_1",y="PCA_2",col = "IC11_RNA",size=5,palette =  IC11_palette)
dev.off()

### CREATE MOSAIC PLOTS TO ASSESS CORRELATIONS BETWEEN CLUSTERS ###################################
# check correlation with ER status
tiff(paste0(date, '_nmf_7_ER_status_mosaic.tiff'), compression = 'lzw', res = 300, width = 7, height = 7, units = 'in')
mosaic(~NMF_7+ER_status, design, shade=TRUE)
dev.off()
# check correlation with PAM50
tiff(paste0(date, '_nmf_8_PAM50_mosaic.tiff'), compression = 'lzw', res = 300, width = 8, height = 7, units = 'in')
mosaic(~NMF_8+PAM50, design, shade=TRUE, 
	spacing = spacing_increase(start = 0.75, rate = 1.5),
	rot_labels=c(45,90,45,0) 
	)
dev.off()
# check correlation with IC10
tiff(paste0(date, '_nmf_8_IC11_mosaic.tiff'), compression = 'lzw', res = 300, width = 9, height = 7, units = 'in')
mosaic(~NMF_8+IC11_RNA, design, shade=TRUE,
	spacing = spacing_increase(start = 0.75, rate = 1.5),
	rot_labels=c(45,90,45,0) 
	)
dev.off()

# check correlation with diagnostic group
tiff(paste0(date, '_nmf_8_diagnostic_group_mosaic.tiff'), compression = 'lzw', res = 300, width = 9, height = 7, units = 'in')
mosaic(~NMF_8+Diagnostic_Group, design, shade=TRUE,
	spacing = spacing_increase(start = 0.75, rate = 1.5),
	rot_labels=c(45,90,45,0) 
	)
dev.off()


### RUN DIFFERENTIAL ABUNDANCE ANALYSIS ACROSS CLUSTERS ###########################################
design <- design[!is.na(design$icluster_K5),]
design <- design[design$icluster_K5 != 1,]
#use full matrix
colData <- data.frame(condition = factor(design$NMF_3))
idna <-  which(is.na(colData$condition))
# set up parallel environment
#BiocParallel::register(MulticoreParam(10))

# set up differential abundance experiment
dds <- DESeqDataSetFromMatrix(countData = rna_matrix[coding_genes,design$Sample_ID],
                              colData = colData,
                              design= ~ 0 +condition)
# run differential abundance
dds <- DESeq(dds,betaPrior =F)#, BPPARAM = MulticoreParam(10),parallel = T,betaPrior = F)
conditions <- resultsNames(dds)

# set number of clusters total 
clusters <- length(conditions)

# extract results
res <- list()
for (i in 1:length(conditions)) {
	current <- paste0('condition', i)
	others <- conditions[conditions != current]
	# extract results
	res[[i]] <- results(dds,contrast = list(current,others),listValues=c(1,-1/length(others))) 
	# extract sig and write to file 
	#extract_sig_and_write_to_file(res[[i]], filename = paste0(date, "_DE_RNA_", args$cohort, "_C", clusters, "_protCodGen_v16_C", i+1, ".txt"))
	# extract up and down regulated genes for gprofiler analysis 
	#extract_up_and_downregulated_genes(res = res[[i]], nmf = i, clusters = clusters)
}


res.lf.1 <- lfcShrink(dds, contrast = list(c("condition1"),c("condition2","condition3","condition4")),listValues=c(1,-1/5),parallel = T,BPPARAM = MulticoreParam(10),type="ashr")
res.lf.2 <- lfcShrink(dds, contrast = list(c("condition2"),c("condition1","condition3","condition4","condition5","condition6")),listValues=c(1,-1/5),parallel = T,BPPARAM = MulticoreParam(10))
res.lf.3 <- lfcShrink(dds, contrast = list(c("condition3"),c("condition2","condition1","condition4","condition5","condition6")),listValues=c(1,-1/5),parallel = T,BPPARAM = MulticoreParam(10))
res.lf.4 <- lfcShrink(dds, contrast = list(c("condition4"),c("condition2","condition3","condition1","condition5","condition6")),listValues=c(1,-1/5),parallel = T,BPPARAM = MulticoreParam(10))
res.lf.5 <- lfcShrink(dds, contrast = list(c("condition5"),c("condition2","condition3","condition4","condition1","condition6")),listValues=c(1,-1/5),parallel = T,BPPARAM = MulticoreParam(10))


### PATHWAY ENRICHMENT ANALYSIS ###################################################################
# create list of reuslts
ranks_list <- lapply(res, function(l) {
  adjpval <- -log10(l$padj)
  adjpval <- ifelse(sign(l$log2FoldChange)==1,adjpval,-1*adjpval)
  names(adjpval)=rownames(l)
  return(adjpval[!is.na(adjpval)])
} )

# download the GMT files from msigdb site (use last version)
Msig_pw <- HTSanalyzeR2::MSigDBGeneSets(species="Hs",collection = c("H"))
all_gmts <- GSEABase::getGmt("~/pathway/msigdb.v7.2.symbols.gmt")
pathways.h <- gmtPathways("~/pathway/h.all.v7.2.symbols.gmt")
pathways.cp <- gmtPathways("~/pathway/c2.cp.v7.2.symbols.gmt")
pathways.imm <- gmtPathways("~/pathway/c7.all.v7.2.symbols.gmt")
pathways.go <- gmtPathways("~/pathway/c5.all.v7.2.symbols.gmt")
pathways.pos <- gmtPathways("~/pathway/c1.all.v7.2.symbols.gmt")

# create lists for pathway output 
gsea <- list()
for (i in 1:length(ranks_list)) {
	# run gsea for cluster considering all above pathways
	gsea[[paste0('hallmark_cluster', i)]] 		<- run_gsea(pathway = pathways.h, ranks_list[[i]], 
		filename = paste0(date, '_nmf', clusters, '_C', i, '_hallmark_gsea_table.tiff'), plot = FALSE)
	gsea[[paste0('canonical_cluster', i)]] 		<- run_gsea(pathway = pathways.cp, ranks_list[[i]], 
		filename = paste0(date, '_nmf', clusters, '_C', i, '_curated_canonical_gsea_table.tiff'), plot = FALSE)
	gsea[[paste0('immunologic_cluster', i)]] 	<- run_gsea(pathway = pathways.imm, ranks_list[[i]], 
		filename = paste0(date, '_nmf', clusters, '_C', i, '_immunologic_gsea_table.tiff'), plot = FALSE)
	gsea[[paste0('ontology_cluster', i)]] 		<- run_gsea(pathway = pathways.go, ranks_list[[i]], 
		filename = paste0(date, '_nmf', clusters, '_C', i, '_ontology_gsea_table.tiff'), plot = FALSE)
	gsea[[paste0('positional_cluster', i)]] 	<- run_gsea(pathway = pathways.pos, ranks_list[[i]], 
		filename = paste0(date, '_nmf', clusters, '_C', i, '_positional_gsea_table.tiff'), plot = FALSE)
}

# create gsea dotmap for each msigdb dataset
create_gsea_dotmap(pathway = 'hallmark', gsea = gsea, clusters = clusters, 
	remove_from_label = 'HALLMARK_')
create_gsea_dotmap(pathway = 'canonical', gsea = gsea, clusters = clusters, width = 15)
create_gsea_dotmap(pathway = 'immunologic', gsea = gsea, clusters = clusters, width = 15)
create_gsea_dotmap(pathway = 'ontology', gsea = gsea, clusters = clusters, width = 15)
create_gsea_dotmap(pathway = 'positional', gsea = gsea, clusters = clusters, width = 9, height = 12)

### GPROFILER #####################################################################################
# run gprofiler 
gprofiler_res <- list()
for (i in 1:length(conditions)) {
	print(i)
	gprofiler_res[[i]] <- run_gprofiler(res = res[[i]], nmf = i+1, clusters = clusters)
}

# create gprofiler dotmap
create_gprofiler_dotmap(database = 'GO:BP', direction = 'downregulated', clusters = clusters) 
create_gprofiler_dotmap(database = 'REAC', direction = 'downregulated', clusters = clusters)

create_gprofiler_dotmap(database = 'GO:BP', direction = 'upregulated', clusters = clusters) 
create_gprofiler_dotmap(database = 'REAC', direction = 'upregulated', clusters = clusters)

### CREATE BOXPLOTS OF GENES OF INTEREST ##########################################################
# normalized by library size factors
rna_matrix_norm <- normalize_by_size_factors(rna_matrix)
# create boxplot of ribosomal genes that are downregulated in cluster 1 
clust1_down <- res[[1]][which(res[[1]]$padj < 0.05 & res[[1]]$log2FoldChange < -2),]
ribosomal <- c('RPL7','RPS13','RPS18','RPL23A','RPL10','RPS27')
for (ribo in ribosomal) {#grep('^RP', rownames(clust1_down), value = TRUE)) {
	create_cluster_boxplot(
		gene = ribo,
		rna_matrix = rna_matrix_norm,
		design = design,
		#ylimits = c(0,11),
		#yat = seq(0,10,2),
		nmf = 6
		)
}

for (hormone in c('ESR1','PGR','ERBB2')) {
	create_cluster_boxplot(
		gene = hormone,
		rna_matrix = rna_matrix_norm,
		design = design,
		#ylimits = c(0,14),
		#yat = seq(0,14,2),
		nmf = 3
		)
}

for (immunegene in c('IL12B','CCL19','CCR7', 'CD80')) {
	create_cluster_boxplot(
		gene = immunegene,
		rna_matrix = rna_matrix_norm,
		design = design,
		#ylimits = c(0,10),
		#yat = seq(0,10,2),
		nmf = 7
		)
}

### CREATE DOTMAP OF DIFFERENTIALLY EXPRESSED GENES ###############################################
# create effect size plot data
es <- lapply(res,'[[', 2)
es <- do.call(cbind, es)
colnames(es) <- paste0('C', 1:length(res))
rownames(es) <- rownames(res[[1]])
# create pvalue plot data 
p <- lapply(res,'[[', 6)
p <- do.call(cbind, p)
colnames(p) <- paste0('C', 1:length(res))
rownames(p) <- rownames(res[[1]])

# only keep genes significant in any cluster 
keep_p <- apply(
	p,
	1,
	function(x) {
		tmp <- which(x < 0.05)
		ifelse(length(tmp) > 0, TRUE, FALSE)
		})
keep_es <- apply(
	es,
	1,
	function(x) {
		tmp <- which(abs(x) >= 2)
		ifelse(length(tmp) > 0, TRUE, FALSE)
		})
keep <- which(keep_p+keep_es == 2)

# only keep siginficant genes 
es <- es[keep,]
p <- p[keep,]
# order plot data by increasing p value
es <- es[order(p[,1],p[,2], p[,3], p[,4], p[,5], p[,6], p[,7]),]
p <- p[order(p[,1],p[,2], p[,3], p[,4], p[,5], p[,6], p[,7]),]


# set spot size and colour functions
spot.size.function <- function(x) {abs(x)*0.4}
spot.colour.function <- function(x) {
    colours <- rep('white', length(x));
    colours[x > 0] <- 'firebrick3';
    colours[x < 0] <- 'dodgerblue3';
    colours[x == 0] <- 'transparent';
    return(colours);
    }

create.dotmap(
	es,
	bg.data = -log10(p),
	filename = paste0(date, '_icluster7_DE_dotmap.tiff'),
	colour.scheme = c('white','black'),
	spot.size.function = spot.size.function,
	spot.colour.function = spot.colour.function,
	yaxis.cex = 0.8,
	na.spot.size = 0,
	bg.alpha = 1,
	colourkey = TRUE,
	colourkey.labels = c(
		expression(1),
		expression(10^-1),
		expression(10^-2),
		expression(10^-3)
		),
	at = c(0, seq(1, 4, 0.1)),
	colourkey.labels.at = 0:3,
	colourkey.cex = 1.3,
	height = 15,
	resolution = 300
	)


# compared to 6 cluster solution 
er_positive <- read.delim(
	'../ER_positive_only/2021-02-10_TBCRC_design_matrix_with_NMF_clusters.txt',
	as.is = TRUE
	)
er_positive$NMF_3_code <- paste('ERpos', er_positive$NMF_3, sep = '_')
er_negative <- read.delim(
	'2021-02-10_TBCRC_design_matrix_with_NMF_clusters.txt',
	as.is = TRUE
	)
er_negative$NMF_3_code <- paste('ERneg', er_negative$NMF_3, sep = '_')
full_clusters <- read.delim(
	'../../../2021-02-02_TBCRC_design_RNA_CNA_icluster.tsv',
	as.is = TRUE
	)

er_specific <- rbind(
	er_positive[,c('Sample_ID','NMF_3_code')],
	er_negative[,c('Sample_ID','NMF_3_code')]
	)
colnames(er_specific) <- c('Sample_ID', 'RNA_3_ER')

full_clusters <- merge(full_clusters, er_specific, by = 'Sample_ID', all.x = TRUE)


tiff(paste0(date, '_rna_nmf6_er_specific_nmf3_mosaic.tiff'), compression = 'lzw', res = 300, width = 9, height = 7, units = 'in')
mosaic(~RNA_6+RNA_3_ER, full_clusters, shade=TRUE,
	spacing = spacing_increase(start = 0.75, rate = 1.5),
	rot_labels=c(45,90,45,0) 
	)
dev.off()

ibc_recurrence <- full_clusters[full_clusters$Diagnostic_Group %in% c('DCIS_with_IBC_recurrence','DCIS_no_recurrence'),]
dcis_recurrence <- full_clusters[full_clusters$Diagnostic_Group %in% c('DCIS_with_DCIS_recurrence','DCIS_no_recurrence'),]

# create survival plots, apply coxph model and anova
# test NMF 6
tiff(paste0(date, '_survival_ER_specific_NMF_3.tiff'), compression = 'lzw', res = 300, width = 7, height = 8, units = 'in')
ggsurvplot(survfit(Surv(Mo_FU,Sample_Type!="DCIS_only")~RNA_3_ER,full_clusters),full_clusters,risk.table = TRUE,pval = T,palette =cluster_10_palette )
dev.off()

tiff(paste0(date, '_survival_DCIS_recurrence_ER_specific_NMF_3.tiff'), compression = 'lzw', res = 300, width = 7, height = 8, units = 'in')
ggsurvplot(survfit(Surv(Mo_FU,Sample_Type!="DCIS_only")~RNA_3_ER,dcis_recurrence),dcis_recurrence,risk.table = TRUE,pval = T,palette =cluster_10_palette )
dev.off()

tiff(paste0(date, '_survival_IBC_recurrence_ER_specific_NMF_3.tiff'), compression = 'lzw', res = 300, width = 7, height = 8, units = 'in')
ggsurvplot(survfit(Surv(Mo_FU,Sample_Type!="DCIS_only")~RNA_3_ER,ibc_recurrence),ibc_recurrence,risk.table = TRUE,pval = T,palette =cluster_10_palette )
dev.off()