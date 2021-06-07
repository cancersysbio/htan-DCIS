### create_supplementary_figure3mn.R ##############################################################
# Create boxplots of mRNA abundance of metabolic enzymes in TBCRC and RAHBT

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

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

### CALCULATE TPM #################################################################################
calculate_tpm <- function(rna_matrix) {
	# convert rna matrix to DESeq data set
	dds <- DESeqDataSetFromMatrix(
		countData = rna_matrix, 
		data.frame(cond=rep(1,ncol(rna_matrix))),
		design=~1
		)
	# extract counts 
	counts <- counts(dds, normalized = F)
	# estimate library size
	dds <- estimateSizeFactors(dds)
	# load gencode
	load('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data/gencode.v33.annotation.RData')
	gtf_df <- gtf_df[gtf_df$type == 'gene',]
	# keep max per gene name 
	gtf_df_median <- aggregate(gtf_df$width, list(gtf_df$gene_name), median)
	colnames(gtf_df_median) <- c('gene_name','width')
	# only keep genes in gencode 
	counts <- counts[rownames(counts) %in% gtf_df_median$gene_name,]
	# reads per kilobase 
	rpk <- counts/(gtf_df_median[match(gtf_df_median$gene_name, rownames(counts)),'width']/1000)
	scaling <- colSums(rpk)/1e6
	tpm <- rpk/scaling
	return(tpm)
	}

### CREATE TBCRC METABOLISM GENE BOXPLOTS #########################################################
# create boxplot of each gene to combine into multiplot 
create_tbcrc_metabolism_gene_boxplot <- function(rna_matrix, design, gene, fc, pvalue, ylimit, yat) {
	# reformat plot data
	plot_data <- data.frame(
		rna = unlist(rna_matrix[gene,]),
		Sample_ID = colnames(rna_matrix)
		)
	plot_data <- merge(plot_data, design[,c('Sample_ID','NMF_3')], by = 'Sample_ID')
	plot_data$NMF_3 <- factor(plot_data$NMF_3)
	# format pvalue as scientific notation
	pvalue_sci <- scientific.notation(pvalue, digits = 2, type = 'list');
	# create ylabel 
	ylabel <- ifelse(gene != 'FH', paste(gene, 'TPM'), paste(gene, 'TPM'))
	# create boxplot 
	boxplot <- create.boxplot(
		rna ~ NMF_3,
		data = plot_data,
		add.stripplot = TRUE,
		ylimit = ylimit,
		yat = yat,
		ylab.axis.padding = 2.5,
		xaxis.lab = c(
			expression(bold('ER'['low'])),
			expression(bold('quiescent')),
			expression(bold('ER'['high']))
			),
		xaxis.cex = 1.5,
		xaxis.rot = 45,
		ylab.cex = 1.8,
		xlab.cex = 1.8,
		yaxis.cex = 1.5,
		ylab.label = ylabel,
		xlab.label = '',
		key = list(
				text = list(
					lab = c(
						as.expression(substitute(
							'log'[2]*foldchange,
							list(foldchange = paste('Fold Change:', round(fc, 2)))
							)),
						as.expression(substitute(
							base %*% 10^exponent, 
							list(base = paste0('FDR: ', pvalue_sci[[1]]), exponent = pvalue_sci[[2]])
							))
						),
					cex = 1.5
					),
				x = 0.01,
				y = 0.99,
				corner = c(0,1)
				),
		resolution = 300
		)
	return(boxplot)
}

### CREATE RAHBT METABOLIC PLOTS ##################################################################
# create boxplot of each gene to combine into multiplot 
create_rahbt_metabolism_gene_boxplot <- function(rna_matrix, design, gene, fc, pvalue, ylimit, yat) {
	# reformat plot data
	plot_data <- data.frame(
		rna = unlist(rna_matrix[gene,]),
		Sample_ID = colnames(rna_matrix)
		)
	plot_data <- merge(plot_data, design[,c('Sample_ID','NMF_3')], by = 'Sample_ID')
	plot_data$NMF_3 <- factor(plot_data$NMF_3)
	# format pvalue as scientific notation
	pvalue_sci <- scientific.notation(pvalue, digits = 2, type = 'list');
	# create ylabel 
	ylabel <- ifelse(gene != 'FH', paste(gene, 'TPM'), paste(gene, 'TPM'))
	# create boxplot 
	boxplot <- create.boxplot(
		rna ~ NMF_3,
		data = plot_data,
		add.stripplot = TRUE,
		ylimit = ylimit,
		yat = yat,
		ylab.axis.padding = 2.5,
		xaxis.lab = c(
			expression(bold('ER'['low'])),
			expression(bold('quiescent')),
			expression(bold('ER'['high'])),
			expression(bold('normal'))
			),
		xaxis.cex = 1.5,
		xaxis.rot = 45,
		ylab.cex = 1.8,
		xlab.cex = 1.8,
		yaxis.cex = 1.5,
		ylab.label = ylabel,
		xlab.label = '',
		key = list(
				text = list(
					lab = c(
						as.expression(substitute(
							'log'[2]*foldchange,
							list(foldchange = paste('Fold Change:', round(fc, 2)))
							)),
						as.expression(substitute(
							base %*% 10^exponent, 
							list(base = paste0('FDR: ', pvalue_sci[[1]]), exponent = pvalue_sci[[2]])
							))
						),
					cex = 1.5
					),
				x = 0.01,
				y = 0.99,
				corner = c(0,1)
				),
		resolution = 300
		)
	return(boxplot)
}

### TBCRC ##########################################################################################
#load original reads (more data coming next week)
rna <- load_rnaseq(
	filename = '/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data/TBCRC_resequenced_merged_n227.raw_reads.csv',
	design = design
	)
# extract mRNA abundance matrix as well as dictionary that lists gene symbol and gene names
rna_matrix <- rna$abundance
rna_dict <- rna$dictionary

# calculate tpm 
rna_matrix_1 <- calculate_tpm(rna_matrix)

# read in de results 
de <- read.delim(
	'../NMF/rna_nmf/216samples_with_ribosomal/2020-12-29_DE_RNA_TBCRC_C3_protCodGen_v16_C2.txt',
	as.is = TRUE
	)
de <- de[de$padj < 0.05,]
de <- de[sign(de$log2FoldChange) == -1,]
de_glycolysis <- de[rownames(de) %in% pathways.h[[35]],]
de_op <- de[rownames(de) %in% pathways.h[[34]],]

# create boxplots
pkm <- create_tbcrc_metabolism_gene_boxplot(
	rna_matrix = rna_matrix_1, 
	design = design, 
	gene = 'PKM', 
	fc = de['PKM','log2FoldChange'], 
	pvalue = de['PKM','padj'], 
	ylimit = c(0,5000), 
	yat = seq(0,5000,2500)
	)
pgk1 <- create_tbcrc_metabolism_gene_boxplot(
	rna_matrix = rna_matrix_1, 
	design = design, 
	gene = 'PGK1', 
	fc = de['PGK1','log2FoldChange'], 
	pvalue = de['PGK1','padj'], 
	ylimit = c(0,1000), 
	yat = seq(0,1000,500)
	)
idh1 <- create_tbcrc_metabolism_gene_boxplot(
	rna_matrix = rna_matrix_1, 
	design = design, 
	gene = 'IDH1', 
	fc = de['IDH1','log2FoldChange'], 
	pvalue = de['IDH1','padj'], 
	ylimit = c(0,1000), 
	yat = seq(0,1000,500)
	)
fh <- create_tbcrc_metabolism_gene_boxplot(
	rna_matrix = rna_matrix_1, 
	design = design, 
	gene = 'FH', 
	fc = de['FH','log2FoldChange'], 
	pvalue = de['FH','padj'], 
	ylimit = c(0,200), 
	yat = seq(0,200,100)
	)
glt8 <- create_tbcrc_metabolism_gene_boxplot(
	rna_matrix = rna_matrix_1, 
	design = design, 
	gene = 'SLC2A8', 
	fc = de['SLC2A8','log2FoldChange'], 
	pvalue =  de['SLC2A8','padj'], 
	ylimit = c(0,100), 
	yat = seq(0,100,50)
	)

create.multipanelplot(
	list(pkm, pgk1, idh1, fh, glt8),
	layout.width = 5,
	layout.height = 1,
	width = 24,
	height = 8,
	filename = 'figures/SupplementaryFigure3N_metabolism_genes_boxplot.png',
	xlab.label = 'RNA Subtypes',
	resolution = 300
	)

## RAHBT ##########################################################################################
# read in de results 
rde <- read.delim(
	'../NMF/rna_nmf/RAHBT_265samples/2021-05-03_DE_RNA_C2.txt',
	as.is = TRUE
	)

# read in design file 
rdesign <- read.delim(
	'../NMF/rna_nmf/RAHBT_265samples/2021-01-04_RAHBT_design_matrix_with_NMF_clusters.txt',
	as.is = TRUE
	)
normal <- read.csv(
	'../data/RAHBT_TBCRC_combined_design_file_SHS_011321.csv',
	as.is = TRUE,
	header = TRUE
	)
normal <- normal[normal$Sample_Type == 'Normal',]

#load original reads (more data coming next week)
rrna <- load_rnaseq(
	filename = '/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data/SU_WU_Read_Matrix_n933_20200505.csv',
	design = rbind(rdesign[,colnames(normal)], normal)
	)
# extract mRNA abundance matrix as well as dictionary that lists gene symbol and gene names
rrna_matrix <- rrna$abundance
rrna_dict <- rrna$dictionary
# calculate tpm 
rrna_matrix_1 <- calculate_tpm(rrna_matrix)

# create design file with normal 
design_normal <- rbind(
	rdesign[,c('Sample_ID','NMF_3')], 
	data.frame(
		Sample_ID = normal$Sample_ID,
		NMF_3 = 'normal'
		)
	)

pkm <- create_rahbt_metabolism_gene_boxplot(
	rna_matrix = rrna_matrix_1, 
	design = design_normal, 
	gene = 'PKM', 
	fc = rde['PKM','log2FoldChange'], 
	pvalue = rde['PKM','padj'], 
	ylimit = c(0,5000), 
	yat = seq(0,5000,2500)
	)
pgk1 <- create_rahbt_metabolism_gene_boxplot(
	rna_matrix = rrna_matrix_1, 
	design = design_normal, 
	gene = 'PGK1', 
	fc = rde['PGK1','log2FoldChange'], 
	pvalue = rde['PGK1','padj'], 
	ylimit = c(0,2000), 
	yat = seq(0,2000,1000)
	)
idh1 <- create_rahbt_metabolism_gene_boxplot(
	rna_matrix = rrna_matrix_1, 
	design = design_normal, 
	gene = 'IDH1', 
	fc = rde['IDH1','log2FoldChange'], 
	pvalue = rde['IDH1','padj'], 
	ylimit = c(0,2500), 
	yat = seq(0,2500,1000)
	)
fh <- create_rahbt_metabolism_gene_boxplot(
	rna_matrix = rrna_matrix_1, 
	design = design_normal, 
	gene = 'FH', 
	fc = rde['FH','log2FoldChange'], 
	pvalue = rde['FH','padj'], 
	ylimit = c(0,200), 
	yat = seq(0,200,100)
	)
glt8 <- create_rahbt_metabolism_gene_boxplot(
	rna_matrix = rrna_matrix_1, 
	design = design_normal, 
	gene = 'SLC2A8', 
	fc = rde['SLC2A8','log2FoldChange'], 
	pvalue =  rde['SLC2A8','padj'], 
	ylimit = c(0,200), 
	yat = seq(0,200,100)
	)
create.multipanelplot(
	list(pkm, pgk1, idh1, fh, glt8),
	layout.width = 5,
	layout.height = 1,
	width = 24,
	height = 8,
	filename = 'figures/SupplementaryFigure3N_RAHBT_metabolism_genes_boxplot.png',
	xlab.label = 'RNA Subtypes',
	resolution = 300
	)


