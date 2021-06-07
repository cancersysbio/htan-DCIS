## create_supplementary_figure3l.R ################################################################
# create boxplot of PGR mRNA abundance stratified by RNA subtypes

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(DESeq2)

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

date <- Sys.Date()
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

### MAIN ##########################################################################################
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
# reformat plot data
pgr_data <- data.frame(
	rna = unlist(rna_matrix_1['PGR',]),
	Sample_ID = colnames(rna_matrix_1)
	)
pgr_data <- merge(pgr_data, design[,c('Sample_ID','NMF_3')], by = 'Sample_ID')

stats <- kruskal.test(rna ~ NMF_3, data = pgr_data)
pvalue_sci <- scientific.notation(stats$p.value, digits = 2, type = 'list');

pgr_data$NMF_3 <- factor(pgr_data$NMF_3)

# create boxplot 
create.boxplot(
	rna ~ NMF_3,
	data = pgr_data,
	add.stripplot = TRUE,
	#ylimit = c(0,1.05),
	#yat = seq(0, 1,0.2),
	ylab.axis.padding = 2,
	xaxis.lab = c(
		expression(bold('ER'['low'])),
		expression(bold('quiescent')),
		expression(bold('ER'['high']))
		),
	xaxis.cex = 1.5,
	#xaxis.rot = 45,
	ylab.cex = 1.8,
	xlab.cex = 1.8,
	yaxis.cex = 1.5,
	ylimits = c(0, 80),
	yat = seq(0, 80, 20),
	#xaxis.cex = 1.2,
	#top.padding = 3,
	ylab.label = 'PGR TPM',
	#ylab.label = 'PGR',
	xlab.label = 'RNA Subtypes',
	filename = paste0('figures/SupplementaryFigure3l_PGR_boxplot.png'),
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
			x = 0.01,
			y = 0.99,
			corner = c(0,1)
			),
	resolution = 300
	)