### plot_supplementary_figures_2gh.R ##############################################################
# apply TBCRC centroids to RAHBT

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(argparse)
library(genefu)
library(DESeq2)
library(pamr)
library(NMF)
library(vcd)

set.seed(1204)

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

# set date
date <- Sys.Date()
### OBTAIN COMMAND LINE ARGUMENTS #################################################################
parser <- ArgumentParser();

parser$add_argument('-n', '--nmf', type = 'integer', help = 'Number of nmf clusters');

args <- parser$parse_args();

### DEBUG PAMR LISTGENES ##########################################################################
list_genes <- function (fit, data, threshold, fitcv = NULL, genenames = FALSE) {
    x <- data$x[fit$gene.subset, fit$sample.subset]
    if (genenames) {
        gnames <- data$genenames[fit$gene.subset]
    }
    if (!genenames) {
        gnames <- NULL
    }
    geneid <- data$geneid[fit$gene.subset]
    if (!is.null(fit$y)) {
        nc <- length(fit$y)
    }
    if (is.null(fit$y)) {
        nc <- ncol(fit$proby)
    }
    clabs <- colnames(fit$centroids)
    aa <- pamr.predict(fit, x, threshold = threshold, type = "nonzero")
    cen <- pamr.predict(fit, x, threshold = threshold, type = "centroid")
    d <- (cen - fit$centroid.overall)[aa, , drop = FALSE]/fit$sd[aa]
    gene.order <- order(-apply(abs(d), 1, max))
    d <- round(d, 4)
    g <- gnames[aa]
    g1 <- geneid[aa]
    if (is.null(gnames)) {
        gnhdr <- NULL
    }
    if (!is.null(gnames)) {
        gnhdr <- "name"
    }
    if (!is.null(fitcv)) {
        nfold = length(fitcv$cv.objects)
        ind = matrix(F, nrow = nrow(x), ncol = nfold)
        ranks = NULL
        for (ii in 1:nfold) {
            cen = pamr.predict(fitcv$cv.objects[[ii]], x[, -fitcv$folds[[ii]]],
                threshold = 0, type = "centroid")
            dtemp <- (cen - fitcv$cv.objects[[ii]]$centroid.overall)[,
                 , drop = FALSE]/fitcv$cv.objects[[ii]]$sd
            r <- apply(abs(dtemp), 1, max)
            ranks = cbind(ranks, rank(-abs(r)))
            junk = pamr.predict(fitcv$cv.objects[[ii]], x[, -fitcv$folds[[ii]]],
                threshold = threshold, type = "nonzero")
            ind[junk, ii] = T
        }
        av.rank = apply(ranks, 1, mean)
        av.rank = round(av.rank[aa], 2)
        prop = apply(ind[aa, , drop = F], 1, sum)/nfold
    }
    options(width = 500)
    schdr <- paste(clabs, "score", sep = "-")
    if (is.null(fitcv)) {
        res <- cbind(as.character(g1), g, d)[gene.order, , drop = F]
        dimnames(res) <- list(NULL, c(gnhdr, schdr))
    }
    if (!is.null(fitcv)) {
        res <- cbind(as.character(g1), g, d, av.rank, prop)[gene.order,
            , drop = F]
        dimnames(res) <- list(NULL, c("id", gnhdr, schdr, "av-rank-in-CV",
            "prop-selected-in-CV"))
    }
    return(res)
}
### MAIN ##########################################################################################
# set rna nmf dir 
nmfdir <- '/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/NMF/rna_nmf/216samples_with_ribosomal'
datadir <- '/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data'
# load NMF data
load(file.path(nmfdir, '2020-12-17_nmf_results_with_ribosomal.RData'))
#extract basis 
basisdf <- basis(res_dcis$fit[[as.character(args$nmf)]])

# read in design file 
tdesign <- read.delim(
	file.path(nmfdir, '..', '..', '2021-02-02_TBCRC_design_RNA_CNA_icluster.tsv'),
	as.is = TRUE
	)

# load RNA matrix
trna <- load_rnaseq(
	filename = file.path(datadir,'TBCRC_resequenced_merged_n227.raw_reads.csv'),
	design = tdesign
	)
# extract mRNA abundance matrix as well as dictionary that lists gene symbol and gene names
trna_matrix <- trna$abundance
trna_dict <- trna$dictionary

# apply variance stabilizing transformation 
trna_matrix_1 <- apply_variance_stabilizing_transformation(trna_matrix)
# standirize data (not center, matrix has to be positive)
trna_matrix_2 <- t(scale(t(trna_matrix_1),center = FALSE))

# extract genes used in NMF 
trna_matrix_3 <- trna_matrix_2[rownames(basisdf),]

if (args$nmf == 6) {
	tdesign <- tdesign[tdesign$RNA_6 != 1,]
}

# train nearest shruficken centroid classifier
mydata <- list()
mydata[['x']] <- trna_matrix_3[,tdesign$Sample_ID]
mydata[['y']] <- tdesign[,paste0('RNA_', args$nmf)]
mydata[['genenames']] <- rownames(trna_matrix_3)
fit <- pamr.train(mydata)

# cross validate 
fitcv <- pamr.cv(fit, mydata, nfold = 4)
# create confusion matrix
pamr.confusion(fit, threshold = 2.537)

# write to centroid 
thres_genes <- list_genes(fit = fit, data = mydata, threshold = 4.096, genenames = TRUE)
write.table(
	thres_genes,
	file = file.path('..', 'NMF', 'rna_nmf', paste0(date, '_TBCRC_RNA3_centroids_47genes.txt')),
	sep = '\t',
	quote = FALSE
	)

### ADD CENTROIDS TO UMAP #########################################################################
trna_matrix_4 <- trna_matrix_2[thres_genes[,'name'],]
# reformat centroids 
trna3_centroids <- apply(
	thres_genes[,-1],
	2,
	as.numeric
	)
rownames(trna3_centroids) <- thres_genes[,1]
# read in design file with UMAP 
tdesign_umap <- read.delim(
	'../NMF/rna_nmf/216samples_with_ribosomal/2021_01_06_TBCRC_design_matrix_with_NMF_clusters.txt',
	as.is = TRUE
	)
tdesign_umap <- tdesign_umap[match(tdesign$Sample_ID, tdesign_umap$Sample_ID),]
# building a map using a svm, mapping the umap coordinates with the expression
pam_um1 <- kernlab::ksvm(tdesign_umap$UMAP_1~t(trna_matrix_4),scaled=F,kernel="vanilladot" )
pam_um2 <- kernlab::ksvm(tdesign_umap$UMAP_2~t(trna_matrix_4),scaled=F,kernel="vanilladot" )

# here we map the centroid using the trained svm
pam_um1.proj <- kernlab::predict(pam_um1,t(trna3_centroids))
pam_um2.proj <- kernlab::predict(pam_um2,t(trna3_centroids))

# create data frame
centroids <- data.frame(
	UMAP_1=pam_um1.proj,
	UMAP_2=pam_um2.proj,
	RNA3=as.factor(c(1,2,3))
	)

create.scatterplot(
	UMAP_2 ~ UMAP_1,
	data = tdesign_umap,
	filename = 'figures/Figure2e_TBCRC_NMF3_centroids_test.png',
	groups = tdesign_umap$NMF_3,
	col = default.colours(3),
	yaxis.tck = 0,
	xaxis.tck = 0,
	ylimit = c(-10,10),
	xlimit = c(-10,10),
	ylab.label = 'UMAP 2',
	xlab.label = 'UMAP 1',
	add.points = TRUE,
	points.x = centroids$UMAP_1,
	points.y = centroids$UMAP_2,
	points.col = default.colours(3),
	points.pch = 1,
	points.cex = 3,
	key = list(
		text = list(
			c('Cluster 1','Cluster 2','Cluster 3'),
			cex = 1,
			col = 'black'
			),
		points = list(
			pch = 19,
			col = default.colours(3),
			cex = 1
			),
		x = 0.01,
		y = 0.99,
		padding.text = 2
		),
	resolution = 300
	)

### APPLY TO RAHBT ################################################################################
# read in raabit design file 
rdesign <- read.delim(
 	file.path(nmfdir, '..', '..', '2021-02-10_RAHBT_design_RNA_CNA_icluster.tsv'),
  	as.is = TRUE
  	)

# load RAHBT rna
rrna <- load_rnaseq(
  filename = file.path(datadir,'SU_WU_Read_Matrix_n933_20200505.csv'),
  design = rdesign
  )
# extract mRNA abundance matrix as well as dictionary that lists gene symbol and gene names
rrna_matrix <- rrna$abundance
rrna_dict <- rrna$dictionary

# apply variance stabilizing transformation 
rrna_matrix_1 <- apply_variance_stabilizing_transformation(rrna_matrix)
# standirize data (not center, matrix has to be positive)
rrna_matrix_2 <- t(scale(t(rrna_matrix_1),center = FALSE))
# extract genes used in NMF 
rrna_matrix_3 <- rrna_matrix_2[rownames(basisdf),]
# make predictions on rahbt
rahbt_k <- pamr.predict(fit, rrna_matrix_3, threshold = 4.096)

# create data frame
# reformat 
rahbt_k <- data.frame(
  Sample_ID = colnames(rrna_matrix_3),
  k = rahbt_k
  )

# compare to de novo subtypes
plot_data <- merge(
	rdesign[,c('Sample_ID',paste0('RNA_', args$nmf))], 
	rahbt_k, 
	by = 'Sample_ID'
	)
colnames(plot_data) <- c('Sample_ID','RAHBT','TBCRC')

# create mosaic plot
png(paste0('figures/Supplementary_figure_2g_nmf', args$nmf, '_mosaic.png'), res = 300, width = 8, height = 7, units = 'in')
mosaic(~RAHBT+TBCRC, plot_data, shade=TRUE, margins = c(0,3,0,3))
dev.off()

### CALCULATE CORRELATION TO TBCRC CENTROIDS ######################################################
# calculate correlations 
rrna_matrix_4 <- rrna_matrix_3[thres_genes[,'name'],]
rna3_cor <- do.call(rbind, apply(
	rdesign,
	1,
	function(x) {
		id <- x['Sample_ID']
		rna3 <- x['RNA_3']
		# calculate correlation 
		cortmp <- cor(
			unlist(rrna_matrix_4[,id]),
			as.numeric(thres_genes[,paste0(rna3, '-score')]),
			method = 'spearman'
			)
		# return correlation 
		data.frame(
			sample = id,
			rna3 = rna3,
			rna3_cor = cortmp
			)
		}))

# calculate correlation with pam50 
# extract pam50 centroids
pam50_centroid <- pam50$centroids
# find intersecting genes
zpam <- intersect(
	rownames(pam50_centroid),
	rownames(rrna_matrix_2)
	)
rrna_matrix_5 <- rrna_matrix_2[zpam,]
pam_2_centroid <- pam50_centroid[zpam,]
# calculate correlation 
pam50_cor <- apply(
	rrna_matrix_5,
	2,
	function(x) {
		apply(pam_2_centroid, 2, 
			function(y) cor(x, y, method='spearman'))
		})
# plot number of each pam50 
dcis_cor <- do.call(rbind, sapply(
	c('Basal','Her2','LumA','LumB','Normal'),
	function(x) {
		samples <- rdesign[rdesign$PAM50 == x,'Sample_ID']
		data.frame(
			pam50 = x,
			sample = samples,
			pam50_cor = pam50_cor[x, samples]
			)
		},
	simplify = FALSE
	))
rna3_cor <- merge(rna3_cor, dcis_cor, by = 'sample')

# write to file 
write.table(
	rna3_cor, 
	file = file.path('..', 'NMF', 'rna_nmf', paste0(date, '_RAHBT_RNA3_TBCRC_centroid_correlations.txt')),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

create.densityplot(
	list(
		PAM50 = rna3_cor$pam50_cor,
		RNA3 = rna3_cor$rna3_cor
		),
	filename = paste0(date, '_PAM50_vs_RNA3_correlations_densityplot.tiff'),
	col = default.colours(3)[2:3],
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = default.colours(3)[2:3],
                             pch = 21,
                             cex = 1.5,
                             fill = default.colours(3)[2:3]
                             ),
                         text = list(
                             lab = c('PAM50','RNA3')
                             ),
                         padding.text = c(0,5,0),
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

### INSTRINSIC SUBTYPES ###########################################################################
library(biomaRt)

# get entrex ids for gene names 
ensembl <- useDataset("hsapiens_gene_ensembl",mart=useMart("ensembl"))
ids <- getBM(
	attributes = c('entrezgene_id', 'external_gene_name'),
	filters = 'external_gene_name',
	values = rownames(rrna_matrix_1),
	mart = ensembl,
	useCache = FALSE
	)

# reformat intrinsic centroids 
centroids 			<- ssp2003$centroids
i 					<- which(ssp2003$centroids.map$EntrezGene.ID %in% ids[,1])
centroids 			<- centroids[i,]
rownames(centroids) <- ssp2003$centroids.map[i,'EntrezGene.ID']
centroids 			<- centroids[!is.na(rownames(centroids)),]
rownames(centroids) <- ids[match(rownames(centroids), ids[,1]),2]

# keep instrinsic genes 
rrna_matrix_6 <- t(scale(t(rrna_matrix_1[rownames(centroids),])))

# calculate correlations 
intrinsic_cor <- apply(
	rrna_matrix_6,
	2,
	function(x) {
		apply(centroids, 2, function(y) cor(x, y, method = 'spearman'))
		})

# plot number of each pam50 
dcis_cor <- do.call(rbind, sapply(
	c('Basal','Her2','LumA','LumB','Normal'),
	function(x) {
		data.frame(
			sample = design[design$PAM50 == x,'Sample_ID'],
			pam50 = x,
			cor = intrinsic_cor[x, design[design$PAM50 == x,'Sample_ID']]
			)
		},
	simplify = FALSE
	))

#combine rna3 and instrinsic correlations
colnames(rna3_cor) <- c('sample','pam50', 'cor')
cor_data <- rbind(dcis_cor, rna3_cor)


col <- rep('black', nrow(cor_data))
col[cor_data$pam50 == 'Basal'] <- "#7D26CD"
col[cor_data$pam50 == 'Her2'] <- "#8B0000"
col[cor_data$pam50 == 'LumA'] <- "#00C5CD"
col[cor_data$pam50 == 'LumB'] <- "#0000ff"
col[cor_data$pam50 == 'Normal'] <- "gray"

create.boxplot(
	cor ~ pam50,
	data = cor_data,
	add.stripplot = TRUE,
	filename = paste0(date, '_DCIS_RNA3_PAM50_correlation_boxplot.png'),
	ylab.label = 'Correlation',
	#xaxis.lab = rep('', nrow(median_cor)),
	width = 9,
	points.col = col,
	xaxis.cex = 1.2,
	yaxis.cex = 1.2,
	#xat = seq(1.5, 11.5, 2),
	#xaxis.lab = unique(cor_data$pam50),
	height = 7,
	xlab.label = 'Subtype',
	yaxis.tck = 0,
	xaxis.tck = 0,
	ylimits = c(0, 1),
	yat = seq(0, 1, 0.5),
	resolution = 300
	)

### CALCULATE SILHOUETTE SCORE ####################################################################
# calculate correlation between 
instrinsic_dist <- dist(t(rrna_matrix_6), method = 'euclidean')
rna3_dist <- dist(t(rrna_matrix_4), method = 'euclidean')

# convert PAM50 to code 
rdesign$PAM50_code <- rdesign$PAM50
rdesign[rdesign$PAM50 == 'Basal', 'PAM50_code'] <- 1
rdesign[rdesign$PAM50 == 'Her2', 'PAM50_code'] 	<- 2
rdesign[rdesign$PAM50 == 'LumA', 'PAM50_code'] 	<- 3
rdesign[rdesign$PAM50 == 'LumB', 'PAM50_code'] 	<- 4
rdesign[rdesign$PAM50 == 'Normal', 'PAM50_code'] <- 5


# calculate silhouette value 
intrinsic_sil <- silhouette(
	as.numeric(rdesign$PAM50_code),
	dist = instrinsic_dist
	)
rna3_sil <- silhouette(
	as.numeric(rdesign$RNA_3),
	dist = rna3_dist
	)

# reformat plot data 
plot_data <- rbind(
	data.frame(cluster = intrinsic_sil[,1], sil = intrinsic_sil[,3], type = 'PAM50'),
	data.frame(cluster = rna3_sil[,1]+5, sil = rna3_sil[,3], type = 'RNA3')
	)
plot_data$cluster_name <- plot_data$cluster
plot_data[plot_data$cluster == 1,'cluster_name'] <- 'Basal'
plot_data[plot_data$cluster == 2,'cluster_name'] <- 'Her2'
plot_data[plot_data$cluster == 3,'cluster_name'] <- 'LumA'
plot_data[plot_data$cluster == 4,'cluster_name'] <- 'LumB'
plot_data[plot_data$cluster == 5,'cluster_name'] <- 'Normal'
plot_data[plot_data$cluster == 6,'cluster_name'] <- 'ER-'
plot_data[plot_data$cluster == 7,'cluster_name'] <- 'quiet'
plot_data[plot_data$cluster == 8,'cluster_name'] <- 'ER+'

# create boxplot 
create.boxplot(
	sil ~ cluster_name,
	plot_data,
	add.stripplot = TRUE,
	xaxis.rot = 90,
	ylab.label = 'Silhouette Width',
	xlab.label = 'Subtype',
	filename = paste0(date, '_intrinsic_vs_RNA3_silhouette_values_boxplot.tiff'),
	resolution = 300
	)

dtf <- data.frame(PAM50 = intrinsic_sil[,3], RNA3 = rna3_sil[,3])
stats <- wilcox.test(dtf$PAM50, dtf$RNA3, paired = TRUE)
es <- median(dtf$RNA3)-median(dtf$PAM50)

pvalue_sci <- scientific.notation(stats$p.value, type = 'list')

create.boxplot(
	sil ~ type,
	plot_data,
	add.stripplot = TRUE,
	#xaxis.rot = 90,
	ylab.label = 'Silhouette Width',
	xlab.label = 'Subtype Scheme',
	filename = paste0(date, '_intrinsic_vs_RNA3_silhouette_values_boxplot.png'),
	resolution = 300,
	legend = list(
		inside = list(
			fun = draw.key,
			args = list(
				key = list(
					text = list(
						lab = c(
							paste('Effect Size:', round(es, digits = 2)),
							as.expression(substitute(
								base %*% 10^exponent,
								list(base = paste0('P: ', pvalue_sci[[1]]), exponent = pvalue_sci[[2]])
                             	))
							),
                         cex = 1.5
                         )
                     ),
                 x = 0.01,
                 y = 0.99,
                 corner = c(0,1),
                 draw = FALSE
                 )
             )
		)
	)

