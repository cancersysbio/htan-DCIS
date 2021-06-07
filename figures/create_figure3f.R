### create_figure3f.R ##############################################################################
# plot GSEA pathways enriched in RNA subtypes in TBCRC and RAHBT

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(fgsea)

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

# set date
date <- Sys.Date()

### MAIN ##########################################################################################
# set input and output file names depending on specified cohort
tbcrc_ranks_list <- 'plotdata/Figure2g_TBCRC_ranks.RData'
rahbt_ranks_list <- 'plotdata/Figure2g_RAHBT_ranks.RData'
# load ranks list 
load(tbcrc_ranks_list)
tbcrc_ranks_list <- ranks_list
load(rahbt_ranks_list)
rahbt_ranks_list <- ranks_list

# load hallmark pathway
pathways.h <- gmtPathways("~/pathway/h.all.v7.2.symbols.gmt")

# set spot size and colour functions
spot.size.function <- function(x) {abs(x)*3}
spot.colour.function <- function(x) {
    colours <- rep('white', length(x));
    colours[x > 0] <- 'firebrick3';
    colours[x < 0] <- 'dodgerblue3';
    colours[x == 0] <- 'transparent';
    return(colours);
    }

# run gsea
pathways_union <- list()
pathways_same <- list()
gsea_res <- list()
for (i in 1:length(ranks_list)) {
	tmp_t <- fgsea(pathways.h, tbcrc_ranks_list[[i]], minSize=15, maxSize = 500, nperm=1000)
	tmp_t_sig <- tmp_t[tmp_t$padj < 0.05,]
	tmp_r <- fgsea(pathways.h, rahbt_ranks_list[[i]], minSize=15, maxSize = 500, nperm=1000)
	tmp_r_sig <- tmp_r[tmp_r$padj < 0.05,]

	pathways_union[[i]] <- unique(c(tmp_t$pathway, tmp_r$pathway))
	pathways_same[[i]] <- intersect(tmp_t$pathway, tmp_r$pathway)
	gsea_res[[paste0('t',i)]] <- tmp_t
	gsea_res[[paste0('r',i)]] <- tmp_r	
}
allpathways <- unique(unlist(pathways_union))

# reformat plot data
plot_data_es <- do.call(cbind, lapply(gsea_res, '[[', 4))
rownames(plot_data_es) <- gsea_res[[1]]$pathway
plot_data_fdr <- do.call(cbind, lapply(gsea_res, '[[', 3))
rownames(plot_data_fdr) <- gsea_res[[1]]$pathway

# only keep pathways sig in at least one comparison
plot_data_es <- plot_data_es[allpathways,]
plot_data_fdr <- plot_data_fdr[allpathways,]

# cluster es to try and group pathways characteristic of each cluster 
# d <- dist(plot_data_es, method = 'euclidean')
# fit <- hclust(d, method = "ward.D2")
# plot_data_es <- plot_data_es[fit$order,]
# plot_data_fdr <- plot_data_fdr[fit$order,]
i <- order(plot_data_fdr[,1], plot_data_fdr[,2], plot_data_fdr[,3],
	plot_data_fdr[,4], plot_data_fdr[,5], plot_data_fdr[,6])
plot_data_es <- plot_data_es[i,]
plot_data_fdr <- plot_data_fdr[i,]

# set yaxis lab 
yaxis.lab <- rownames(plot_data_es)
yaxis.lab <- gsub('HALLMARK_', '', yaxis.lab)
yaxis.lab <- gsub('_',' ', yaxis.lab)

# set covariate to show clusters and cohort
clust_cov <- rep(default.colours(3), each = 2)
cohort_cov <- rep(c('violetred3','turquoise3'), 3)

covariate <- list(
	rect = list(
		col = 'black',
		fill = clust_cov,
		lwd = 1.5
		),
	rect = list(
		col = 'black',
		fill = cohort_cov,
		lwd = 1.5
		)
	)

cov_grob <- covariates.grob(
	covariates = covariate,
	ord = 1:6,
	side = 'top'
	)

cov_legend <- list(
	legend = list(
		colours = default.colours(3),
		labels = c(
			expression('ER'['low']),
			expression('quiescent'),
			expression('ER'['high'])
			),
		title = 'RNA\nSubtypes'
		),
	legend = list(
		colours = c('violetred3','turquoise3'),
		labels = c('TBCRC','RAHBT'),
		title = 'Cohort'
		)
	)

legend_grob <- legend.grob(
	legends = cov_legend,
	label.cex = 1.8,
	title.cex = 1.8
	)

# create dotmap
key_sizes <- seq(1, -1, -0.5);
create.dotmap(
		x = plot_data_es,
		bg.data = -log10(plot_data_fdr),
		filename = 'figures/Figure2g_gsea_dotmap_fdr_ordered.svg',
		xaxis.rot = 90,
		yaxis.cex = 1.2,
		yaxis.tck = 0,
		xaxis.tck = 0,
		spot.colour.function = spot.colour.function,
		spot.size.function = spot.size.function,
		yaxis.lab = yaxis.lab,
		xaxis.lab = rep('', ncol(plot_data_es)),
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
		at = c(0, seq(1, 3, 0.1)),
		colourkey.labels.at = 0:3,
		colourkey.cex = 1.8,
		axis.top = 1.5,
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
				cex = 1.5,
				adj = 1.5,
				fontface = 'bold'
				),
			title = expression(bold("Effect\nSize:")),
			cex.title = 1.8,
			background = 'white',
			padding = 8
			),
		legend = list(
			top = list(
				fun = cov_grob
				),
			inside = list(
				fun = legend_grob,
				x = 1.05,
				y = 1.05
				)
			),
		width = 11,
		height = 12,
		bottom.padding = 3,
		left.padding = 5,
		right.padding = 20,
		top.padding = 5,
		resolution = 300
		)

### ENRICHMENT ANALYSIS ###########################################################################
# calculate statistics 
res <- list()
for (i in 1:ncol(plot_data_fdr)) {
	res[[i]] <- which(plot_data_fdr[,i] < 0.05)
}

overlap1 <- intersect(names(res[[1]]), names(res[[2]]))
overlap1_dir <- sign(plot_data_es[overlap1,1]) == sign(plot_data_es[overlap1,2])
p1 <- phyper(
	q = sum(overlap1_dir),
	m = length(res[[1]]),
	k = length(res[[2]]),
	n = length(pathways.h)-length(res[[1]]),
	lower.tail = FALSE
	)

overlap2 <- intersect(names(res[[3]]), names(res[[4]]))
overlap2_dir <- sign(plot_data_es[overlap2,3]) == sign(plot_data_es[overlap2,4])
p2 <- phyper(
	q = sum(overlap2_dir),
	m = length(res[[3]]),
	k = length(res[[4]]),
	n = length(pathways.h)-length(res[[3]]),
	lower.tail = FALSE
	)


overlap3 <- intersect(names(res[[5]]), names(res[[6]]))
overlap3_dir <- sign(plot_data_es[overlap3,5]) == sign(plot_data_es[overlap3,6])
p3 <- phyper(
	q = sum(overlap3_dir),
	m = length(res[[5]]),
	k = length(res[[6]]),
	n = length(pathways.h)-length(res[[5]]),
	lower.tail = FALSE
	)


overlap_res <- data.frame(
	C1 = p1,
	C2 = p2,
	C3 = p3
	)

write.table(
	overlap_res,
	file = 'figures/Figure2g_overlap_pvalue.txt',
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)


