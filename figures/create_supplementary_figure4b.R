### create_supplmentary_figure4b.R ################################################################
# create supplementary figure 4b - forest plot showing HRs of CNAs associated with recurrence

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(survival)

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript/figures')

# set date
date <- Sys.Date()
### TEST CNA SURVIVAL ASSOCIATIONS ################################################################
test_cna_survival_association <- function(segment, calls, design, month_cutoff = NULL) {
	# reformat cna calls 
	calls_df <- data.frame(
		DNA_sample = colnames(calls)[-c(1:9, ncol(calls))],
		cna = unlist(calls[calls$'Unique Name' == segment, -c(1:9, ncol(calls))])
		)
	# binarize cna call 
	calls_df$cna <- (calls_df$cna > 0)*1
	# merge with design 
	plot_data <- merge(design[,c('DNA_sample','Diagnostic_Group','Sample_Type','Treatment','Mo_FU')], calls_df, by = 'DNA_sample')
	# add recurrence column 
	if (is.null(month_cutoff)) {
		# add recurrence column 
		plot_data$event <- (plot_data$Sample_Type != 'DCIS_only')*1
	} else {
		plot_data$event <- (plot_data$Mo_FU < month_cutoff & plot_data$Sample_Type != 'DCIS_only')*1
		plot_data[plot_data$Mo_FU >= month_cutoff,'Mo_FU'] <- month_cutoff
		}
	# test all recurrence 
	all_cox <- coxph(
		Surv(Mo_FU, event) ~ cna + Treatment,
		data = plot_data
		)
	# test dcis recurrence 
	dcis_cox <- coxph(
		Surv(Mo_FU, event) ~ cna + Treatment,
		data = plot_data[plot_data$Diagnostic_Group %in% c('DCIS_with_DCIS_recurrence','DCIS_no_recurrence'),]
		)
	# test ibc recurrence 
	ibc_cox <- coxph(
		Surv(Mo_FU, event) ~ cna + Treatment,
		data = plot_data[plot_data$Diagnostic_Group %in% c('DCIS_with_IBC_recurrence','DCIS_no_recurrence'),]
		)
	# return results 
	results <- data.frame(
		segment = segment,
		descriptor = calls[calls$'Unique Name' == segment,'Descriptor'],
		num = sum(calls_df$cna),
		all_hr = exp(coef(all_cox))['cna'],
		all_l95 = summary(all_cox)$conf.int['cna',3],
		all_u95 = summary(all_cox)$conf.int['cna',4],
		all_pvalue = summary(all_cox)$coefficients['cna',5],
		dcis_hr = exp(coef(dcis_cox))['cna'],
		dcis_l95 = summary(dcis_cox)$conf.int['cna',3],
		dcis_u95 = summary(dcis_cox)$conf.int['cna',4],
		dcis_pvalue = summary(dcis_cox)$coefficients['cna',5],
		ibc_hr = exp(coef(ibc_cox))['cna'],
		ibc_l95 = summary(ibc_cox)$conf.int['cna',3],
		ibc_u95 = summary(ibc_cox)$conf.int['cna',4],
		ibc_pvalue = summary(ibc_cox)$coefficients['cna',5]
		)
	return(results)
	}

### CREATE FOREST PLOT ############################################################################
create_survival_forestplot <- function(type, tbcrc, rahbt, cov = TRUE) {
	# extract type
	tbcrc <- tbcrc[,c(1, 2, grep(type, colnames(tbcrc)))]
	colnames(tbcrc) <- gsub(type, 'tbcrc', colnames(tbcrc))
	rahbt <- rahbt[,c(1, 2, grep(type, colnames(rahbt)))]
	colnames(rahbt) <- gsub(type, 'rahbt', colnames(rahbt))
	# merge stats from both cohorts 
	plot_data <- merge(tbcrc, rahbt, by = c('segment','descriptor'))
	# order by tbcrc hr
	plot_data <- plot_data[order(-plot_data$tbcrc_hr),]
	plot_data$index <- 1:nrow(plot_data)

	# create covariate
	if (cov) {
		cov_col <- rep('firebrick', nrow(plot_data))
		cov_col[grepl('Del', plot_data$segment)] <- 'dodgerblue'

		cov_list <- list(
			rect = list(
				col = 'transparent',
				fill = cov_col
				)
			)
		cov_grob <- covariates.grob(
			covariates = cov_list,
			ord = nrow(plot_data):1,
			side = 'right',
			size = 2
			)
		cov_legend <- list(
			legend = list(
				colours = c('firebrick', 'dodgerblue'),
				labels = c('Gain','Loss')
				)
			)
		cov_legend_grob <- legend.grob(
			legends = cov_legend,
			label.cex = 2,
			title.cex = 2
			)
	} 

	if (cov) {
		# create plots
		tbcrc_plot <- create.scatterplot(
			index ~ log2(tbcrc_hr),
			data = plot_data,
			main = 'TBCRC',
			horizontal = TRUE,
		    x.error.right = log2(plot_data$tbcrc_u95)-log2(plot_data$tbcrc_hr),
		    x.error.left = log2(plot_data$tbcrc_hr)-log2(plot_data$tbcrc_l95),
		    x.error.bar.col = 'black',
		    yat = 1:nrow(plot_data),
		    yaxis.lab = plot_data$descriptor,
		    ylimits = c(0,nrow(plot_data)+1),
		    yaxis.tck = 0,
		    xaxis.tck = 0,
		    yaxis.cex = 1.75,
			xaxis.cex = 1.75,
		    xlab.label = 'Hazard Ratio',
		    ylab.label = 'Recurrent CNAs',
		    xlimits = c(-5.5, 5.5),
		    xat = log2(c(0.05, 0.3, 1, 4, 20)),
		    xaxis.lab = c('0.05','0.3','1.0','4.0','20.0'),
		    legend = list(
		    	inside = list(fun = cov_legend_grob, x = 0.99, y = 0.99, corner = c(1,1))
		    	),
		    abline.v = 0,
		    abline.lty = 2,
		    resolution = 300
			)
		rahbt_plot <- create.scatterplot(
			index ~ log2(rahbt_hr),
			data = plot_data,
			main = 'RAHBT',
			horizontal = TRUE,
		    x.error.right = log2(plot_data$rahbt_u95)-log2(plot_data$rahbt_hr),
		    x.error.left = log2(plot_data$rahbt_hr)-log2(plot_data$rahbt_l95),
		    x.error.bar.col = 'black',
		    yat = 1:nrow(plot_data),
		    yaxis.lab = rep('', nrow(plot_data)),
		    ylimits = c(0,nrow(plot_data)+1),
		    yaxis.tck = 0,
		    xaxis.tck = 0,
		    yaxis.cex = 1.75,
			xaxis.cex = 1.75,
		    ylab.label = '',
		    xlab.label = 'Hazard Ratio',
		    xlimits = c(-5.5, 5.5),
		    xat = log2(c(0.05, 0.3, 1, 4, 20)),
		    xaxis.lab = c('0.05','0.3','1.0','4.0','20.0'),
		    legend = list(
		    	right = list(fun = cov_grob)
		    	),
		    abline.v = 0,
		    abline.lty = 2,
		    resolution = 300
			)
		} else {
			# create plots
			tbcrc_plot <- create.scatterplot(
				index ~ log2(tbcrc_hr),
				data = plot_data,
				main = 'TBCRC',
				horizontal = TRUE,
			    x.error.right = log2(plot_data$tbcrc_u95)-log2(plot_data$tbcrc_hr),
			    x.error.left = log2(plot_data$tbcrc_hr)-log2(plot_data$tbcrc_l95),
			    x.error.bar.col = 'black',
			    yat = 1:nrow(plot_data),
			    yaxis.lab = plot_data$descriptor,
			    ylimits = c(0,nrow(plot_data)+1),
			    yaxis.tck = 0,
			    xaxis.tck = 0,
			    yaxis.cex = 1.75,
			    xaxis.cex = 1.75,
			    xlab.label = 'Hazard Ratio',
			    ylab.label = 'Recurrent CNAs',
			    xlimits = c(-5.5, 5.5),
			    xat = log2(c(0.05, 0.3, 1, 4, 20)),
			    xaxis.lab = c('0.05','0.3','1.0','4.0','20.0'),
			    abline.v = 0,
			    abline.lty = 2,
			    resolution = 300
				)
			rahbt_plot <- create.scatterplot(
				index ~ log2(rahbt_hr),
				data = plot_data,
				main = 'RAHBT',
				horizontal = TRUE,
			    x.error.right = log2(plot_data$rahbt_u95)-log2(plot_data$rahbt_hr),
			    x.error.left = log2(plot_data$rahbt_hr)-log2(plot_data$rahbt_l95),
			    x.error.bar.col = 'black',
			    yat = 1:nrow(plot_data),
			    yaxis.lab = rep('', nrow(plot_data)),
			    ylimits = c(0,nrow(plot_data)+1),
			    yaxis.tck = 0,
			    xaxis.tck = 0,
			     yaxis.cex = 1.75,
			    xaxis.cex = 1.75,
			    ylab.label = '',
			    xlab.label = 'Hazard Ratio',
			    xlimits = c(-5.5, 5.5),
			    xat = log2(c(0.05, 0.3, 1, 4, 20)),
			    xaxis.lab = c('0.05','0.3','1.0','4.0','20.0'),
			    abline.v = 0,
			    abline.lty = 2,
			    resolution = 300
				)
		}

	create.multipanelplot(
		list(tbcrc_plot, rahbt_plot),
		filename = paste0(date, '_', type, '_horizontal_recurrence_forest_plot.png'),
		layout.width = 2,
		layout.height = 1,
		plot.objects.width = c(0.57,0.43),
		height = 11,
		width = 11,
		resolution = 300
		)
	}

### MAIN ##########################################################################################
# read in design 
design <- read.csv(
	'/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/data/RAHBT_TBCRC_combined_design_file_SHS_021121.csv',
	as.is = TRUE,
	header = TRUE
	)
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
calls <- peaks[which(peaks$'Amplitude Threshold' != 'Actual Copy Change Given'),]

# read in cnas associated with coverage 
coverage <- read.delim(
	file.path('..', '..', 'NMF', 'cna_nmf', '2021-03-31_recurrent_cna_associated_with_coverage.txt'),
	as.is = TRUE
	)
coverage$cna <- gsub(' ', '', coverage$cna)
coverage_sig <- coverage[coverage$fdr < 0.05,]

# only keep cnas not associated with coverage 
i <- paste(substr(calls[,1], 1, 3), gsub(' ', '', calls[,2]), sep = '_')
calls <- calls[!i %in% coverage_sig$cna,]

# test cnas - only TBCRC samples
tcna_surv <- do.call(rbind, sapply(
	calls$'Unique Name',
	test_cna_survival_association,
	calls = calls,
	design = design[design$Clustering_TBCRC == 1,],
	month_cutoff = 60,
	simplify = FALSE
	))
tcna_surv$all_fdr <- p.adjust(tcna_surv$all_pvalue, method = 'fdr')
tcna_surv$dcis_fdr <- p.adjust(tcna_surv$dcis_pvalue, method = 'fdr')
tcna_surv$ibc_fdr <- p.adjust(tcna_surv$ibc_pvalue, method = 'fdr')

# test cnas - only RAHBT samples
rcna_surv <- do.call(rbind, sapply(
	calls$'Unique Name',
	test_cna_survival_association,
	calls = calls,
	design = design[design$NEW_Clustering_RAHBT_021121 == 1 & design$Mo_FU > 12,],
	month_cutoff = 60,
	simplify = FALSE
	))
rcna_surv$all_fdr <- p.adjust(rcna_surv$all_pvalue, method = 'fdr')
rcna_surv$dcis_fdr <- p.adjust(rcna_surv$dcis_pvalue, method = 'fdr')
rcna_surv$ibc_fdr <- p.adjust(rcna_surv$ibc_pvalue, method = 'fdr')

create_survival_forestplot(type = 'all', tbcrc = tcna_surv, rahbt = rcna_surv)
create_survival_forestplot(type = 'ibc', tbcrc = tcna_surv, rahbt = rcna_surv)


