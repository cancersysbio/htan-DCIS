### create.km.plot.R ######################################################################
# Description:
#       Produces a single Kaplan-Meier (KM) plot for any number of risk groups, as well as a
#       corresponding risk table.
# Input variables:
#       
# Output variables: If no filename is specified returns a trellis object.  If a filename is
#       specified, prints the trellis object to a .tiff file
#       
create.km.plot <- function(
	survival.object, 
	patient.groups = NA, 
	filename = NULL, 
	xat = NA, 
	yat = seq(0,1,0.2),
	xlimits = NA, 
	ylimits = c(0,1.03),
	xaxis.cex = 2,
	yaxis.cex = 2,
	xlab.label = 'Time (Months)',
	xlab.cex = 2.75,
	ylab.label = 'Estimated Proportion',
	ylab.cex = 2.75,
        xaxis.fontface = 'bold',
        yaxis.fontface = 'bold',
	risk.labels = NA,
	key.groups.labels = levels(as.factor(patient.groups)),
	key.groups.cex = 2,
	risk.label.pos = NA,
	risk.label.fontface = 'bold',
	key.groups.title = NULL, 
	key.groups.title.cex = 1.4,
	key.stats.cex = 1.5, 
	explicit.HR.label = TRUE,
	main = NULL, 
	main.cex = 3.0,
	covariates = NULL,
	stratification.factor = NULL,
	stratification.value = NULL,
	lwd = 2,
	lty = 1,
	censoring.pch.cex = 1.1,
	digits = 2, 
	line.colours = NA, 
	statistical.method = NA, 
	predefined.p = NULL,
	predefined.hr = NULL,
	predefined.hr.ci = NULL,
	predefined.p.statistic.type = 'P',
	ph.assumption.check = "warning",
	cox.zph.threshold = 0.1,
	cox.zph.truncation.thresholds = c(5, 10),
	show.key.groups = NA, 
	show.risktable = TRUE,
	risktable.fontsize = NULL, 
	key.groups.corner = c(0,0), 
	key.groups.x.pos = 0, 
	key.groups.y.pos = 0.01, 
	key.stats.corner = c(1,0), 
	key.stats.x.pos = 1, 
	key.stats.y.pos = 0.01, 
	ylab.axis.padding = 1,
	bottom.padding = 2,  
	top.padding = 0.1,
	right.padding = 0.1,
	left.padding = 0.5,
	return.statistics = FALSE,
	height = 7, 
	width = 7,
	style = 'BoutrosLab', 
	resolution = 1000, 
	size.units = 'in', 
	enable.warnings = TRUE,
	description = NULL,
	use.legacy.settings = FALSE
	) {

	### INPUT VALIDATION ##########################################################################
	# if patient.groups is not specified, set a dummy value by defining a single group to which all
	# subjects belong
	if (all(is.na(patient.groups))) {
		patient.groups <- factor(
			x = rep(1, times = dim(survival.object)[1]),
			levels = '1'
			);
		}

	# if patient.groups is passed as a factor, check that none of the levels of the factor are empty
	if (is.factor(patient.groups) & !(all(levels(patient.groups) %in% unique(patient.groups)))) {
		empty.levels <- list();
		for (this.level in levels(patient.groups)) {
			if (! this.level %in% patient.groups) {
				empty.levels <- c(empty.levels, this.level);
				levels(patient.groups)[levels(patient.groups) == this.level] <- NA;
				}
			}
		warning('The following levels of patient.groups are empty and have been removed: ', empty.levels); 
		}

	# localize the number of risk groups for readability
	ngroups <- length(levels(as.factor(patient.groups)));

	# if xat was not specified, set reasonable values using pretty()
	if (all(is.na(xat))) {
		xat <- pretty(c(0, as.numeric(survival.object)[1:(dim(survival.object)[1])]));
		}

	# verify xat is numeric
	if (!is.numeric(xat)) {
		stop('Invalid value of xat: ', xat);
		}

	# if xlimits was not specified, set them to min(xat) and max(xat)
	if (all(is.na(xlimits))) {
		xlimits <- c(min(xat), max(xat));
		}

	# ensure that the x-axis bounds are correctly ordered
	if (!is.numeric(xlimits) | xlimits[1] < 0 | xlimits[2] < xlimits[1]) {
		stop('Invalid value of xlimits: ', xlimits);
		}

	# ensure that the y-axis bounds are correctly ordered
	if (!is.numeric(ylimits) | ylimits[1] < 0 | ylimits[2] < ylimits[1]) {
		stop('Invalid value of ylimits: ', ylimits);
		}

	# check that xlab.cex is numeric and non-negative
	if (!is.numeric(xlab.cex) | xlab.cex < 0) {
		stop('Invalid value of xlab.cex: ', xlab.cex);
		}

	# check that ylab.cex is numeric and non-negative
	if (!is.numeric(ylab.cex) | ylab.cex < 0) {
		stop('Invalid value of ylab.cex: ', ylab.cex);
		}

	# if risk.labels was not specified, set reasonable value for it
	if (all(is.na(risk.labels))) {
		risk.labels <- ifelse(
			test = rep(
				x = all(1 == ngroups & '1' == levels(as.factor(patient.groups))),
				times = ngroups
				),
			yes = '', 
			no = levels(as.factor(patient.groups))
			);
		}

	# if one of the risklabels is blank, replace it by "na" (as long as it is not the case that '1' is the only level of patient.groups,
	# in which case risklabel should be left blank)
	if (any('' == gsub(' ', '', risk.labels)) && '1' != levels(patient.groups)) {
		for (i in 1:length(risk.labels)){
			if ('' == gsub(' ', '', risk.labels)[i]) {
				risk.labels[i] <- 'na';
				}
			}
		}

	# if there are duplicates in risk.labels, throw an error
	if (any(duplicated(risk.labels))){
		stop('Two or more risk labels have the same name.  Please provide a different label for each risk group.');
		}

	# if one of the key.groups.labels is blank, replace it by "na"
	if (any("" == gsub(" ", "", key.groups.labels))) {
		for (i in 1:length(key.groups.labels)){
			if ("" == gsub(" ", "", key.groups.labels)[i]) {
				key.groups.labels[i] <- 'na';
				}
			}
		}

	# if there are duplicates in risklabels, throw an error
	if (any(duplicated(key.groups.labels))) {
		stop('Two or more key.groups labels have the same name.  Please provide a different label for each risk group to be displayed in the legend.');
		}

	# if risk.label.pos was not specified, set reasonable value for it
	if (is.na(risk.label.pos)) {
		risk.label.pos <- -(xlimits[2]-xlimits[1])/5;
		}

	# check that explicit.HR.label is interpretable as logical
	if (!is.logical(explicit.HR.label)) {
		stop('Invalid value of explicit.HR.label - it must be interpretable as logical: ', explicit.HR.label);
		}

	# check that main.cex is numeric and non-negative
	if (!is.numeric(main.cex) | main.cex < 0) {
		stop('Invalid value of main.cex - it must be numeric: ', main.cex);
		}

	# check that lwd is numeric and non-negative
	if (!is.numeric(lwd) | lwd < 0) {
		stop('Invalid value of lwd - it must be numeric and non-negative: ', lwd);
		}

	# check that censoring.pch.cex is numeric and non-negative
	if (!is.numeric(censoring.pch.cex) | censoring.pch.cex < 0) {
		stop('Invalid value of censoring.pch.cex - it must be numeric and non-negative: ', censoring.pch.cex);
		}

	# check that digits is numeric, non-negative, and an integer
	if (!is.numeric(digits)) {
		stop('Invalid value of digits - it must be numeric: ', digits);
		}
	else if (digits <= 0 | digits != round(digits)) {
		stop('Invalid value of digits - it must a non-negative integer: ', digits);
		}

	# if line.colours is NA, set default colours
	if (any(is.na(line.colours))) {
		if (1 == ngroups) {
			line.colours <- 'black';
			}
		else {
			line.colours <- BoutrosLab.plotting.general::default.colours(number.of.colours = ngroups, palette.type = 'survival');
			}
		}

	# If there are more risk groups than there are colours in line.colours, throw an error
	if (ngroups > length(line.colours)) {
		stop('There are not enough colours to plot each risk group in a different colour. Add more colours to the line.colours parameter or to the survival palette in default.colours.');
		}
  
	# if the user hasn't set it, decide which statistical method to use
	if (is.na(statistical.method) & is.null(predefined.p)) { 

		# if only one group, cannot do any statistical testing
		if (1 == ngroups) { statistical.method <- 'none'; }

		# if more than two risk groups, must use a log-rank test
		else if (ngroups > 2) { statistical.method <- 'logrank'; } 

		# if no censored observations, use a t-test
		else if (all(1 == survival.object[, 'status'])) { statistical.method <- 'ttest'; }
    
		# If no events in one of groups (or both groups), use log-rank
		# Already down to case with only two groups, only need to check first two levels of factor
		else if (
			all(0 == survival.object[1 == as.numeric(as.factor(patient.groups)), 'status']) ||
			all(0 == survival.object[2 == as.numeric(as.factor(patient.groups)), 'status']) ) {
			statistical.method <- 'logrank';
			}

		# otherwise try a Cox model
		else { statistical.method <- 'cox'; }

		}

	# verify that a valid statistical method was selected (predefined.p overrides any statistical tests)
	if (is.null(predefined.p) & !tolower(statistical.method) %in% c('cox', 'logrank', 'ttest', 'none')) {
		stop('Invalid statistical method selected: ', statistical.method);
		}

	# if a statistical method was specified but there is only one risk group, set statistical.method = 'none'
	if (1 == ngroups & (tolower(statistical.method) %in% c('cox', 'logrank', 'ttest'))) {
		statistical.method <- 'none';
		warning('Since a single risk group was specified, no statistical analysis can be performed');
		}

	# if  the 'cox' statistical method was specified but ngroups > 2, set
	# statistical.method = 'none' because cox statistical analysis is only implemented for two risk groups.
	if (is.null(predefined.hr) & 2 < ngroups & 'cox' == tolower(statistical.method)) {
		statistical.method <- 'none';
		warning('Cox statistical analysis is only implemented for 2 risk groups, and you have specified more than 2 risk groups.');
		}

	# if ph.assumption.check does not take on a valid value, throw an error.
	if (! ph.assumption.check %in% c('ignore', 'warning', 'logrank', 'warning.and.plot', 'residual.plot')) {
		stop('The value of ph.assumption.check is invalid');
		}

	# if patient.group is a factor, inform the user of which group is used as baseline
	if (is.factor(patient.groups) & 'cox' == statistical.method) {
		warning('Note that the risk group labelled ', levels(patient.groups)[1], ' will be used as the baseline risk group. If this is not what is desired, please add a "levels" argument to the factor passed to patient.groups, where the first factor refers to the baseline group');
		}

	# if patient.groups is not a factor, coerce it into one.  Additionally, if
	# statistical.method == 'cox', inform the user of which level will be used as baseline
	if (!is.factor(patient.groups)) {
		patient.groups <- as.factor(patient.groups);
		if ('cox' == statistical.method) {
			warning(paste('The argument you passed to patient.groups was not a factor. The risk group labelled ', levels(patient.groups)[1], ' will be used as the baseline risk group'));
			}
		}

	# Set the contrast attribute of patient.groups to contr.treatment when there are two or more levels, to ensure each group is compared 
	# to the baseline group, in the statistical analysis, if appropriate. Note that C() is the contrast setting function, to be distinguished
	# from the concatenation argument, c().
	if (1 < ngroups) {
		patient.groups <- C(patient.groups, contr.treatment);
		}

	# if show.key.groups is NA, set default value
	if (is.na(show.key.groups)) {
		if (1 == ngroups) {
			show.key.groups <- FALSE;
			}
		else {
			show.key.groups <- TRUE;
			}
		}

	# verify that show.key.groups is logical
	if (!is.logical(show.key.groups)) {
		stop('Invalid value of show.key.groups - it must be interpretable as logical: ', show.key.groups);
		}

	# verify that key.groups.corner, key.groups.x.pos and key.groups.y.pos are numeric
	if (!is.numeric(key.groups.corner) | !is.numeric(key.groups.x.pos) | !is.numeric(key.groups.y.pos)){
		stop('Invalid value for location of legend (key.groups) - values must be numeric');
		}

	# verify that key.stats.corner, key.stats.x.pos and key.stats.y.pos are numeric
	if (!is.numeric(key.stats.corner) | !is.numeric(key.stats.x.pos) | !is.numeric(key.stats.y.pos)){
		stop('Invalid value for location of statistical results key (key.stats)');
		}

	# check that ylab.axis.padding is numeric
	if (!is.numeric(ylab.axis.padding)) {
		stop('Invalid value of ylab.axis.padding - it must be numeric: ', ylab.axis.padding);
		}

	# check that bottom.padding is numeric
	if (!is.numeric(bottom.padding)) {
		stop('Invalid value of bottom.padding - it must be numeric: ', bottom.padding);
		}

	# check that top.padding is numeric
	if (!is.numeric(top.padding)) {
		stop('Invalid value of top.padding - it must be numeric:', top.padding);
		}

	# check that right.padding is numeric
	if (!is.numeric(right.padding)) {
		stop('Invalid value of right.padding - it must be numeric:', right.padding);
		}

	# check that left.padding is numeric
	if (!is.numeric(left.padding)) {
		stop('Invalid value of left.padding - it must be numeric:', left.padding);
		}

	# check that return.statistics is logical
	if(!is.logical(return.statistics)) {
		stop('Invalid value of return.statistics - it must be interpretable as logical: ', return.statistics);
		}

	# check that height is numeric
	if (!is.numeric(height)) {
		stop('Invalid value of height - it must be numeric:' , height);
		}

	# check that width is numeric
	if (!is.numeric(width)) {
		stop('Invalid value of width - it must be numeric: ', width);
		}

	# check that resolution is numeric
	if (!is.numeric(resolution)) {
		stop('Invalid value of resolution - it must be numeric: ', resolution);
		}

	# check that enable.warnings is logical
	if(!is.logical(enable.warnings)) {
		stop('Invalid value of enable.warnings - it must be interpretable as logical: ', enable.warnings);
		}

	### CONSTRUCTION OF SURVIVAL DATA FRAME FOR PLOTTING ##########################################
	# create survfit object to be used in plotting
	survfit.object <- survfit(survival.object ~ patient.groups);

	# Create factor g.  This factor will be used to construct the data frame
	# that will contain survival data to be plotted.
	# The levels of g are the same as the levels of patient.groups
	# Suppose there are n1 distinct times at which subjects in the first level 
	# of patient.groups have failures or censorings, n2 for the second level, 
	# ... , and nk for the kth level.
	# In this case, the first level of patient.groups (a character vector) 
	# is repeated in the first n1 positions of g, the second level of patient 
	# groups is repeated in the following n2 positions of g, etc.
	g <- c();

	if (ngroups > 1) {
		for (i in 1:ngroups){
			g <- factor(
				x = c(
					as.character(g),
					rep(
						levels(patient.groups)[i],
						times = survfit.object$strata[[i]]
						)
					),
				levels = levels(patient.groups)
				);
			}
		}

	# If there is only one risk group, g consists of the level of the single
	# patient group repeated n times, where n is the number of distinct times
	# at which an individual either fails or is censored.  The case ngroups == 1
	# is different from ngroups > 1 because when ngroups == 1, survfit.object$strata
	# is not defined.
	else if (1 == ngroups) { 
		g <- factor(
			rep(
				levels(patient.groups),
				times = length(survfit.object$surv)
				),
			levels = levels(patient.groups)
			);
		}

	# This should never happen, but catches weird errors in the number of patient.groups
	else {
		stop('Invalid patient group number detected');
		}

	# Defining data frame that will be used to plot survival curves
	# strata : name of risk group corresponding to failure/censoring occuring at corresponding 'time'
	# surv : survival probability corresponding to risk group named in 'strata', at time point in 'time'
	# time : time at which a given failure or censoring occured
	# n.event : number of events that occured at time 'time' in risk group 'strata'
	# event.indicator : takes value 1 if at least one failure occured at time 'time' in group 'strata', 0 otherwise
	# Note that the last 3*ngroups rows are dummy rows.  In order for plotting to work properly,
	# it is necessary to ensure that for each risk group there are at least two distinct time points at which there
	# are 'failures' and two distinct time points at which there are 'censorings'.  This is due to the fact that 
	# we will want to plot a line for each level of 'new.group', where new.group is a concatenation of the
	# failure/censoring indicator and the risk group name (note that these lines will be visible for failures and
	# invisible for censorings.)  For censorings, the two selected time points are -0.5 and -1.  They were both
	# selected to be negative to ensure that no stray 'censoring' tick marks appear on the plot.   For failures,
	# the first selected time point is 0, to ensure that the survival curve line starts at time 0.  The second
	# selected time point for failures will be defined below - for each risk group, it will correspond to the last
	# time a failure or censoring occured in that risk group.
	x <- data.frame(
		strata = factor(x = c(as.character(g), rep(levels(g), 3)), levels = levels(patient.groups)),
		surv = c(survfit.object$surv, rep(1, 3*length(levels(g)))),
		time = c(survfit.object$time, rep(0, length(levels(g))), rep(-0.5, length(levels(g))), rep(-1, length(levels(g)))),
		n.event = c(survfit.object$n.event, rep(1, length(levels(g))), rep(0, 2*length(levels(g)))),
		event.indicator = c(
			ifelse(0 == survfit.object$n.event, 0, 1),
			rep(1, length(levels(g))),
			rep(0, 2 * length(levels(g)))
			)
		);

	# Define the vector last.time containing the last time at which a failure or censoring occured in each risk group
	# and the vector last.surv containing the survival probability at the time of the last failure or censoring in each group
	if (ngroups > 1) {
		last.time <- as.vector(tapply(x$time, x$strata, max));
		last.surv <- as.vector(tapply(x$surv, x$strata, min));
		}

	# When there is only one risk group, let last.time be the time at which the last failure or censoring occured and
	# last.surv be the survival probability at the time of the last failure or censoring.
	else if (1 == length(levels(patient.groups))) {
		last.time <- max(x$time);
		last.surv <- min(x$surv);
		}

	else {
		stop('Invalid patient group number detected');
		}

	# Adding ngroups dummy rows to the data frame x, corresponding to the dummy rows for 'failures' occuring at the last
	# time a true failure or censoring occured for each risk group.
	x <- data.frame(
		strata = factor(x = c(as.character(g), rep(levels(g), 4)), levels = levels(patient.groups)),
		surv = c(x$surv, last.surv),
		time = c(x$time, last.time),
		n.event = c(x$n.event, rep(1, length(levels(g)))),
		event.indicator = c(x$event.indicator, rep(1, length(levels(g)))),
		new.groups = factor(
			x = c(paste(x$event.indicator, x$strata), paste(rep(1,length(levels(g))), levels(g))), 
			levels = paste(c(rep(0,ngroups),rep(1,ngroups)), levels(x$strata))
			)
		);

	### Statistical Analysis ###########################################################################
	### Predefined p-value (e.g. from log likelihood test), just format the pvalue properly
	
	# Define default value of result.zph.  If the PH assumption fails and a warning on the plot
	# is requested, the value of result.zph will be updated to contain the warning.

	result.zph <- "";
	if (! is.null(predefined.p)||!is.null(predefined.hr)) {
		if (!is.null(predefined.p)) {
			statistical.result.pvalue <- BoutrosLab.plotting.general::display.statistical.result(
				x = predefined.p,
				digits = digits,
				statistic.type = predefined.p.statistic.type
				);
			}

		if (!is.null(predefined.hr)) {
			if (statistical.method != "cox") {
				warning("The statistical method used was not cox, and hence the provided hr value will not be printed");
				}

		        # double check that the predefined.hr.ci is provided
		        if (is.null(predefined.hr.ci)) {
				stop("If the hazard ratio is predefined, then predefined.hr.ci cannot be NULL.");
				}
			else if (length(predefined.hr.ci) != 2) {
				stop("The hazard ratio must be provided with exactly two CI bound values.");
				}

			statistical.result.hr <- substitute( 
				expr = paste('HR'[comparison.group], ': ', this.hr, ' (', this.95l, ',', this.95u, ')', sep=''),
				env = list(
					comparison.group = key.groups.labels[2],
					this.hr = signif(predefined.hr,digits=2),
					this.95l = signif(min(predefined.hr.ci[1],predefined.hr.ci[2]),digits=2),
					this.95u = signif(max(predefined.hr.ci[1],predefined.hr.ci[2]),digits=2)
					)
				);
			}
		}
	
	### Cox modelling analysis ###
	else if ('cox' == statistical.method) {  

		stats <- fit.coxmodel(
			groups = patient.groups,
			survobj = survival.object,
			rounding = 1000,
			other.data = covariates,
			stratification.factor = stratification.factor,
			stratification.value = stratification.value,
			return.cox.model = FALSE
			);

		# define statistical parameters
		this.hr   <- round(stats[1], digits = digits);
		this.95l  <- round(stats[2], digits = digits);
		this.95u  <- round(stats[3], digits = digits);
		this.pval <- stats[4];

		# fit coxmodel using coxph, and then run cox.zph
		if (ph.assumption.check != "ignore") {  
			if (length(levels(patient.groups)) == 1) {
				warning("Only one patient group. Can't check proportional hazard assumptions!!!");
				}
			else {
				cox.model <- fit.coxmodel(
					groups = patient.groups,
					survobj = survival.object,
					rounding = 1000,
					other.data = covariates,
					stratification.factor = stratification.factor,
					stratification.value = stratification.value,
					return.cox.model = TRUE
					);

				ph.failed <- BoutrosLab.statistics.survival::ph.fails(cox.model = cox.model, cox.zph.threshold = cox.zph.threshold, pvalues = FALSE);
				pvalues.zph <- BoutrosLab.statistics.survival::ph.fails(cox.model = cox.model, cox.zph.threshold = cox.zph.threshold, pvalues = TRUE);

				# If cox.zph suggests that PH assumption fails, print warning to the screen			
				if (ph.failed) {
					warning(paste0(
						"The cox.zph test yielded a small pvalue for the following factors: ", 
						names(pvalues.zph[pvalues.zph<cox.zph.threshold]),
						", so the PH assumption may not be valid. Use a non-parametric test (log-rank) to get the pvalue.",
						" If Cox modelling is desired, use stratified model or time-varying covariate.",
						" Talk to Paul if you are not sure what to do."
						));
					
					# print the multi-point HR table
					warning("SEE BELOW FOR MULTI-POINT HR TABLE:");
					if (is.null(covariates)) {
						print(BoutrosLab.statistics.survival::multi.point.HR.table(
							all.groups = patient.groups, 
							all.survtime = as.vector(as.matrix(survival.object)[,'time']), 
							all.survstat = as.vector(as.matrix(survival.object)[,'status']),
							truncation.thresholds = cox.zph.truncation.thresholds
							)); 
						}
					else {
						print(BoutrosLab.statistics.survival::multi.point.HR.table(
							all.groups = patient.groups, 
							all.survtime = as.vector(as.matrix(survival.object)[,'time']), 
							all.survstat = as.vector(as.matrix(survival.object)[,'status']),
							truncation.thresholds = cox.zph.truncation.thresholds,
							covariates = as.data.frame(covariates)
							)); 
						}
					}

				# Print a warning on the plot if it was requested and at least one of the cox.zph pvalues is smaller than cox.zph.threshold
				if (ph.assumption.check %in% c("warning", "warning.and.plot") & ph.failed) {
					result.zph <-  c("WARNING: Small cox.zph", "pvalue for one or more factors.");
					}

				# Produce Schoenfeld residual plots for each beta if they were requested and at least one of the 
				# cox.zph pvalues is smaller than cox.zph.threshold
				if (ph.assumption.check %in% c("residual.plot", "warning.and.plot") & ph.failed) {
					BoutrosLab.plotting.survival::schoenfeld.residual.plots(cox.model = cox.model, filename = filename);	
					} 

				if (ph.assumption.check %in% c('logrank') & ph.failed) {
					statistical.method <- 'logrank';
					statistical.result <- '';
					warning("The cox.zph test yielded a small pvalue so the PH assumption may not be valid. The statistical method was changed from 'cox' to 'logrank', because the logrank test is nonparametric. Talk to Paul if you are not sure what to do.");
					}
				}
			}

		# prepare summary statistical values to return, if requested
		ret.stats <- data.frame(
			statistical.method = statistical.method,
			baseline.group = levels(patient.groups)[1],
			comparison.group = levels(patient.groups)[2],
			pvalue = this.pval,
			hr = this.hr,
			lower.95 = this.95l,
			upper.95 = this.95u
			);

		# prepare values to display results of statistical test
		if (explicit.HR.label) {
			statistical.result.hr <- substitute( 
				expr = paste('HR'[comparison.group], ': ', this.hr, ' (', this.95l, ',', this.95u, ')', sep=''),
				env = list(
					comparison.group = key.groups.labels[2],
					this.hr = ret.stats$hr[1],
					this.95l = ret.stats$lower.95,
					this.95u = ret.stats$upper.95
					)
				);
			statistical.result.hr <- as.expression(statistical.result.hr);
			}
		else {
			statistical.result.hr <- paste('HR: ', ret.stats$hr[1], ' (', ret.stats$lower.95[1], ',', ret.stats$upper.95[1], ')', sep='');
			}
	
		statistical.result.pvalue <- BoutrosLab.plotting.general::display.statistical.result(x = ret.stats$pvalue[1], digits = digits);
		}	

	### t-test analysis ###
	if ('ttest' == statistical.method) { 

		ret.stats <- BoutrosLab.statistics.general::ttest.analysis(values = survival.object[,1], groups = patient.groups, alternative = "two.sided");

		# prepare values to display results of statistical test
		statistical.result <- BoutrosLab.plotting.general::display.statistical.result(ret.stats$pvalue[1], digits = digits)
		}

	### Log rank analysis ###
	if ('logrank' == statistical.method) { 

		# get summary statistics from logrank test (to assess survival differences)
		ret.stats <- BoutrosLab.statistics.survival::logrank.analysis(survival.object, patient.groups);

		# format results of statistical test for display (i.e. convert pvalues into
		# scientific notation and obtain expression with result and descriptive words)
		statistical.result <- BoutrosLab.plotting.general::display.statistical.result(ret.stats$pvalue[1], digits = digits)
		}

	### No statistical analysis ###
	if ('none' == statistical.method) {
		statistical.result <- c(' ');
		ret.stats <- c();
		}


	### Displaying Statistical Results (key.stats) ################################################
	# Note that if the statistical method is 'cox', there are two lines of output for statistical
	# results. Hence, the case where statistical.method == 'cox' is dealt with separately from the
	# case where statistical.method is anything but 'cox'
	if ('cox' == statistical.method) {
		key.stats <- list(
			text = list(
				lab = c(result.zph, statistical.result.hr, statistical.result.pvalue),
				col = 'black',
				cex = key.stats.cex 
				)
			);
		}

	else if (!is.na(statistical.method)) {
		key.stats <- list(
			text = list(
				lab = statistical.result,
				col = 'black',
				cex = key.stats.cex
				)
			);
		}

	else if (!is.null(predefined.p)) {
		key.stats <- list(
			text = list(
				lab = statistical.result.pvalue,
				col = 'black',
				cex = key.stats.cex
				)
			);
		}

	### Invalid statistical method - this case should never occur ###
	else {
		stop('Invalid value of statistical method was specified for key.stats processing');
		}

	### Displaying legend (key.groups) ############################################################
	# If show.key.groups is TRUE, set text and line colours corresponding to each risk group
	if (show.key.groups) {
		key.groups <- list(
			text = list(
				lab = key.groups.labels, 
				col = 'black',
				cex = key.groups.cex
				),
			lines = list(
				lty = lty,
				col=line.colours[1:ngroups],
				lwd = lwd
				),
			title = key.groups.title,
			cex = key.groups.title.cex
			);
		}

	# If show.key.groups is FALSE, set dummy parameter values to ensure key is hidden
	else {
		key.groups <- list(
			text = list(
				lab='', 
				col = 'black',
				cex = 1.5
				)
			);
		}

	### Set-up of risk table parameters ###########################################################
	# Generate vector containing the numbers that will appear in the risk table, moving row
	# by row, left to right, top to bottom.  Read in appropriate risktable display parameters
	# (fontsize, positions, etc.) from text file.  If no parameters available in text file
	# for the corresponding number of patient groups, give warning and turn risk table off
	if (show.risktable){

		number.at.risk <- c(); 

		for (i in 1:ngroups) {
			number.at.risk <- c(
				number.at.risk,
				BoutrosLab.statistics.survival::calculate.number.at.risk(
					survobj = survival.object[patient.groups == levels(as.factor(patient.groups))[i],], 
					cut.points = xat
					)
				);
			};

		# Add dummy value to the end of number.at.risk vector.  This dummy value will be "printed"
		# (invisibly) in the risktable textGrob at position y = 0 to ensure that the size of the
		# textGrob behaves as expected.
		number.at.risk <- c(number.at.risk, '');

		# generate vector specifying y positions of items in risk table.  Note that the '0' at the end
		# of y.pos is included to ensure that the size of the risktable textGrob behaves as expected.
		# If ngroups == 1, y.pos contains only one value (apart from the dummy 0 at the end), and this
		# value is repeated length(xat) times.
		if (1 == ngroups) {
			y.pos <- c(rep(1, times = length(xat)), 0); 
			}

		# If ngroups >= 2, y.pos contains ngroups distinct values (apart from the dummy 0 at the end)
		# and these values are each repeated length(xat) times.
		else if (2 <= ngroups) { 
			y.pos <- c();
			for (i in ngroups:1) {
				y.pos <- c(y.pos, rep((i-1) / (ngroups-1), times = length(xat)));
				}
			}

		# This case should never be reached, as ngroups should always be greater or equal to 1
		else {
			stop('Unexpected number of risk groups');
			}

		# Depending on the number of patient groups, set values of parameters for risktable
		# (i.e. fontsize, risktable.height, risktable.y) by reading in from a file

		# make a list of potential locations for the reference cex file
		data.directories <- paste(.libPaths(), '/BoutrosLab.plotting.survival/', sep = '');
		data.directories <- c('./', data.directories);

		# then search all locations
		file.checks <- file.exists( paste(data.directories, "optimal.risktable.parameters.txt", sep = '/') );

		# check to see if the file was actually found
		if (any(file.checks)) {
			data.directory <- data.directories[ order(file.checks, decreasing = TRUE)[1] ];
			}
		else {
			stop('Unable to find reference risktable parameter file');
			}
		optimal.risktable.parameters <- read.table(
			file = paste(data.directory, "optimal.risktable.parameters.txt", sep = '/'),
			header = TRUE,
			sep = "\t",
			row.names = NULL,
			as.is = TRUE
			);

		# If ngroups is not listed in the optimal.risktable.parameters table, give a warning and 
		# turn off risktable display
		if (!(ngroups %in% optimal.risktable.parameters$ngroups)) {
			warning('Note that show.risktable was set to FALSE because optimal risktable parameters not available for ', ngroups, ' patient groups.  Optimal risktable parameters for this number of patient groups can be added to the file optimal.risktable.parameters.txt, in the inst directory');
			show.risktable <- FALSE;
			}

		# If we reach this point, it means that there are optimal risktable parameter values available
		# for ngroups patient groups, so we will set variables values accordingly
		else {

			# localize the parameterization vector for readability
			rows.to.use <- optimal.risktable.parameters$ngroups == ngroups;

			# extract pre-set parameters overridden if input parameter is not NULL
			if (!is.null(risktable.fontsize)) { 
				fontsize <- risktable.fontsize;
				}
			else {
				fontsize <- optimal.risktable.parameters$fontsize[rows.to.use];
				}
			risktable.height <- grid::unit(ngroups*15, 'points');
			risktable.y <- grid::unit(0.5, 'npc');
			}
		}

	# if show.risktable == FALSE, must set dummy parameter values to ensure risktable is hidden
	else {
		number.at.risk <- c();
		y.pos <- 1;
		risktable.height <- 0;
		risktable.y <- 0;
		fontsize <- 12;
		risk.labels <- ''; # to suppress printing the risk labels when the risktable is turned off
		}

	### Kaplan-Meier Plot #########################################################################
	# Call to BoutrosLab.ploting.general::create.scatterplot to produce KM plot, with a line for 
	# each group (as specified in call to create.km.plot).  Returns a trellis object.
	km.plot <- BoutrosLab.plotting.general::create.scatterplot(
		formula = surv ~ time,
		data = x, 
		groups = x$new.groups, 
		main = main, 
		xlab.label = xlab.label,
		ylab.label = ylab.label,
		main.cex = main.cex,
		xlab.cex = xlab.cex,
		ylab.cex = ylab.cex,
		xlimits = xlimits,
		ylimits = ylimits, 
		xat = xat, 
		yat = yat,
		xaxis.cex = xaxis.cex,
		yaxis.cex = yaxis.cex,
		type = c('s','p'),
		cex = censoring.pch.cex, 
		pch = c(rep('|', times = ngroups),rep('', times = ngroups)), 
		col = line.colours[1:ngroups], 
		lty = c(rep(0, times = ngroups), rep(lty, times = ngroups)), 
		lwd = lwd,
		add.axes = FALSE,
		axis.lwd = 1,
		bottom.padding = 2,
		ylab.axis.padding = ylab.axis.padding,
		left.padding = left.padding,
		right.padding = right.padding,
		top.padding = top.padding,
		xaxis.fontface = xaxis.fontface,
		yaxis.fontface = yaxis.fontface,
		width = width, 
		height = height, 
		style = style,
		legend = list(
			inside = list(
				fun = draw.key,
				args = list(key.groups),
				x = key.groups.x.pos,
				y = key.groups.y.pos,
				corner = key.groups.corner,
				draw = FALSE
				),
			inside = list(
				fun = draw.key,
				args = list(key.stats),
				x=key.stats.x.pos,
				y=key.stats.y.pos,
				corner = key.stats.corner,
				draw = FALSE
				),
			bottom = list(
				fun = grid::textGrob,
				args = list(
					x = c(rep(xat, times = ngroups), xat[1]),
					y = y.pos,
					label = number.at.risk,
					just = 'top',
					default.units = 'native',
					gp = grid::gpar(fontsize=fontsize),
					vp = grid::viewport(
						xscale = xlimits,
						yscale = c(0,1),
						height = risktable.height,
						x = grid::unit(0.5,'native'),
						y = risktable.y 
						)
					)
				),
			bottom = list(
				fun = grid::textGrob,
				args = list(
					x = rep(risk.label.pos, times = ngroups),
					y = unique(y.pos),
					label = c(risk.labels, ''),
					just = c('left','top'),
					default.units = 'native',
					gp = grid::gpar(fontsize = fontsize, fontface = risk.label.fontface),
					vp = grid::viewport(
						xscale = xlimits,
						yscale = c(0,1),
						height = risktable.height,
						x = grid::unit(0.5,'native'),  # unit(0.5,'native'),
						y = risktable.y 
						)
					)
	        		)
			),
		description = description,
		use.legacy.settings = use.legacy.settings,
		);
	km.plot$par.settings$layout.heights$key.top <- 1;
        km.plot$par.settings$layout.heights$key.bottom <- 1;
	### RETURN VALUE ##############################################################################
	# Return value of the function make.survival.plot()
	# If no filename is provided in the call to make.survival.plot, a trellis object is returned,
	# otherwise, a .tiff file is produced and function returns 0/1.
	# If return.statistics is TRUE, returns a list whose first element is a trellis object or 0/1
	# indicator (if a filename was provided) and whose other elements are relevant statistics.  The
	# statistics returned differ depending on the statistical method used.
	if (return.statistics) {
		return( 
			list(
				plot = BoutrosLab.plotting.general::write.plot(
					trellis.object = km.plot, 
					filename = filename, 
					height = height, 
					width = width, 
					size.units = size.units, 
					resolution = resolution, 
					enable.warnings = enable.warnings,
					description = description
					),
				statistics = ret.stats
				)
			);
		}

	else {
		return(
			BoutrosLab.plotting.general::write.plot(
				trellis.object = km.plot, 
				filename = filename, 
				height = height, 
				width = width, 
				size.units = size.units, 
				resolution = resolution, 
				enable.warnings = enable.warnings,
				description = description
				)
			);
		}
	}