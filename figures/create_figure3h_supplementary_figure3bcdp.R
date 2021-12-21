### plot_figure2h.R ##############################################################################
# create KM plot

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(survival)

source("functions/create_KM_plot.R")
setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

# set date
setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/DCIS/manuscript')

### MAIN ##########################################################################################
# read in tbcrc design file 
tdesign <- read.delim(
	'../NMF/2021-02-02_TBCRC_design_RNA_CNA_icluster.tsv',
	as.is = TRUE
	)
tdesign$event <- (tdesign$Sample_Type != 'DCIS_only')*1
# add subtype names
tdesign$Subtype <- tdesign$RNA_3
tdesign[tdesign$RNA_3 == 1,'Subtype'] <- 'ER-'
tdesign[tdesign$RNA_3 == 2,'Subtype'] <- 'quiescent'
tdesign[tdesign$RNA_3 == 3,'Subtype'] <- 'ER+'

# create KM plot 
create.km.plot(
	survival.object = Surv(tdesign$Mo_FU, tdesign$event),
	patient.groups = tdesign$Subtype,
	covariates = tdesign$Treatment,
	filename = file.path('figures', 'Figure2h_TBCRC_RNA3_KMplot.svg'),
	line.colours = default.colours(3)[c(1,3,2)],
	# risk.labels = c(
	# 	expression(bold('ER'['low'])),
	# 	expression(bold('ER'['high'])),
	# 	expression(bold('quiet'))
	# 	),
	# key.groups.labels = c(
	# 	expression(bold('ER'['low'])),
	# 	expression(bold('ER'['high'])),
	# 	expression(bold('quiescent'))
	# 	),
	#statistical.method = 'cox',
	#predefined.p = '',
	#predefined.hr = 2.6,
	#predefined.p = expression("HR"['ER- vs quiescent']*'='*4.5%*%10^-3),
	#predefined.hr = expression("P"['ER- vs quiescent']*'=2.6'),
	key.stats.corner = c(1,1),
	key.stats.x.pos = 1, 
	key.stats.y.pos = 0.98, 
	resolution = 300
	)

# create KM plot 
create.km.plot(
	survival.object = Surv(tdesign$Mo_FU, tdesign$event),
	patient.groups = tdesign$PAM50,
	covariates = tdesign$Treatment,
	filename = file.path('figures', 'SupplementaryFigure3B_TBCRC_PAM50_KMplot.png'),
	line.colours = c("#7D26CD",  "#8B0000", "#00C5CD", "#0000ff", "gray"),
	ylab.cex = 2.5,
	key.stats.corner = c(1,1),
	key.stats.x.pos = 1, 
	key.stats.y.pos = 0.98, 
	resolution = 300
	)

# create KM plot 
create.km.plot(
	survival.object = Surv(tdesign$Mo_FU, tdesign$event),
	patient.groups = tdesign$IC11_RNA,
	covariates = tdesign$Treatment,
	filename = file.path('figures', 'SupplementaryFigure3D_TBCRC_IC11_KMplot.png'),
	line.colours = c("#FF5500","#00EE76","#CD3278","#00C5CD","cyan",
				"#8B0000", "#FFFF40", "#0000CD", "#FFAA00",
				"#EE82EE", "#7D26CD"),
	height = 12,
	width = 8,
	#key.groups.title.cex = 1,
	#ylab.cex = 2.2,
	key.stats.corner = c(1,1),
	key.stats.x.pos = 1, 
	key.stats.y.pos = 0.98, 
	resolution = 300
	)

create.km.plot(
	survival.object = Surv(tdesign$Mo_FU, tdesign$event),
	patient.groups = tdesign$CNA_6,
	covariates = tdesign$Treatment,
	filename = file.path('figures', 'SupplementaryFigure4G_TBCRC_CNA6_KMplot.png'),
	ylab.cex = 2.3,
	line.colours = default.colours(12)[-c(1:4,6,7)],
	key.stats.corner = c(1,1),
	key.stats.x.pos = 1, 
	key.stats.y.pos = 0.98, 
	resolution = 300
	)

# read iin pga 
pga <- read.delim(
	'../data/2021-03-26_TBCRC_RAHBT_pga_del_gain.txt',
	as.is = TRUE
	)
colnames(pga) <- gsub('sample','DNA_sample', colnames(pga))
tdesign_pga <- merge(tdesign, pga, by = 'DNA_sample')
tdesign_pga$pga_bin <- tdesign_pga$pga >= median(tdesign_pga$pga)

create.km.plot(
	survival.object = Surv(tdesign_pga$Mo_FU, tdesign_pga$event),
	patient.groups = tdesign_pga$pga_bin,
	covariates = tdesign_pga$Treatment,
	filename = file.path('figures', 'SupplementaryFigure4C_TBCRC_PGA_KMplot.png'),
	statistical.method = 'logrank',
	line.colours = c('dodgerblue','firebrick'),
	risk.labels = c('Low','High'),
	key.groups.labels = c('Low','High'),
	key.stats.corner = c(1,1),
	key.stats.x.pos = 1, 
	key.stats.y.pos = 0.98, 
	resolution = 300
	)


# read in tbcrc design file 
rdesign <- read.delim(
	'../NMF/2021-02-12_RAHBT_design_RNA_CNA_icluster.tsv',
	as.is = TRUE
	)
rdesign$event <- (rdesign$Sample_Type != 'DCIS_only')*1
# remove recurrences before 12 months and only go until 250 months 
rdesign <- rdesign[rdesign$Mo_FU > 12,]
rdesign[rdesign$Mo_FU > 250,'event'] <- 0
rdesign[rdesign$Mo_FU > 250,'Mo_FU'] <- 250
# add subtype names
rdesign$Subtype <- rdesign$RNA_3
rdesign[rdesign$RNA_3 == 1,'Subtype'] <- 'ER-'
rdesign[rdesign$RNA_3 == 2,'Subtype'] <- 'quiescent'
rdesign[rdesign$RNA_3 == 3,'Subtype'] <- 'ER+'

# create KM plot 
create.km.plot(
	survival.object = Surv(rdesign$Mo_FU, rdesign$event),
	patient.groups = rdesign$Subtype,
	covariates = rdesign$Treatment,
	filename = file.path('figures', 'SupplementaryFigure3N_RAHBT_RNA3_KMplot.png'),
	line.colours = default.colours(3)[c(1,3,2)],
	risk.labels = c('ER-','ER+','quiet'),
	resolution = 300
	)

# create KM plot 
create.km.plot(
	survival.object = Surv(rdesign$Mo_FU, rdesign$event),
	patient.groups = rdesign$PAM50,
	covariates = rdesign$Treatment,
	filename = file.path('figures', 'SupplementaryFigure3C_RAHBT_PAM50_KMplot.png'),
	line.colours = c("#7D26CD",  "#8B0000", "#00C5CD", "#0000ff", "gray"),
	ylab.cex = 2.5,
	resolution = 300
	)
cols <- default.colours(12)[-c(1:4,6,7)][c(1:4,6,5)]
cols[5] <- 'gold'
create.km.plot(
	survival.object = Surv(rdesign$Mo_FU, rdesign$event),
	patient.groups = rdesign$CNA_6,
	covariates = rdesign$Treatment,
	filename = file.path('figures', 'SupplementaryFigure4H_RAHBT_CNA6_KMplot.png'),
	ylab.cex = 2.3,
	line.colours = cols,
	resolution = 300
	)
