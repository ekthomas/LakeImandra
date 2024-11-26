########################################################################
############### Plot lake model output and calibrate #######################
########################################################################

# by H Holtzman, 2023, adapted from A Cluett, 6/28/2021

## Load packages

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(seas)
library(Rmisc)

## Make custom ggplot theme

theme_ac <- function (base_size = 12, base_family = "Helvetica") {
  theme_bw() %+replace% 
    theme(
      panel.grid.major.x = element_line(color = NA, linewidth = 0.5),
      panel.grid.major.y = element_line(color = NA, linewidth = 0.5),
      panel.grid.minor = element_line(color = NA, linewidth=0.5),
      panel.background = element_rect(fill = NA),
      panel.border = element_rect(color = "black", fill = NA, size = 2),
      axis.line = element_line(color = NA, linewidth = 1, lineend = "square"),
      axis.ticks = element_line(color = "black", size = 0.8, lineend = "square"),
      axis.ticks.length=unit(0.1, "cm"),
      axis.text = element_text(color = "black"),
      legend.position = "none"
    )
}

########################################################################
############################# Format Data ##############################
########################################################################

## Load parameters
setwd("C:/Users/hlhol/Documents/Graduate_School/PRYSM_Masterclass/Imandra_hypercube/Parameter-files/Parameter files")
params = read.delim('lake-params-100-calibration.txt', sep ="",  header = F, stringsAsFactors =F)
params['cdrn'] = 2.e-3*params[1] + 1.e-3
params['eta'] = 0.5*params[2] + 0.2
# params['albslush'] = 0.3*params[4] + 0.4
# params['albsnow'] = 0.2*params[3] + 0.7
params['albslush'] = 0.15*params[4] + 0.4
params['albsnow'] = 0.09*params[3] + 0.7
params['albsed'] = 0.15*params[7] + 0.05
params['condsed'] = 2.0*params[6] + 0.5
params['csed'] = 2.e6*params[5] + 2.e6
params['d18Oa'] = 24*params[8] - 42.1
params['d2ha'] = 188.9*params[9] - 322.8   
params['f'] = 0.9* params[10] + 0.0      
params['melt_ratio'] = 0.95*params[11] + 0.05     
params['rp_ratio_summer'] = 0.95*params[12] + 0.05     
params['rp_ratio_winter'] = 0.95*params[13] + 0.05    
params['rsm_ratio'] = 0.95*params[14] + 0.05    
params['p'] = 100*params[15] + 0
params['s'] = 10*params[16] + 0
params['thresh_spring'] = 20 * params[17]  
params['thresh_fall'] = 20 * params[18]
params['trial']=1:100 # number of simulations

## Set WD to folder of surface files from set of simulations
setwd("C:/Users/hlhol/Documents/Graduate_School/PRYSM_Masterclass/Imandra_hypercube/Final Calibration/100 calibrated simulations")

# Make list of all runs in folder

runfiles <- list.files(path = ".", pattern = "surface", all.files = FALSE,
                       full.names = FALSE, recursive = FALSE,
                       ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

Runs <- lapply(runfiles, read.delim, sep ="",  header = T, stringsAsFactors =F)
Runs <- setNames(Runs, substr(runfiles, 9,nchar(runfiles)-4))
RunOrder <- as.numeric(names(Runs))

## Adjust for leap years to calculate daily averages, calculate dx

Runs <- lapply(Runs, function(x) {x$YEAR <- as.numeric(x$YEAR); x })
Runs <- lapply(Runs, function(x) {x$DAY_leapcrx = ifelse(x$YEAR%%4 != 0 & x$DAY > 59, x$DAY +1, x$DAY); x})
Runs <- lapply(Runs, function(x) {x$dx <- x$d2H - 8 * x$d18O; x})

## Calculate median and mean daily values

doy_median <- function(x) {group_by(x, DAY_leapcrx) %>% dplyr::summarize(tlake=median(tlake), fice=median(fice), hice=median(hice), hsnow=median(hsnow),
                                                                         evap=median(evap), lakelev=median(lakelev), discharge=median(discharge), mix=median(mix), 
                                                                         d18O=median(d18O), d2H=median(d2H),dx=median(dx))}

doy_mean = function(x) {group_by(x, DAY_leapcrx) %>% dplyr::summarize(tlake=mean(tlake), fice=mean(fice), hice=mean(hice), hsnow=mean(hsnow),
                                                                      evap=mean(evap), lakelev=mean(lakelev), discharge=mean(discharge), mix=mean(mix), 
                                                                      d18O=mean(d18O), d2H=mean(d2H),dx=mean(dx))}

doy_mean_CI = function(x) {group_by(x, DAY_leapcrx) %>% dplyr::summarize(tlakelow=as.numeric(CI(tlake, ci = 0.95)[3]),tlakehigh=as.numeric(CI(tlake, ci = 0.95)[1]), ficelow=as.numeric(CI(fice, ci = 0.95)[3]), ficehigh=as.numeric(CI(fice, ci = 0.95)[1]), evaplow=as.numeric(CI(evap, ci = 0.95)[3]), evaphigh=as.numeric(CI(evap, ci = 0.95)[1]), lakelevlow=as.numeric(CI(lakelev, ci = 0.95)[3]), lakelevhigh=as.numeric(CI(lakelev, ci = 0.95)[1]), 
                                                                         dischargelow=as.numeric(CI(discharge, ci = 0.95)[3]), dischargehigh=as.numeric(CI(discharge, ci = 0.95)[1]), mixlow=as.numeric(CI(mix, ci = 0.95)[3]), mixhigh=as.numeric(CI(mix, ci = 0.95)[1]),
                                                                         d18Olow=as.numeric(CI(d18O, ci = 0.95)[3]), d18Ohigh=as.numeric(CI(d18O, ci = 0.95)[1]),
                                                                         d2Hlow=as.numeric(CI(d2H, ci = 0.95)[3]), d2Hhigh=as.numeric(CI(d2H, ci = 0.95)[1]), dxlow=as.numeric(CI(dx, ci = 0.95)[3]), dxhigh=as.numeric(CI(dx, ci = 0.95)[1]))}

Medians <- lapply(Runs, doy_median)
Means <- lapply(Runs, doy_mean)
MeanCI <- lapply(Runs, doy_mean_CI)

## Calculate average ice free season isotope values

# This function will calculate average values for days with no lake ice (all months)
icefree_lakeiso_average  <- function(x) {x[which(x$fice == 0),] %>% group_by(YEAR) %>% dplyr::summarize(d18O=mean(d18O), d2H=mean(d2H),dx=mean(dx))}

# # All year
# lakeiso_average  <- function(x) {x %>% group_by(YEAR) %>% dplyr::summarize(d18O=mean(d18O), d2H=mean(d2H),dx=mean(dx))}

# Apply average functions
icefree_isotopes <- lapply(Runs, icefree_lakeiso_average)
# year_isotopes <- lapply(Runs, lakeiso_average)

## Calibrate to ice-off date range

## Calculate ice-off days for median runs

fmed = bind_rows(Medians, .id="id")
fmed = fmed[c('id','DAY_leapcrx','fice')]
iceoff = data.frame(1:100, DOY=double(100))
for (ii in 1:100) {
  yearf = fmed[which(fmed['id'] == ii),]
  iceoff[ii,2] = yearf[min(which(yearf['fice']==0)),'DAY_leapcrx']
}
good_runs = which(iceoff['DOY'] <= 165 & iceoff['DOY'] >= 145)


# Plot good runs
par(mfrow=c(6,3),
    mar = c(4, 4, 1, 1),
    oma = c(0,0,0,0))

plot(params[, 'cdrn'],iceoff[, 'DOY'], pch=19, ylab="Ice-off DOY", xlab="cdnr", col='grey', cex=0.5, las=1)
points(params[good_runs, 'cdrn'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'eta'],iceoff[, 'DOY'], pch=19, ylab="", yaxt="n", xlab="eta", col='grey', cex=0.5)
points(params[good_runs, 'eta'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'albsnow'],iceoff[, 'DOY'], pch=19, ylab="", yaxt="n", xlab="alb snow", col='grey', cex=0.5)
points(params[good_runs, 'albsnow'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'albslush'],iceoff[, 'DOY'], pch=19, ylab="Ice-off DOY", xlab="alb slush", col='grey', cex=0.5, las=1)
points(params[good_runs, 'albslush'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'csed'],iceoff[, 'DOY'], pch=19, ylab="", yaxt="n", xlab="csed", col='grey', cex=0.5)
points(params[good_runs, 'csed'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'condsed'],iceoff[, 'DOY'], pch=19, ylab="", yaxt="n", xlab="condsed", col='grey', cex=0.5)
points(params[good_runs, 'condsed'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'albsed'],iceoff[, 'DOY'], pch=19, ylab="Ice-off DOY", xlab="alb sed", col='grey', cex=0.5, las=1)
points(params[good_runs, 'albsed'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'd18Oa'],iceoff[, 'DOY'], pch=19, ylab="", yaxt="n", xlab="d18Oa", col='grey', cex=0.5)
points(params[good_runs, 'd18Oa'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'd2ha'],iceoff[, 'DOY'], pch=19, ylab="", yaxt="n", xlab="d2Ha", col='grey', cex=0.5)
points(params[good_runs, 'd2ha'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'f'],iceoff[, 'DOY'], pch=19, ylab="Ice-off DOY", xlab="f", col='grey', cex=0.5, las=1)
points(params[good_runs, 'f'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'melt_ratio'],iceoff[, 'DOY'], pch=19, ylab="", yaxt="n", xlab="melt ratio", col='grey', cex=0.5)
points(params[good_runs, 'melt_ratio'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'rp_ratio_summer'],iceoff[, 'DOY'], pch=19, ylab="", yaxt="n", xlab="rp ratio summer", col='grey', cex=0.5)
points(params[good_runs, 'rp_ratio_summer'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'rp_ratio_winter'],iceoff[, 'DOY'], pch=19, ylab="Ice-off DOY", xlab="rp ratio winter", col='grey', cex=0.5, las=1)
points(params[good_runs, 'rp_ratio_winter'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'rsm_ratio'],iceoff[, 'DOY'], pch=19, ylab="", yaxt="n", xlab="rsm ratio", col='grey', cex=0.5)
points(params[good_runs, 'rsm_ratio'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'p'],iceoff[, 'DOY'], pch=19, ylab="", yaxt="n", xlab="p", col='grey', cex=0.5)
points(params[good_runs, 'p'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 's'],iceoff[, 'DOY'], pch=19, ylab="Ice-off DOY", xlab="s", col='grey', cex=0.5, las=1)
points(params[good_runs, 's'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'thresh_spring'],iceoff[, 'DOY'], pch=19, ylab="", yaxt="n", xlab="thresh spring", col='grey', cex=0.5)
points(params[good_runs, 'thresh_spring'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

plot(params[, 'thresh_fall'],iceoff[, 'DOY'], pch=19, ylab="", yaxt="n", xlab="thresh fall", col='grey', cex=0.5)
points(params[good_runs, 'thresh_fall'],iceoff[good_runs, 'DOY'], col='black', pch=19, cex=0.5)

# Albedo vs. albedo
par(mfrow=c(1,1), mar=c(4,4,0,0))
plot(params[, 'albslush'],params[, 'albsnow'], pch=19, xlab="alb slush", ylab="alb snow", col='grey', las=1, asp=1)
points(params[good_runs, 'albslush'],params[good_runs, 'albsnow'], col='black', pch=19)


###### Calculate ice-free mean isotopes for each simulation

mean_icefreeiso = bind_rows(icefree_isotopes, .id="id")
iso_icefree = data.frame('id'=1:100, 'd18O'=double(100), 'd2H'=double(100), 'dx'=double(100))
for (ii in 1:100) {
  yeari = mean_icefreeiso[which(mean_icefreeiso['id'] == ii),]
  iso_icefree[ii,2] = mean(yeari$d18O)
  iso_icefree[ii,3] = mean(yeari$d2H)
  iso_icefree[ii,4] = mean(yeari$dx)
}

good_runs_1 = iso_icefree[which(abs(iso_icefree['d18O']+12.0) <= 1 & abs(iso_icefree['d2H']+89.7) <= 8),] # ice-free isotopes for good simulations
good_runs_id = good_runs_1$id # simulation number of good simulations

# Simulations that meet both ice-off and isotope targets

comb_ids = c(good_runs, good_runs_id)
best_runs = comb_ids[duplicated(comb_ids)]


##### Calculate annual mean isotopes for each simulation

# mean_yeariso = bind_rows(year_isotopes, .id="id")
# iso_year = data.frame('id'=1:100, 'd18O'=double(100), 'd2H'=double(100), 'dx'=double(100))
# for (ii in 1:100) {
#   yeari = mean_yeariso[which(mean_yeariso['id'] == ii),]
#   iso_year[ii,2] = mean(yeari$d18O)
#   iso_year[ii,3] = mean(yeari$d2H)
#   iso_year[ii,4] = mean(yeari$dx)
# }
# 
# good_runs_2 = iso_year[which(abs(iso_year['d18O']+12.0) <= 1 & abs(iso_year['d2H']+89.7) <= 8),]
# good_runs_id = good_runs_2$id

# both ice-off and isotopes
# comb_ids = c(good_runs, good_runs_id)
# best_runs = comb_ids[duplicated(comb_ids)]

######
## Boxplots of calibrated lake water

par(mfrow=c(1,3), mar = c(1, 6.25, 1, 1), oma=c(1,1,1,1))
boxplot(iso_icefree$d2H[best_runs], las=1, cex.axis=1.5, ylim = range(iso_icefree$d2H[best_runs]))
title(ylab='Lake water isotopes (per mil)', line=4, cex.lab=1.5)
points(mean(iso_icefree$d2H[best_runs]), col='black', pch=19, ylim = range(iso_icefree$d2H[best_runs]))
points(-89.68, col='blue', pch=19, ylim = range(iso_icefree$d2H[best_runs]))
boxplot(iso_icefree$d18O[best_runs], las=1, cex.axis=1.5, ylim=c(-19.25, -5.15))
points(mean(iso_icefree$d18O[best_runs]), col='black', pch=19, ylim=c(-19.25, -5.15))
points(-12.03, col='blue', pch=19, ylim=c(-19.25, -5.15))
boxplot(iso_icefree$dx[best_runs], las=1, cex.axis=1.5, ylim=c(1.3,15.7))
points(mean(iso_icefree$dx[best_runs]), col='black', pch=19, ylim=c(1.3,15.7))
points(6.56, col='blue', pch=19, ylim=c(1.3,15.7))

###### Plot isotopes vs. parameters

param_names = names(params)[19:36]

par(mfrow=c(1,3), oma=c(1,1,1,1), mar=c(4,5,1,1))
for (ii in 1:1) {
  plot(params[,param_names[ii]], iso_icefree$d2H, pch=19, col='grey', xlab=param_names[ii], ylab='d2H')
  points(params[good_runs_id, param_names[ii]],iso_icefree[good_runs_id, 'd2H'], col='black', pch=19)
  
  plot(params[,param_names[ii]], iso_icefree$d18O, pch=19, col='grey', xlab=param_names[ii], ylab='d18O')
  points(params[good_runs_id, param_names[ii]],iso_icefree[good_runs_id, 'd18O'], col='black', pch=19)
  
  plot(params[,param_names[ii]], iso_icefree$dx, pch=19, col='grey', xlab=param_names[ii], ylab='dx')
  points(params[good_runs_id, param_names[ii]],iso_icefree[good_runs_id, 'dx'], col='black', pch=19)
}

###############
#### The parameters that matter for isotopes

par(mfrow=c(2,2), oma=c(1,1,1,1), mar=c(4,5,1,1))

plot(params[,'d2ha'], iso_icefree$d2H, pch=19, col='grey', xlab='d2ha', ylab='d2H', las=1)
points(params[good_runs_id, 'd2ha'],iso_icefree[good_runs_id, 'd2H'], col='black', pch=19)
  
plot(params[,'d18Oa'], iso_icefree$d18O, pch=19, col='grey', xlab='d18oa', ylab='d18O', las=1)
points(params[good_runs_id, 'd18Oa'],iso_icefree[good_runs_id, 'd18O'], col='black', pch=19)

plot(params[,'f'], iso_icefree$d2H, pch=19, col='grey', xlab='f', ylab='d2H', las=1)
points(params[good_runs_id, 'f'],iso_icefree[good_runs_id, 'd2H'], col='black', pch=19)

plot(params[,'rp_ratio_summer'], iso_icefree$d2H, pch=19, col='grey', xlab='summer rp ratio', ylab='d2H', las=1)
points(params[good_runs_id, 'rp_ratio_summer'],iso_icefree[good_runs_id, 'd2H'], col='black', pch=19)


########################################################################
######################## Plot All Runs ###########################
########################################################################

## Plot temperature medians only

BestMedians=Medians[as.character(good_runs)]

Runs <- lapply(Runs, function(x) {x$YEAR <- as.factor(x$YEAR); x })
mycolors = rep("grey",100)
mon <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
mon_breaks <- c(1, 32, 61,92,122,153,183,214, 245, 275, 306, 336)

library(ggnewscale)

t <-ggplot() + geom_line(aes(DAY_leapcrx, tlake, colour=id), data = bind_rows(Medians, .id="id"), alpha=0.1, linewidth=0.8) +
  scale_color_manual(values = rep("grey",100)) +
  new_scale_color() +
  geom_line(aes(DAY_leapcrx, tlake, colour=id), data = bind_rows(BestMedians, .id="id"), alpha=0.99, linewidth=0.8) + 
  scale_color_manual(values = rep("black",100)) +
  labs(x = "Month", y = "Lake surface temperature (°C)") +
  scale_x_continuous(breaks = mon_breaks , labels = mon) + 
  ylim(0,20) +
  geom_rect(aes(xmin = 145, xmax = 165, ymin=-Inf, ymax=Inf), fill = "#FFE0AD", color = NA, alpha = 0.5)

plot_grid(t+theme_ac())

## Plot ice cover medians only

Runs <- lapply(Runs, function(x) {x$YEAR <- as.factor(x$YEAR); x })
mycolors = rep("gray88",100)
mon <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
mon_breaks <- c(1, 32, 61,92,122,153,183,214, 245, 275, 306, 336)

fi <- ggplot(bind_rows(Runs, .id="id"), aes(DAY_leapcrx, fice, colour=id)) + geom_rect(aes(xmin = 135, xmax = 166, ymin=0, ymax=1), fill = "#FFE0AD", color = NA, alpha = 0.9) +
  geom_line(aes(DAY_leapcrx, fice, colour=id), data = bind_rows(Medians, .id="id"), alpha=1, linewidth=0.8) + labs(x = "Month", y = "Fraction Lake Ice") +
  scale_x_continuous(breaks = mon_breaks , labels = mon)  +
  scale_color_manual(values = mycolors)

plot_grid(fi+theme_ac())