########################################################################
############### Plotting output of multiple sensitivity tests #######################
########################################################################

# by A Cluett, 6/28/2021
# adapted by H Holtzman, Spring 2023

## Load packages

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(seas)
library(Rmisc)


## Start function

plotsimulations <- function(simulation){
  
  setwd(simulation)
  
  # Make list of all runs in folder
  
  runfiles <- list.files(path = ".", pattern = "surface", all.files = FALSE,
                         full.names = FALSE, recursive = FALSE,
                         ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  Runs <- lapply(runfiles, read.delim, sep ="",  header = T, stringsAsFactors =F)
  Runs <- setNames(Runs, substr(runfiles, 9,nchar(runfiles)-4))
  
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
  
  Medians <- lapply(Runs, doy_median)
  Means <- lapply(Runs, doy_mean)

  
  ## Calculate average ice free season isotope values
  
  # This function will calculate average values for days with no lake ice (all months)
  icefree_lakeiso_average  <- function(x) {x[which(x$fice == 0),] %>% group_by(YEAR) %>% dplyr::summarize(d18O=mean(d18O), d2H=mean(d2H),dx=mean(dx))}
  
  # Apply average functions
  icefree_isotopes <- lapply(Runs, icefree_lakeiso_average)

  ### Plot boxplots of average summer lake water isotope values
  
  allisos = bind_rows(icefree_isotopes, .id="id")$d2H
  boxplot(allisos, ylim=c(-110,-40), las=1, axes=F)
  points(mean(allisos), col='black', pch=19, cex=1.5)
  
  ## Calculate overall ice-free season d2H
  print(paste("The mean of",substring(simulation,3,nchar(simulation)),"is",mean(allisos)))
  # print(paste("The median of",substring(simulation,3,nchar(simulation)),"is",median(allisos)))

  # Move back up to parent directory
  
  setwd("..")
}

## Set WD to folder of files from multiple runs
setwd("C:/Users/hlhol/Documents/Graduate_School/PRYSM_Masterclass/Imandra_hypercube/Final Sensitivity Tests")
tests = list.dirs(recursive=FALSE)
tests = tests[c(1,3,4,6,7,8,9,12,16,18)]
tests = tests[c(7,10,9,3,8,1,2,6,5,4)]

# Set up boxplots
simulations = tests
par(mfrow=c(1,length(simulations)), mar = c(0.5, 0.5, 0.5, 0.5), oma=c(1,1,1,4))

## LOOP THROUGH ALL TESTS
for (tt in 1:length(simulations)){
  plotsimulations(simulations[tt])
}

# Add y axis (on right side)
axis(side=4, las=1, cex.axis=1.25)
