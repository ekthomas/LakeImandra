# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 13:34:38 2022

@author: hlhol
"""

# Calibrating parameters to temperature observations with Latin Hypercube, adapted from Carrie Morrill 
# (https://github.com/carriemorrill/lake-model-utilities/blob/master/Visualize-Latin-Hypercube-Calibration.ipynb)
# by Hannah Holtzman (2023)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import glob
import os

os.chdir('C:\\Users\\hlhol\\Documents\\Graduate_School\\PRYSM_Masterclass\\Imandra_hypercube')

# Read in temperature observations for 2015
obs = pd.read_csv('Imandra-Temp-Obs.txt', delim_whitespace=True) 
obs['datetime'] = pd.to_datetime(obs['datetime'], format='%m/%d/%Y')  # convert to datetime format
obs['Month'] = obs['datetime'].dt.month      # create Month column from datetime
obs['Year'] = obs['datetime'].dt.year    # create Year column from datetime
obs['Day'] = obs['datetime'].dt.day
obs_mon =obs.groupby(['Month','Day','Year']).mean()  # groupby Year, Month, Day and take average of groups --> daily averages
obs_mon.reset_index(inplace=True)
obs_mon['Julian'] = [117, 215, 296]
obs_mon

# Read in parameter values used for 1000 calibration simulations. Values for parameters chosen from specified ranges using latin
# hypercube.
params = pd.read_csv('lake-params-1000.txt', delim_whitespace=True, header=None) 
#params.drop(params.columns[18], axis=1, inplace=True)
params['cdrn'] = 2.e-3*params[0] + 1.e-3  #  neutral drag coefficient from 1 to 3 e-3
params['eta'] = 0.5*params[1] + 0.2       # shortwave extinction coefficient from 0.2 to 0.5
params['albslush'] = 0.3*params[3] + 0.4       # Slush albedo from 0.4 to 0.7
params['albsnow'] = 0.2*params[2] + 0.7       # Snow albedo from 0.7 to 0.9
params['albsed'] = 0.15*params[6] + 0.05     # Sediment albedo from 0.05 to 0.2
params['condsed'] = 2.0*params[5] + 0.5       # Sediment thermal conductivity from 0.5 to 2.5
params['csed'] = 2.e6*params[4] + 2.e6     # Sediment specific heat capacity from 2e6 to 4e6
params['d18Oa'] = 24*params[7] - 42.1
params['d2ha'] = 188.9*params[8] - 322.8   
params['f'] = 1.0* params[9] + 0.0      
params['melt_ratio'] = 0.95*params[10] + 0.05     
params['rp_ratio_summer'] = 0.95*params[11] + 0.05     
params['rp_ratio_winter'] = 0.95*params[12] + 0.05    
params['rsm_ratio'] = 0.95*params[13] + 0.05    
params['p'] = 100*params[14] + 0
params['s'] = 10*params[15] + 0
params['thresh_spring'] = 20 * params[16]  
params['thresh_fall'] = 20 * params[17]
trials = range(1,1001)
params['trial']=trials
params[['cdrn','eta','albslush','albsnow','albsed','condsed','csed','d18Oa','d2ha','f','melt_ratio','rp_ratio_summer','rp_ratio_winter','rsm_ratio','p','s','thresh_spring','thresh_fall','trial']].head(5)

# read in lake temperatures from 1000 calibration simulations
path = r'C:\Users\hlhol\Documents\Graduate_School\PRYSM_Masterclass\Imandra_hypercube\Final Calibration\1000 simulations'
all_files = glob.glob(os.path.join(path, "surface*.txt"))    
df_from_each_file = (pd.read_csv(f, usecols=range(3,4), header=0, delim_whitespace=True, skiprows=1461, nrows=365) for f in all_files) # Julian day and lake temp
LST = pd.concat(df_from_each_file, axis=1, ignore_index=True)
LST['Julian'] = range(1,366)

fnames = list()   # note: Python does not read in trials in numeric order, so get Python order here
for s in all_files:
    fnames.append(int(s[112:-4]))
LST.head(5)

# read in calibration simulations second water layer (1 meter down) from temp profile docs
all_files_2 = glob.glob(os.path.join(path, "profile-laketemp*.txt"))
df_from_each_file = (pd.read_csv(f, usecols=range(4,5), header=0, delim_whitespace=True, skiprows=1461, nrows=365) for f in all_files_2) # Julian day and lake temp
LST_2 = pd.concat(df_from_each_file, axis=1, ignore_index=True)
LST_2['Julian'] = range(1,366)

fnames_2 = list()   # note: Python does not read in trials in numeric order, so get Python order here
for s in all_files_2:
    fnames_2.append(int(s[129:-4]))
LST_2.head(5)

# Calculate values for objective functions to choose best calibration simulations
def bias(predictions, targets):  # note: I am not using bias as a criteria, just curious
    return (targets.mean() - predictions.mean())
def NashSut(predictions, targets):
    return (1. - ((predictions.values - targets) ** 2).sum()/((targets - targets.mean()) ** 2).sum())
def rsr(predictions, targets):
    return np.sqrt(((predictions.values - targets) ** 2).sum())/np.sqrt(((targets - targets.mean()) ** 2).sum())

plines = ['NSE','rsr','bias']
index = range(0,1000)  
LST_stat = pd.DataFrame(index=index, columns=plines)
# for i in range(0,1000):
#     LST_stat.iloc[i,0]=NashSut(LST.iloc[[214,295],i],obs_mon.iloc[:,3])  # want > 0.75
#     LST_stat.iloc[i,1]=rsr(LST.iloc[[214,295],i],obs_mon.iloc[:,3])  # want <= 0.5
#     LST_stat.iloc[i,2]=bias(LST.iloc[[214,295],i],obs_mon.iloc[:,3])     
for i in range(0,1000):
    LST_stat.iloc[i,0]=NashSut(LST_2.iloc[[116,214,295],i],obs_mon.iloc[:,3])  # want > 0.75
    LST_stat.iloc[i,1]=rsr(LST_2.iloc[[116,214,295],i],obs_mon.iloc[:,3])  # want <= 0.5
    LST_stat.iloc[i,2]=bias(LST_2.iloc[[116,214,295],i],obs_mon.iloc[:,3]) 

LST_stat['trial'] = fnames_2  # note that Python reads in trials not in numeric order
topNSE = np.argsort(LST_stat['NSE'])

LST_stat = pd.merge(LST_stat,params[['trial','cdrn','eta','albslush','albsnow','albsed','condsed','csed','d18Oa','d2ha','f','melt_ratio','rp_ratio_summer','rp_ratio_winter','rsm_ratio','p','s','thresh_spring','thresh_fall']],right_on="trial",
                    left_on="trial",how="outer")
#LST_good = LST_stat[(LST_stat['NSE'] >= 0.8) & (LST_stat['rsr'] <= 0.5)] # these were criteria used for *monthly* streamflow
LST_good = LST_stat[(LST_stat['NSE'] >= 0.85)] # these were criteria used for *monthly* streamflow
#LST_good = LST_stat.sort_values(by=['NSE']).tail(50)
LST_good.head(5)


# plot modeled lake temperatures for 1000 calibration simulations with observations
# plot second layer temps

days = (0,50,100,150,200,250,300,350)
temps= (0,5,10,15,20)
trials = topNSE[-50:]

plt.rcParams["figure.dpi"] = 400
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 4))
for i in range(0,1000):
    ax.plot(LST_2.Julian,LST_2.iloc[:,i],linewidth=0.1,color="black")
ax.plot(obs_mon.Julian,obs_mon.temp,color="blue",linestyle="",marker="o",markersize=10)
ax.set_xlim(0,365)
ax.set_ylim(0,20)
ax.set_xlabel("Day of year in 1995",fontname='Arial', fontsize=16)
ax.set_ylabel("Lake surface temperature ($\mathregular{^o}$C)",fontname='Arial', fontsize=16)
ax.set_title("All parameter sets",fontname='Arial', fontsize=16)
#ax.set_xticklabels(days,fontname= 'Arial',fontsize=14)
#ax.set_yticklabels(temps,fontname= 'Arial',fontsize=14)

#plt.savefig("graph.pdf")


# make scatter plots of parameter values by NSE; each dot is one of the 1000 simulation ensemble

plt.rcParams["figure.dpi"] = 400
fig, axes = plt.subplots(nrows=6, ncols=3, figsize=(9,10))
axes[0,0].scatter(LST_stat['cdrn'], LST_stat['NSE'], color='gray', s=5)
axes[0,0].scatter(LST_good['cdrn'], LST_good['NSE'], color='black', s=5)
axes[0,1].scatter(LST_stat['eta'], LST_stat['NSE'], color='gray', s=5)
axes[0,1].scatter(LST_good['eta'], LST_good['NSE'], color='black', s=5)
axes[0,2].scatter(LST_stat['albsnow'], LST_stat['NSE'], color='gray', s=5)
axes[0,2].scatter(LST_good['albsnow'], LST_good['NSE'], color='black', s=5)
axes[1,0].scatter(LST_stat['albslush'], LST_stat['NSE'], color='grey', s=5)
axes[1,0].scatter(LST_good['albslush'], LST_good['NSE'], color='black', s=5)
axes[1,1].scatter(LST_stat['csed'], LST_stat['NSE'], color='grey', s=5)
axes[1,1].scatter(LST_good['csed'], LST_good['NSE'], color='black', s=5)
axes[1,2].scatter(LST_stat['condsed'], LST_stat['NSE'], color='grey', s=5)
axes[1,2].scatter(LST_good['condsed'], LST_good['NSE'], color='black', s=5)
axes[2,0].scatter(LST_stat['albsed'],LST_stat['NSE'], color='grey', s=5)
axes[2,0].scatter(LST_good['albsed'],LST_good['NSE'], color='black', s=5)
axes[2,1].scatter(LST_stat['d18Oa'],LST_stat['NSE'], color='grey', s=5)
axes[2,1].scatter(LST_good['d18Oa'],LST_good['NSE'], color='black', s=5)
axes[2,2].scatter(LST_stat['d2ha'],LST_stat['NSE'], color='grey', s=5)
axes[2,2].scatter(LST_good['d2ha'],LST_good['NSE'], color='black', s=5)
axes[3,0].scatter(LST_stat['f'],LST_stat['NSE'], color='grey', s=5)
axes[3,0].scatter(LST_good['f'],LST_good['NSE'], color='black', s=5)
axes[3,1].scatter(LST_stat['melt_ratio'],LST_stat['NSE'], color='grey', s=5)
axes[3,1].scatter(LST_good['melt_ratio'],LST_good['NSE'], color='black', s=5)
axes[3,2].scatter(LST_stat['rp_ratio_summer'],LST_stat['NSE'], color='grey', s=5)
axes[3,2].scatter(LST_good['rp_ratio_summer'],LST_good['NSE'], color='black', s=5)
axes[4,0].scatter(LST_stat['rp_ratio_winter'],LST_stat['NSE'], color='grey', s=5)
axes[4,0].scatter(LST_good['rp_ratio_winter'],LST_good['NSE'], color='black', s=5)
axes[4,1].scatter(LST_stat['rsm_ratio'],LST_stat['NSE'], color='grey', s=5)
axes[4,1].scatter(LST_good['rsm_ratio'],LST_good['NSE'], color='black', s=5)
axes[4,2].scatter(LST_stat['p'],LST_stat['NSE'], color='grey', s=5)
axes[4,2].scatter(LST_good['p'],LST_good['NSE'], color='black', s=5)
axes[5,0].scatter(LST_stat['s'],LST_stat['NSE'], color='grey', s=5)
axes[5,0].scatter(LST_good['s'],LST_good['NSE'], color='black', s=5)
axes[5,1].scatter(LST_stat['thresh_spring'],LST_stat['NSE'], color='grey', s=5)
axes[5,1].scatter(LST_good['thresh_spring'],LST_good['NSE'], color='black', s=5)
axes[5,2].scatter(LST_stat['thresh_fall'],LST_stat['NSE'], color='grey', s=5)
axes[5,2].scatter(LST_good['thresh_fall'],LST_good['NSE'], color='black', s=5)

axes[0,0].set_xlabel('cdrn')
axes[0,1].set_xlabel('eta')
axes[0,2].set_xlabel('alb snow')
axes[1,0].set_xlabel('alb slush')
axes[1,1].set_xlabel('csed')
axes[1,2].set_xlabel('condsed')
axes[2,0].set_xlabel('alb sed')
axes[2,1].set_xlabel('d18Oa')
axes[2,2].set_xlabel('d2Ha')
axes[3,0].set_xlabel('f')
axes[3,1].set_xlabel('melt ratio')
axes[3,2].set_xlabel('rp ratio summer')
axes[4,0].set_xlabel('rp ratio winter')
axes[4,1].set_xlabel('rsm ratio')
axes[4,2].set_xlabel('p')
axes[5,0].set_xlabel('s')
axes[5,1].set_xlabel('thresh spring')
axes[5,2].set_xlabel('thresh fall')

axes[0,0].set_xlim(0.001,0.003)
axes[0,1].set_xlim(0.2,0.7)
axes[0,2].set_xlim(0.7,0.9)
axes[1,0].set_xlim(0.4,0.7)
axes[1,1].set_xlim(2e6,4e6)
axes[1,2].set_xlim(0.5,2.5)
axes[2,0].set_xlim(0.05,0.2)
axes[2,1].set_xlim(-42.1,-18.1)
axes[2,2].set_xlim(-322.8,-133.9)
axes[3,0].set_xlim(0,1)
axes[3,1].set_xlim(0.05,1)
axes[3,2].set_xlim(0.05,1)
axes[4,0].set_xlim(0.05,1)
axes[4,1].set_xlim(0.05,1)
axes[4,2].set_xlim(0.100)
axes[5,0].set_xlim(0.100)
axes[5,1].set_xlim(0,20)
axes[5,2].set_xlim(0,20)

axes[0,0].set_ylabel('NSE')
axes[1,0].set_ylabel('NSE')
axes[2,0].set_ylabel('NSE')
axes[3,0].set_ylabel('NSE')
axes[4,0].set_ylabel('NSE')
axes[5,0].set_ylabel('NSE')

plt.tight_layout()
plt.savefig('NSE params.png',bbox_inches='tight', dpi=400)

plt.show()
