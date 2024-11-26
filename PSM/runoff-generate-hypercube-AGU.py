import os
import warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
from math import sqrt, pi, exp

warnings.simplefilter(action='ignore', category=FutureWarning)
np.seterr(divide='ignore')

path= r'/mnt/c/Users/hlhol/Documents/Graduate_School/PRYSM_Masterclass/Imandra_hypercube/Parameter-files'
os.chdir(path)

# Read in bias corrected ERA5 data

met_data = pd.read_csv('met-input-Imandra-MODERN-ERA5-6hrlycorr.txt', sep="\t", header=None)
met_data.columns = ['YEAR', 'MONTH', 'DAY', 'HOUR', 'T2M', 'RH', 'WIND', 'SSRD', 'STRD', 'SP', 'TP']

# Read in monthly precipitation isotope values
piso_data = pd.read_csv('isotope-seasonality.csv', header=None)
piso_data.columns = ['MONTH', 'd2H', 'd18O']

exportfile = 'met-input.txt'
accexportfile = 'acc-met-input.txt'

show_plt = False

# =============================================================================
# =================== Set parameters and scaling factors ======================
# =============================================================================

# Parameters for tuning hydrology
#melt_ratio = 0.1  # fraction of acc precip to runoff per time step
#rp_ratio_summer = 0.75  # fraction of direct precip converted to runoff (vs lost to ET)
#rp_ratio_winter = 0.2  # fraction of direct precip converted to runoff (vs lost to ET)
#rsm_ratio = 0.7  # fraction of snowmelt converted to runoff
glacier_flux = 0.  # mm runoff per ha basin area; set to 0 if no glacier in catchment
#p = 10.  # period for wma calculation; 0 if no smoothing – ??
#s = 4.  # sigma for wma calculation; range of values??
#thresh_spring = 2  # set threshold for number of days above freezing
thresh_spring = round(thresh_spring)
#thresh_fall = 3  # set threshold for number of days above freezing required to allow runoff
thresh_fall = round(thresh_fall)
thresh_avg = False

glacier_2H = min(piso_data['d2H'])  # glacier isotopes default to minimum monthly precip values
glacier_18O = min(piso_data['d18O'])

# Can change months where summer/winter rp ratios apply
rp_ratio_months = [np.repeat(rp_ratio_winter, 4), np.repeat(rp_ratio_summer, 5), np.repeat(rp_ratio_winter, 3)] # based on growing season
rp_ratio_months = np.concatenate(rp_ratio_months)

# =============================================================================
# =========================== Calculate runoff ================================
# =============================================================================

M_tp = met_data['TP']
M_T2M = met_data['T2M']
month = met_data['MONTH'].astype(int)

d2H = piso_data['d2H'].values
d18O = piso_data['d18O'].values

M_d18O = np.zeros(len(month))
M_d2H = np.zeros(len(month))

for i in range(len(month)):
    M_d2H[i] = d2H[(month[i] - 1)]
    M_d18O[i] = d18O[(month[i] - 1)]

M_acc = np.zeros(len(month))
M_runoff = np.zeros(len(month))
M_acc_d2H = np.zeros(len(month))
M_acc_d18O = np.zeros(len(month))
M_runoff_d2H = np.zeros(len(month))
M_runoff_d18O = np.zeros(len(month))



thresh_list_spring = list()
thresh_list_fall = list()

## Make if statement for spring temp threshold
if thresh_avg == False:
    if thresh_spring == 0:
        if_state_spring = 'M_T2M[i] > 273'
    else:
        for j in range(int(thresh_spring)):
            thresh_list_spring.append('M_T2M[i-' + str(j * 4) + ']>273 and M_T2M[i-' + str(j * 4 + 1) + ']>273 and')
        if_state_spring = ' '.join(thresh_list_spring)[0:-4]
else:
    if thresh_spring == 0:
        if_state_spring = 'M_T2M[i] > 273'
    else:
        for j in range(int(thresh_spring)):
            thresh_list_spring.append('mean([M_T2M[i-' + str(j * 4) + '], M_T2M[i-' + str(j * 4 + 1) + ']])>273 and')
        if_state_spring = ' '.join(thresh_list_spring)[0:-4]

## Make if statement for fall temp threshold
if thresh_avg == False:
    if thresh_fall == 0:
        if_state_fall = 'M_T2M[i] < 273'
    else:
        for j in range(int(thresh_fall)):
            thresh_list_fall.append('M_T2M[i-' + str(j * 4) + ']<273 and M_T2M[i-' + str(j * 4 + 1) + ']<273 and')
        if_state_fall = ' '.join(thresh_list_fall)[0:-4]
else:
    if thresh_fall == 0:
        if_state_fall = 'M_T2M[i] < 273'
    else:
        for j in range(int(thresh_fall)):
            thresh_list_fall.append('mean([M_T2M[i-' + str(j * 4) + '], M_T2M[i-' + str(j * 4 + 1) + ']])<273 and')
        if_state_fall = ' '.join(thresh_list_fall)[0:-4]



for i in range(len(month)):
    M_acc[0] = 0
    M_acc_d2H[0] = M_d2H[0]
    M_acc_d18O[0] = M_d18O[0]

    rp_ratio = rp_ratio_months[(month[i] - 1)]

    if ((eval(if_state_spring) and (month[i] in range(0, 7) or month[i] in range(11,13))) or #Spring statement applies Nov-June
         (M_T2M[i] > 273 and month[i] in range(7,9)) or #No threshold for July-Aug
           (month[i] in range(9, 11)) and eval(if_state_fall) == False ): #Fall-winter threshold applies Sept-Dec

        M_runoff[i] = (M_tp[i] * rp_ratio) + (M_acc[i - 1] * melt_ratio * rsm_ratio) + (
                glacier_flux * rsm_ratio)
        if M_runoff[i] < 0.0001:  # to prevent rounding errors with v low run-off amounts
            M_acc[i] = M_acc[i - 1]
            M_runoff_d2H[i] = 'NaN'
            M_acc_d2H[i] = M_acc_d2H[i - 1]
            M_runoff_d18O[i] = 'NaN'
            M_acc_d18O[i] = M_acc_d18O[i - 1]
        else:
            M_acc[i] = M_acc[i - 1] * (1 - melt_ratio)
            M_runoff_d2H[i] = (M_tp[i] * rp_ratio * M_d2H[i] / M_runoff[i]) + \
                              (M_acc[i - 1] * melt_ratio * rsm_ratio * M_acc_d2H[i - 1] / M_runoff[i]) + \
                              (glacier_flux * rsm_ratio * glacier_2H / M_runoff[i])
            M_acc_d2H[i] = M_acc_d2H[i - 1]

            M_runoff_d18O[i] = (M_tp[i] * rp_ratio * M_d18O[i] / M_runoff[i]) + \
                               (M_acc[i - 1] * melt_ratio * rsm_ratio * M_acc_d18O[i - 1] / M_runoff[i]) + \
                               (glacier_flux * rsm_ratio * glacier_18O / M_runoff[i])
            M_acc_d18O[i] = M_acc_d18O[i - 1]

    else:
        M_runoff[i] = 0
        M_acc[i] = M_acc[i - 1] + M_tp[i]
        M_runoff_d2H[i] = 'NaN'
        M_runoff_d18O[i] = 'NaN'

        if M_acc[i - 1] == 0 and M_tp[i] == 0:
            M_acc_d2H[i] = 'NaN'
            M_acc_d18O[i] = 'NaN'
        elif M_acc[i - 1] == 0 and M_tp[i] > 0:
            M_acc_d2H[i] = M_d2H[i]
            M_acc_d18O[i] = M_d18O[i]
        else:
            M_acc_d2H[i] = ((M_acc[i - 1] * M_acc_d2H[i - 1]) / (M_acc[i - 1] + M_tp[i])) + (
                    (M_tp[i] * (M_d2H[i])) / (M_acc[i - 1] + (M_tp[i])))

            M_acc_d18O[i] = ((M_acc[i - 1] * M_acc_d18O[i - 1]) / (M_acc[i - 1] + M_tp[i])) + (
                    (M_tp[i] * (M_d18O[i])) / (M_acc[i - 1] + (M_tp[i])))

    if np.isnan(M_runoff_d18O[i]):
        M_runoff_d18O[i] = d18O[(month[i] - 1)]
        M_runoff_d2H[i] = d2H[(month[i] - 1)]


## Smooth runoff across timesteps

def gtail_wma(arr, period, sigma):
    period = math.ceil(period / 2.) * 2  # period must even, rounded up if odd
    r = range(-int(period / 2), int(period / 2) + 1)
    kernel = np.asarray([1 / (sigma * sqrt(2 * pi)) * exp(-float(x) ** 2 / (2 * sigma ** 2)) for x in r])
    kernel[range(int(period / 2 + 1), int(period + 1))] = 0
    knorm = np.flip(kernel / kernel.sum())
    return np.convolve(arr, knorm, 'same')


M_runoff_wma = gtail_wma(M_runoff, period=p, sigma=s)

M_runoff_d2H_wma = gtail_wma(M_runoff * M_runoff_d2H, period=p, sigma=s) / M_runoff_wma
M_runoff_d18O_wma = gtail_wma(M_runoff * M_runoff_d18O, period=p, sigma=s) / M_runoff_wma

M_runoff_d2H_wma[np.isnan(M_runoff_d2H_wma)] = M_runoff_d2H[np.isnan(M_runoff_d2H_wma)]
M_runoff_d18O_wma[np.isnan(M_runoff_d18O_wma)] = M_runoff_d18O[np.isnan(M_runoff_d18O_wma)]

if show_plt:
    plt.figure(figsize=(20, 12))
    plt.subplot(811)
    plt.plot(np.array(M_T2M))
    plt.axis([0, 15000, 240, 290])
    plt.ylabel("T2M")
    plt.hlines(273, 0, 15000)
    plt.subplot(712)
    plt.plot(M_tp)
    plt.ylabel("Precip")
    plt.axis([0, 15000, 0, max(M_tp)])
    plt.subplot(713)
    plt.plot(M_runoff)
    plt.plot(M_runoff_wma, linestyle='dashed')
    plt.ylabel("Run-off")
    plt.axis([0, 15000, 0, max(M_runoff)])
    plt.subplot(714)
    plt.plot(M_acc)
    plt.ylabel("Acc Precip")
    plt.axis([0, 15000, 0, max(M_acc)])
    plt.subplot(715)
    plt.plot(M_acc_d2H)
    plt.ylabel("Acc d2H")
    plt.axis([0, 15000, -200, -50])
    plt.subplot(716)
    plt.plot(M_runoff_d2H)
    plt.ylabel("Runoff d2H")
    plt.axis([0, 15000, -200, -50])
    plt.subplot(717)
    plt.plot(M_d2H)
    plt.ylabel("Precip d2H")
    plt.axis([0, 15000, -200, -50])

    plt.subplot(818)
    plt.plot(month)
    plt.ylabel("Precip d2H")
    plt.axis([0, 15000, 0, 13])
    plt.show()

# Format and export met-input file #============================================

met_data['RUNOFF'] = M_runoff_wma
met_data['d18OP'] = M_d18O
met_data['d18OR'] = M_runoff_d18O_wma
met_data['d2HP'] = M_d2H
met_data['d2HR'] = M_runoff_d2H_wma

np.savetxt(exportfile, met_data, '%2.3f', delimiter='\t')

# Export other values w/in runoff calc to compare to obs #========================================

met_data['ACC'] = M_acc
met_data['d2HACC'] = M_acc_d2H
met_data['d18OACC'] = M_acc_d18O
np.savetxt(accexportfile, met_data, '%2.3f', delimiter='\t')
