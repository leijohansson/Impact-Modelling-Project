# -*- coding: utf-8 -*-
"""
S2S stats and comparison to ERA data. Question 3 plots.
"""
from S2SEnsemble_class import S2S_one
from AverageERA5 import ERA5_averaged
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
from AverageERA5 import logbin_ro

sns.set_palette(sns.crayon_palette(['Brick Red', 'Blue', 'Yellow Orange']),
                n_colors = 3)
sns.set_context('talk')
sns.set_style('whitegrid')


def year_S2S_data(year):
    '''
    Get S2S data available in that year

    Parameters
    ----------
    year : int
        year.

    Returns
    -------
    dataarr : array
        array of all S2S data from that year.
    dstart : datetime.date
        start date of forecast data
    dend : datetime.date
        end date of forecast data

    '''
    name = 'LFPW'
    filetype = 'PF'

    #setting first days of S2S forecast for specified year
    if year == 2018:
        year, month, day = '2018', '06', '07'
    if year == 2017:
        year, month, day = '2017', '06', '01'
    if year == 2019:
        year, month, day = '2019', '06', '06'
    
    #reding in S2S data
    S2S = S2S_one(name, day, month, year)
    ndays = len(S2S.ts)
    filelist = os.listdir(S2S.filedir)
    
    #find all files containing that year
    filelist = [i for i in filelist  if year in i]
    
    #getting year, month, day from filename
    years = [i[14:-7] for i in filelist]
    months = [i[18:-5] for i in filelist]
    days = [i[20:-3] for i in filelist]
    nforecasts = len(filelist)
    nmembers = 50
    
    #setting up array to store data
    dataarr = np.zeros(ndays*nforecasts*nmembers)
    len1 = nmembers*ndays
    dstart = S2S.dstart
    sind = 0
    
    #looping through all S2S files for specified year
    for i in range(len(filelist)):
        S2S = S2S_one(name, days[i], months[i], years[i])
        length = len(S2S.ro_ts.flatten())
        # adding all to a 1D array
        dataarr[sind: sind+length] = S2S.ro_ts.flatten()
        sind += length
    
    #removing extra 0s at the end of the array 
    dataarr = dataarr[:sind]
    dend = S2S.dend
    return dataarr, dstart, dend




#producing quantiles for qq mapping (using 2017 and 2019 data)
#needed for qq_mapping function
allS2S2017, dstart, dend = year_S2S_data(2017)
ERA2017 = ERA5_averaged(dstart, dend)
print(dstart, dend)
allS2S2019, dstart, dend = year_S2S_data(2019)
print(dstart, dend)
ERA2019 = ERA5_averaged(dstart, dend)
ERA = np.concatenate((ERA2017.runoffarr, ERA2019.runoffarr))       
S2S = np.concatenate((allS2S2017, allS2S2019))
quantiles1 = np.arange(0, 1.01, 0.05)
quantiles2 = np.arange(0, 1.01, 0.01)
quantiles_data1 = np.quantile(ERA, quantiles1)
quantiles_data2 = np.quantile(S2S, quantiles2)
matched_quantiles_data1 = np.interp(quantiles2, quantiles1, quantiles_data1)



def qq_mapping(data):
    '''
    Quantile Quantile mapping to upscaled ERA5 runoff

    Parameters
    ----------
    data : arr
        runoff data to map.

    Returns
    -------
    mapped : arr
        calibrated runoff data.

    '''
    mapped = np.interp(data, quantiles_data2, matched_quantiles_data1)
    return mapped

if __name__ == "__main__":
    #Plotting S2S distribution plots (histogram and logbinned data)
    name = 'LFPW'
    year, month, day = '2017', '06', '01'
    filetype = 'PF'    
    fig6, ax6 = plt.subplots(1,2, figsize = (16,5))
    
    bins = np.linspace(0, 0.1, 100)
    
    x, y = logbin_ro(ERA)
    ax6[1].scatter(x, y, label = f'Upscaled ERA5', color = 'black',
                   marker = 'd')

    ax6[0].hist(ERA, bins=bins, density=True, histtype = 'stepfilled', 
                label = 'Upscaled ERA5', color = 'black')

    
    ax6[1].set_yscale('log')
    ax6[1].set_xscale('log')
    x, y = logbin_ro(S2S)
    ax6[1].scatter(x, y, label = f'S2S', marker = '.')
    ax6[0].hist(S2S, bins = bins, density = True, histtype = 'step', 
                label = 'S2S', linewidth = 2)
    
    ax6[1].set_ylabel('Frequency')
    ax6[1].set_xlabel('Runoff, m')
    ax6[0].set_ylabel('Density')
    ax6[0].set_xlabel('Runoff, m')
    

    mapped = qq_mapping(S2S)
    ax6[0].hist(mapped, bins = bins, density = True, histtype = 'step', 
                label = 'Mapped S2S', linewidth = 2)
    
    
    x, y = logbin_ro(mapped)
    ax6[1].scatter(x, y, label = f'Mapped S2S', marker = '.')
    ax6[0].legend()
    ax6[1].legend()
    # fig6.savefig('S2S_Distrib.pdf', bbox_inches = 'tight')
    
    
        