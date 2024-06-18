# -*- coding: utf-8 -*-
"""
Comparison of ERA5 runoff to upscaled ERA5 runoff, and calibration
"""
from S2SEnsemble_class import S2S_one, read_S2S_cf
from AverageERA5 import ERA5_averaged
from readera_runoff_series import extract_series
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
from Funcs import figures
import seaborn as sns
from Distrib import qq_mapping
import datetime
from damop_model import *
import DefaultParams as D
from matplotlib.dates import DateFormatter
from scipy.optimize import curve_fit


def linear(x, m, c):
    '''
    Linear function

    Parameters
    ----------
    x : array
        x values to calculate linear function for.
    m : float
        gradieint.
    c : float
        y-intercept.

    Returns
    -------
    array
        evaluated function values.

    '''
    return m*x + c

def scale_runoff(ro):
    '''
    Calibrates upscaled ERA5 data to ERA5 data

    Parameters
    ----------
    ro : TYPE
        DESCRIPTION.

    Returns
    -------
    ro : TYPE
        DESCRIPTION.

    '''
    ro = linear(ro, res[0], res[1])
    return ro



sns.set_context('talk')
sns.set_style('whitegrid')

dstart = datetime.date(2017, 6, 1)
dend  = datetime.date(2018, 5, 31)
#reading downscaled ERA5 data
ERA = ERA5_averaged(dstart, dend)

fpath = './dataall/'
fstem = 'japan_ERA5land.'
#reading original ERA5 data
dayarr, timarr, runoffarr = extract_series(fpath, fstem, D.lonpick, D.latpick, 
                                           dstart, dend, plot = 0)
#getting parameters for linear fit 
res_2 = curve_fit(linear, ERA.runoffarr, runoffarr)
res = res_2[0]


if __name__ == '__main__':
    fig1, axs1 = plt.subplots(1, 2, figsize = (14, 5), sharey = True)
    axs1[0].scatter(ERA.runoffarr, runoffarr, color = 'mediumblue', 
                    marker = 'x')
    axs1[1].scatter(scale_runoff(ERA.runoffarr), runoffarr,
                    color = 'mediumblue', marker = 'x')
    for i in axs1:
        i.plot([0, 0.035], [0, 0.035], color = 'black', zorder = -2)
        i.set_xlim(0, 0.035)
        i.set_ylim(0, 0.035)
    axs1[0].set_ylabel('ERA5 Runoff, m')
    axs1[0].set_xlabel('ERA5 Upscaled Runoff, m')
    axs1[1].set_xlabel('Calibrated ERA5 Upscaled Runoff, m')
    fig1.tight_layout()
    # fig1.savefig('ERA_pairedCalib.pdf')
    
    print(res)
    cov = res_2[1]
    print(np.sqrt(cov))

    fig, axs = plt.subplots(1, 2, sharex = True, figsize = (15, 7))
    date_form = DateFormatter("%m-%d")
    axs[0].xaxis.set_major_formatter(date_form)
    axs[1].xaxis.set_major_formatter(date_form)

    axs[0].plot(timarr, runoffarr, label = 'ERA5')
    axs[0].plot(ERA.timarr, ERA.runoffarr, label = 'Upscaled ERA5')
    axs[0].plot(timarr, linear(ERA.runoffarr, res[0], res[1]), color = 'black', 
             linestyle = '--', label = 'Calibrated Upscaled ERA5')
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Runoff, m')
    axs[0].legend(loc = 'upper left')
    
    axs[1].plot(ERA.timarr, runoffarr - ERA.runoffarr, label = 'Upscaled ERA5')
    axs[1].plot(timarr, runoffarr - linear(ERA.runoffarr, res[0], res[1]),
                color = 'black', 
             linestyle = '--', label = 'Calibrated Upscaled ERA5')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Difference in Runoff from ERA5, m')
    axs[1].legend()
    fig.tight_layout()
    # fig.savefig('ERAcomparison_TS.pdf')

    
    
    
    
    
