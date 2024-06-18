# -*- coding: utf-8 -*-
"""
Plots for Question 1
"""

#Question 1
from Funcs import figures
import numpy as np
import matplotlib.pyplot as plt
from readera_runoff_series import *
import DefaultParams as D
from damop_model import *
from matplotlib.dates import DateFormatter
import seaborn as sns

#plots for question 1
if __name__ == '__main__':
    sns.set_palette(sns.crayon_palette(['Brick Red', 'Yellow Orange', 'Green', 'Vivid Violet', 'Wild Strawberry']), n_colors = 6)
    sns.set_context('talk')
    
    
    fpath = './dataall/'
    fstem = 'japan_ERA5land.'
    
    #Extracting time series
    # Picking location
    lonpick = 136.502
    latpick = 35.667
    #Setting date range to extract
    timeframes = ['1 month', '2 months', '3 months', '4 months', '6 months']
    
    plots = figures(2, figsize = (14, 4))
    
    for y in range(2):
        year = 2017+y
        dstart = datetime.date(year, 6, 1)
        dends = [datetime.date(year, 6, 30), datetime.date(year, 7, 31), 
                   datetime.date(year, 8, 30), datetime.date(year, 9, 30),
                  datetime.date(year, 10, 31)]
        #extracting data
        dayarr, timarr, runoffarr = extract_series(fpath, fstem, lonpick, latpick, 
                                                   dstart, dends[-1], plot = 0)
        for i in range(len(dends)):
            plotRO = 0
            dend = dends[i]
            ndays = (dend-dstart).days+1
            if i == len(dends)-1:
                plotRO = 1
            plots.run_and_plot(runoffarr[:ndays+1], timarr[:ndays+1],timeframes[i],
                               plotRO=plotRO, col = y, zorder = 10-i)
    plots.plot_min_max_flow(D.Wmin, D.Wmax)
    plots.plot_min_max_head(D.Hmin, D.Hmax)
    plots.plot_titles(['2017', '2018'])
    plots.plot_legend()
    # plots.savefigs('Q1')
