# -*- coding: utf-8 -*-
"""
Plots for question 2
"""
#Question 2
from Funcs import figures
import numpy as np
import matplotlib.pyplot as plt
from readera_runoff_series import *
import DefaultParams as D
from damop_model import *
from matplotlib.dates import DateFormatter

#plots for question
if __name__ == '__main__':
    palette = sns.crayon_palette(['Brick Red', 'Green', 'Vivid Violet'])
    
    fpath = './dataall/'
    fstem = 'japan_ERA5land.'
    
    plotHmax = figures(1, (9, 5))
    plotHmin = figures(1, (9, 5))
    plotWmax= figures(1, (9, 5))
    
    year = 2017
    dstart = datetime.date(year, 6, 1)
    dend = datetime.date(year, 9, 30)
    
    dayarr, timarr, runoffarr = extract_series(fpath,fstem,D.lonpick,D.latpick, 
                                               dstart, dend, plot = 0)
    
    Hmaxs = np.array([D.Hmax*0.8, D.Hmax, D.Hmax*1.2])
    
    for i in range(len(Hmaxs)):
        color = palette[i]
        plotRO = 0
        ndays = (dend-dstart).days+1
        if i == len(Hmaxs)-1:
            plotRO = 1
        plotHmax.run_and_plot(runoffarr,timarr,
                              '$H_{max} = $'+str(round(Hmaxs[i],1)), 
                              Hmax = Hmaxs[i], plotRO=plotRO, zorder = 10-i, 
                              col = p, color = color)
        plotHmax.plot_min_max_head(D.Hmin, Hmaxs[i], color = color,
                                   labels = False)
    plotHmax.min_max_head_labels()
    plotHmax.plot_min_max_flow(D.Wmin, D.Wmax)
    plotHmax.plot_legend()
    
       
    Hmins = np.array([D.Hmin*0.8, D.Hmin, D.Hmin*1.2])
    for i in range(len(Hmins)):
        color = palette[i]
        plotRO = 0
        ndays = (dend-dstart).days+1
        if i == len(Hmins)-1:
            plotRO = 1
        plotHmin.run_and_plot(runoffarr, timarr, 
                              '$H_{min} = $'+str(round(Hmins[i],1)), 
                              Hmin = Hmins[i], plotRO=plotRO, zorder = 10-i, 
                              col = p, color = color)
        plotHmin.plot_min_max_head(Hmins[i], D.Hmax, color = color, 
                                   labels = False)
    plotHmin.min_max_head_labels()
    plotHmin.plot_min_max_flow(D.Wmin, D.Wmax)
    plotHmin.plot_legend()
    
    Wmaxs = np.array([0.8, 1, 1.2])*D.Wmax
    Wmins = 0.1*Wmaxs #Minimum flow rate through turbines #m^3s^-1
    Rmaxs = 0.2*Wmaxs #Maximum relief flow avoiding turbines #m^3s^-1
    
    for i in range(len(Wmaxs)):
        color = palette[i]
        plotRO = 0
        ndays = (dend-dstart).days+1
        if i == len(Wmaxs)-1:
            plotRO = 1
        plotWmax.run_and_plot(runoffarr, timarr, 
                              '$W_{max} = $'+str(round(Wmaxs[i], 1)), 
                              Wmax = Wmaxs[i], Wmin = Wmins[i], Rmax =Rmaxs[i],
                              plotRO=plotRO, zorder = 10-i, color = color,
                              col = p)
        plotWmax.plot_min_max_flow(Wmins[i], Wmaxs[i], color = color,
                                   labels = False)
    plotWmax.plot_min_max_head(D.Hmin, D.Hmax)
    plotWmax.min_max_flow_labels()
    plotWmax.plot_legend()
    
    # plotHmin.savefigs(f'Q2_Hmin_{year}')
    # 1plotHmax.savefigs(f'Q2_Hmax_{year}')
    # plotWmax.savefigs(f'Q2_Wmax_{year}')
    
    # plt.close('all')    
