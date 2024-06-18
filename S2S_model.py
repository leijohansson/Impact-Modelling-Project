# -*- coding: utf-8 -*-
"""
Plots for question 6
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
from ERAcomparison import scale_runoff

sns.set_context('talk')
name = 'LFPW'
def PositiveColumns(Flows, Inflow, Head, Relief, ROs, Gout):
    '''
    Removes ensemble members for which the model did not converge or
    members with unphysical flow solutions

    Parameters
    ----------
    Flows : array
        Flow array from S2S_one class.
    Inflow : array
        Inflow array from S2S_one class.
    Head : array
        Dam head level array from S2S_one class.
    Relief : array
        relief flow array from S2S_one class.
    ROs : array
        runoff array from S2S_one class.
    Gout : array
        array of Gout for each ensemble.

    Returns
    -------
    Flows_pos : array
        Corrected flow array from S2S_one class.
    Inflow : array
        Corrected Inflow array from S2S_one class.
    Head : array
        Corrected Dam head level array from S2S_one class.
    Relief : array
        Corrected relief flow array from S2S_one class.
    ROs : array
        Corrected runoff array from S2S_one class.
    Gout : array
        corrected Gout array
    nnan : int
        number of ensembles that did not converge.
    nneg : int
        number of ensembles with unphysical flow.

    '''
    #removing nan columns
    nanbool = np.isnan(Flows)
    nanboolsum = np.sum(nanbool, axis = 0)
    nanboolcol = (nanboolsum == 0) 
    Flows_nonan = Flows[:, nanboolcol]
    Inflow = Inflow[:, nanboolcol]
    Relief = Relief[:, nanboolcol]
    Head = Head[:, nanboolcol]
    Gout = Gout[nanboolcol]
    ROs = ROs[:, nanboolcol]
    nnan = Flows.shape[1] - Flows_nonan.shape[1]
    
    #removing columns with negative values
    negbool = (Flows_nonan<0) 
    
    #finding number of negative numbers in each row
    negboolsum = np.sum(negbool, axis = 0) 
    
    #finding rows with no negative numbers
    negboolcol = (negboolsum == 0)
    
    #only keeping ensemble members with no negative flows
    Flows_pos = Flows_nonan[:, negboolcol]
    Inflow = Inflow[:, negboolcol]
    Relief = Relief[:, negboolcol]
    Head = Head[:, negboolcol]
    Gout = Gout[negboolcol]
    ROs = ROs[:, negboolcol]
    nneg = Flows_nonan.shape[1] - Flows_pos.shape[1]

    return Flows_pos, Inflow, Head, Relief, ROs, Gout, nnan, nneg



def S2S_extended(ax, scale_ro = False, n = 0):
    '''
    Runs and plots flow solution for extended S2S forecast

    Parameters
    ----------
    ax : matplotlib ax object
        axis to plot flows on.
    scale_ro : bool, optional
        whether to calibrate runoff to ERA5 data. The default is False.
    n : int, optional
        0 or 1. Determines which forecast start date The default is 0.

    Returns
    -------
    Gout : array
        power generated for each ensemble.
    nnan : int
        number of ensemble solutions that did not converge.
    nneg : int
        number of ensembles that had unphysical flow.

    '''
    #setting dates for the two different scenarios
    if n == 0:
        forecast = S2S_one(name, '28', '06', '2018')
        dstart = datetime.date(2018, 3, 28)
        dend  = datetime.date(2018, 6, 27)
    if n==1:
        forecast = S2S_one(name, '05', '07', '2018')
        dstart = datetime.date(2018, 4, 5)
        dend  = datetime.date(2018, 7, 4)
        
    #reading in downscaled ERA data
    ERA = ERA5_averaged(dstart, dend)
    
    #setting up arrays
    ROs = np.zeros((len(forecast.ts)+len(ERA.timarr), 50))
    Flows = ROs.copy()
    Relief = ROs.copy()
    Inflow = ROs.copy()
    Head = ROs.copy()
    Gout = np.zeros(50)
    
    #looping through all ensemble members
    for i in range(50):
        try:
            #calibrating forecast data
            ro = qq_mapping(forecast.ro_ts[:, i])
            ro = forecast.ro_ts[:, i]
            if scale_ro:
                ro = scale_runoff(ro)
                
            #creating extended forecast using downscaled ERA5 data
            ro = np.concatenate((ERA.runoffarr, ro))
            ts = np.concatenate((ERA.timarr,forecast.ts))
            
            #running damop model
            inflow, x, w, r, gout = damop_model(ro, D.day, D.Acatch, D.kappa,
                                                D.Hmax, D.Hmin, D.Wmax, D.Wmin,
                                                D.Rmax, D.sigma)
            Inflow[:,i] = inflow
            Head[:,i] = x
            Flows[:,i] = w
            Relief[:,i] = r
            ROs[:, i] = ro  
            Gout[i] = gout
        except:
            #unphysical solutions are shown as nans
            Inflow[:,i] = np.nan
            Head[:,i] = np.nan
            Flows[:,i] = np.nan
            Relief[:,i] = np.nan       
            ROs[:, i] = ro      
            Gout[i] = np.nan
    
    #selecting valid solutions
    Flows, Inflow, Head, Relief, ROs, Gout, nnan, nneg \
        = PositiveColumns(Flows, Inflow, Head, Relief, ROs, Gout)
        
    #plotting all valid solutions, mean ensemble solutions and 10-90% range
    ax.plot(ts, Flows, color = 'firebrick', alpha = 0.3)
    ax.plot(ts[0], Flows[0,0], color = 'firebrick', alpha = 0.3, 
            label = 'Ensemble Members')
    ax.plot(ts, Flows.mean(axis = 1), color = 'black', zorder = 5, 
            label = 'Ensemble Mean')
    ax.fill_between(ts, np.percentile(Flows, 10, axis= 1)
                    ,np.percentile(Flows, 90, axis= 1), color = 'black'
                    , alpha = 0.2, ec = None, label = '10%-90%')
    print(str(Flows.shape[1])+' ensemble members had physical results')
    return Gout, nnan, nneg


def only_S2S(ax, scale_ro = False):
    '''
    Runs and plots flow solution for one S2S forecast

    Parameters
    ----------
    ax : matplotlib ax object
        axis to plot flows on.
    scale_ro : bool, optional
        whether to calibrate runoff to ERA5 data. The default is False.

    Returns
    -------
    Gout : array
        power generated for each ensemble.
    nnan : int
        number of ensemble solutions that did not converge.
    nneg : int
        number of ensembles that had unphysical flow.

    '''
    forecast = S2S_one(name, '28', '06', '2018')
    
    ROs = np.zeros((len(forecast.ts), 50))
    Flows = ROs.copy()
    Relief = ROs.copy()
    Inflow = ROs.copy()
    Head = ROs.copy()
    Gout = np.zeros(50)
    ts = forecast.ts
    for i in range(50):
        try:
            ro = qq_mapping(forecast.ro_ts[:, i])
            ro = forecast.ro_ts[:, i]
            if scale_ro:
                ro = scale_runoff(ro)
            inflow, x, w, r, gout = damop_model(ro, D.day, D.Acatch, D.kappa,
                                                D.Hmax, D.Hmin, D.Wmax, D.Wmin,
                                                D.Rmax, D.sigma)
            Inflow[:,i] = inflow
            Head[:,i] = x
            Flows[:,i] = w
            Relief[:,i] = r
            ROs[:, i] = ro  
            Gout[i] = gout
        except:
            Inflow[:,i] = np.nan
            Head[:,i] = np.nan
            Flows[:,i] = np.nan
            Relief[:,i] = np.nan       
            ROs[:, i] = ro      
            Gout[i] = np.nan
    
    Flows, Inflow, Head, Relief, ROs, Gout, nnan, nneg \
        = PositiveColumns(Flows, Inflow, Head, Relief, ROs, Gout)    
    print(str(Flows.shape[1])+' ensemble members had physical results')

    ax.plot(ts, Flows, color = 'blue', alpha = 0.3)
    ax.plot(ts[0], Flows[0,0], color = 'blue', alpha = 0.3,
            label = 'Ensemble Members')
    ax.plot(ts, Flows.mean(axis = 1), color = 'black', zorder = 5,
            label = 'Ensemble Mean')
    ax.fill_between(ts, np.percentile(Flows, 10, axis= 1)
                    ,np.percentile(Flows, 90, axis= 1), color = 'blue'
                    , alpha = 0.3, label = '10%-90%', ec = None, zorder = 6)
    return Gout, nnan, nneg
    
#plotting flow solutions for ERA5, upscaled ERA5, calibrated upscaled ERA5, 
#and the two extended forecasts
if  __name__ == '__main__':
    fig, ax = plt.subplots(1, 1, figsize = (18, 6))
    
    #setting x axis date format
    date_form = DateFormatter("%m-%d")
    ax.xaxis.set_major_formatter(date_form)
    ax.set_ylabel('Flow, $m^3s^{-1}$')
    ax.set_xlabel('Time')
    
    #plotting only S2S solutions
    Gout_daily, nnan, nneg =only_S2S(ax, scale_ro = True)
    ax.legend(loc = 'upper right')
    # fig.savefig('S2S_only.pdf', bbox_inches = 'tight')
    print(f'{nnan} ensembles did not converge, {nneg} had unphysical flow')
    
    #plotting ERA5 solutions and extended forecast solutions
    n_f = 2
    fig, axs = plt.subplots(n_f+1, 1, figsize = (16, (n_f+1)*5), sharex = True,
                            sharey = True)
    for axi in axs:
        date_form = DateFormatter("%m-%d")
        axi.xaxis.set_major_formatter(date_form)
        axi.set_ylabel('Flow, $m^3s^{-1}$')
    axi.set_xlabel('Time')
    
    dstart = datetime.date(2018, 3, 28)
    dend  = datetime.date(2018, 7, 30)
    
    #plotting ERA5, downscaled ERA5 and calibrated downscaled ERA5 solutions
    ERA = ERA5_averaged(dstart, dend)
    inflow, x, w, r, gout = damop_model(ERA.runoffarr, D.day, D.Acatch,D.kappa,
                                        D.Hmax, D.Hmin, D.Wmax, D.Wmin,
                                        D.Rmax, D.sigma)
    axs[0].plot(ERA.timarr, w, color = 'black', label = 'Upscaled ERA5 Data',  
                linestyle = '--', zorder = 5)
    scaleERA5 = scale_runoff(ERA.runoffarr)
    inflow, x, w, r, gout = damop_model(ERA.runoffarr, D.day, D.Acatch,D.kappa,
                                        D.Hmax, D.Hmin, D.Wmax, D.Wmin,
                                        D.Rmax, D.sigma)
    axs[0].plot(ERA.timarr, w, color = 'blue', 
                label = 'Calibrated Upscaled ERA5 Data')
    
    fpath = './dataall/'
    fstem = 'japan_ERA5land.'
    dayarr, timarr, runoffarr=extract_series(fpath, fstem, D.lonpick,D.latpick, 
                                               dstart, dend, plot = 0)
    inflow, x, w, r, gout = damop_model(runoffarr, D.day, D.Acatch, D.kappa,
                                        D.Hmax, D.Hmin, D.Wmax, D.Wmin,
                                        D.Rmax, D.sigma)
    
    axs[0].plot(timarr, w, color = 'red', label = 'ERA5 Data')
    nnans = []
    nnegs = []
    dates = [datetime.date(2018, 6, 28), datetime.date(2018, 7, 5),
             datetime.date(2018, 7, 12)]
    up, down = axs[0].get_ylim()
    axs[0].set_title('ERA5 Data')
    axs[0].legend()
    
    #plotting extended forecast solutions
    for i in range(n_f):
        Gouts_ext, nnan_ext, nneg_ext = S2S_extended(axs[i+1], scale_ro=True, 
                                                     n=i)
        nnans.append(nnan_ext)
        nnegs.append(nneg_ext)
        axs[i+1].vlines(dates[i], down, up, color = 'black',linestyle='dotted', 
                        label = 'Start of Forecast')
        axs[i+1].set_title(f'Calibrated S2S Data: {dates[i]}')
    for axi in axs:
        axi.set_xlim(datetime.date(2018, 3, 28), datetime.date(2018, 7, 30))
    axi.legend()
    # fig.savefig('S2S_extended.pdf', bbox_inches = 'tight')
