# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 14:46:24 2024

@author: Linne
"""
import xarray as xr
from DefaultParams import *
import pandas as pd
from readera_runoff_series import subset_field
import numpy as np
from AverageERA5 import logbin_ro

class S2S_one():
    def __init__(self, name, day, month, year, latlon = [latpick, lonpick]):
        filetype = 'PF'
        self.filedir = f'S2SData/{name}/Data/{filetype}/'
        filename = f'japan_{name}_{filetype}.{year}{month}{day}.nc'
        ds = xr.open_dataset(self.filedir + filename)
        if year == '2019':
            ds = ds.sel(time = slice(None, "2019-07-31"))
        self.ds = ds
        #extracting time and runoff from file
        self.ts = ds.time

        #finding closest gridpoint
        alon = ds.longitude.values
        alat = ds.latitude.values
        intlon, intlat = subset_field(alon, alat, lonpick, latpick)
        #extracting runoff at that gridpoint
        ro_acc = ds.ro[:, :, intlat, intlon].values
        self.ro_ts = ro_acc.copy()
        self.ro_ts = np.nan_to_num(self.ro_ts)
        #converting accumulated runoff into runoff
        self.ro_ts[1:] = ro_acc[1:]-ro_acc[:-1]
        #converting units to m
        # ro_bool = self.ro_ts<0
        # print(np.sum(ro_bool))
        self.ro_ts /= 1000
        self.dstart = pd.to_datetime(self.ts.min().values)
        self.dend = pd.to_datetime(self.ts.max().values)
        self.mean_ro = np.mean(self.ro_ts, axis = 1)


    
    def plot_ensemble_mean(self, ax):
        # std_ro = np.std(self.ro_ts, axis =1)#standard error?
        # bot = mean_ro - std_ro
        # bot[np.where(bot<0)] = 0
        bot = np.percentile(self.ro_ts, 10, axis = 1)
        top = np.percentile(self.ro_ts, 90, axis = 1)

        ax.plot(self.ts, self.mean_ro)#, label = 'S2S mean')
        ax.fill_between(self.ts, top, bot, alpha = 0.1) #, label = 'S2S std')
        ax.set_ylabel('runoff, m')
        ax.set_xlabel('Time')
        return ax
    
    def plot_ro_distrib(self, ax, bins):
        flat = self.ro_ts.flatten()
        ax.hist(flat, bins = bins, density = True, histtype = 'step', label = 'S2S', cumulative = True)
    
    def plot_logbin_ro(self, ax, scale = 1.1, multiplier = 1000):
        x, y = logbin_ro(self.ro_ts.flatten(), scale = 1.1, multiplier = 1000)
        ax.scatter(x, y, label = 'S2S')
    
    # def qq_mapping(self):
        # for 
        
    
def read_S2S_cf(name, day, month, year, latlon = [latpick, lonpick]):
    filetype = 'CF'
    filedir = f'S2SData/{name}/Data/{filetype}/'
    filename = f'japan_{name}_{filetype}.{year}{month}{day}.nc'
    ds = xr.open_dataset(filedir + filename)
    if year == '2019':
        ds = ds.sel(time = slice(None, "2019-07-31"))
    #extracting time and runoff from file
    ts = ds.time

    #finding closest gridpoint
    alon = ds.longitude.values
    alat = ds.latitude.values
    intlon, intlat = subset_field(alon, alat, lonpick, latpick)
    #extracting runoff at that gridpoint
    ro_acc = ds.ro[:, intlat, intlon].values
    ro_ts = ro_acc.copy()
    ro_ts = np.nan_to_num(ro_ts)
    #converting accumulated runoff into runoff
    ro_ts[1:] = ro_acc[1:]-ro_acc[:-1]
    #converting units to m
    # ro_bool = self.ro_ts<0
    # print(np.sum(ro_bool))
    ro_ts /= 1000
    dstart = pd.to_datetime(ts.min().values)
    dend = pd.to_datetime(ts.max().values)
    return ro_ts, ts, dstart, dend


#%%
if __name__ == "__main__":
    name = 'ECMF'
    year, month, day = '2017', '06', '01'
    filetype = 'PF'
    S2S = S2S_one(name, day, month, year)
    ndays = len(S2S.ts)
    filelist = os.listdir(S2S.filedir)
    filelist = [i for i in filelist  if year in i]
    years = [i[14:-7] for i in filelist]
    months = [i[18:-5] for i in filelist]
    days = [i[20:-3] for i in filelist]
    fig, ax = plt.subplots()
    
    for i in range(len(filelist)):
        S2S = S2S_one(name, days[i], months[i], years[i])
        S2S.plot_ensemble_mean(ax)
        
