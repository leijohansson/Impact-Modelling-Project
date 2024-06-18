# -*- coding: utf-8 -*-
"""
Upscaling ERA5 data to same resolution as S2S data
"""
from DefaultParams import *
import xarray as xr
from datetime import date
from datetime import timedelta
import datetime
import numpy as np
from readera_runoff_series import *
from logbin import *



class ERA5_averaged():
    '''
    Class to upscale ERA5 data
    '''
    def __init__(self, dstart, dend, latlon = [36, 136.5]):
        '''

        Parameters
        ----------
        dstart : datetime.date
            Start date.
        dend : datetime.date
            End date
        latlon : [lat, lon], optional
            latitude and longitude of grid to upscale.
            The default is [36, 136.5].

        Returns
        -------
        None.

        '''
        fpath = './dataall/'
        fstem = 'japan_ERA5land.'
        fdate = dstart.strftime('%Y%m%d')
        
        #finding end date
        dendp = dend+timedelta(days=1)
        #finding number of days between start and end
        tinterval = dendp-dstart
        ndays = tinterval.days
        
        #opening files and getting grid point indexes
        filename = str(fpath+fstem+fdate+'.nc')
        ds = xr.open_dataset(filename)
        alon = ds.longitude.values
        alat = ds.latitude.values
        self.intlon, self.intlat = subset_field(alon, alat, lonpick, latpick)
        
    
        # dayarr = np.arange(ndays)
        self.timarr = np.arange(np.datetime64(str(dstart)), np.datetime64(str(dendp)))
        self.runoffarr = np.zeros(ndays)

        dcur = dstart
        #find upscaled runoff for each day in range specified
        for n in range(ndays):
            fdate = dcur.strftime("%Y%m%d")
            filename = str(fpath+fstem+fdate+'.nc')
            alon, alat, time, runoff = read_era(filename, 0)
            self.runoffarr[n] = self.find_ave_ro(runoff)
            dcur=dcur+timedelta(days=1)
        
        
        
    def find_ave_ro(self, ro_arr):
        '''
        Average data from surrounding 1.5 x 1.5 area

        Parameters
        ----------
        ro_arr : array
            original era5 runoff array.

        Returns
        -------
        ave_ro : float
            downscaled runoff at the selected grid point

        '''
        #getting 17 by 17 grid centred on grid point
        ro_arr = ro_arr[:, self.intlat-8:self.intlat+9,
                        self.intlon-8:self.intlon+9].values
        
        #only using half of values at the edge
        weight = np.zeros(ro_arr.shape)
        weight = weight*ro_arr
        weight[:, 1:-1, 1:-1] += 0.5
        weight += 0.5
        #setting corners to 0.25
        weight[0, 0] *= 0.5
        weight[0, -1] *= 0.5
        weight[-1, -1] *= 0.5
        weight[-1, 0] *= 0.5
        #not including nans
        weight = np.nan_to_num(weight)
        ro_arr = np.nan_to_num(ro_arr)
        
        ave_ro = np.sum(weight*ro_arr)/np.sum(weight)
        return ave_ro
    
    def plot_ro_distrib(self, ax, nbins):
        '''
        Plots hisogram of runoff values

        Parameters
        ----------
        ax : matplotlib ax object
            ax object to plot on
        nbins : int or array
            number of bines of array of bin edges

        Returns
        -------
        None.

        '''
        ax.hist(self.runoffarr, bins = nbins, density = True, 
                histtype = 'step', label = 'ERA5', cumulative = True)
    def plot_logbin_ro(self, ax, scale = 1, multiplier = 1000):
        '''
        Plot log-binned run off distribution

        Parameters
        ----------
        ax : ax object
            ax to plot distribution on.
        scale : float, optional
            scale for log binning. The default is 1.
        multiplier : float, optional
            amount to multiply runoff by to get integer values for logbinning.
            The default is 1000.

        Returns
        -------
        None.

        '''
        x, y = logbin_ro(self.runoffarr, scale = 1, multiplier = 1000)
        ax.scatter(x, y, label = 'Upscaled ERA5',  marker = 'd', color = 'black')
        


def logbin_ro(runoff, scale=1, multiplier = 1000):
    '''
    Function to logbin runoff values

    Parameters
    ----------
    runoff : array
        runoff values.
    scale : float, optional
        scale for log binning. The default is 1.
    multiplier : float, optional
        Scale to multiply runoff values by before rounding to integers.
        The default is 1000.

    Returns
    -------
    array
        bin centres.
    y : array
        frequency .

    '''
    #converting runoff values to integers
    data = np.round(runoff*multiplier,0).astype(int)
    x, y =logbin(data, scale = scale)
    return x/multiplier,y




        

    
