# -*- coding: utf-8 -*-
"""
Plotting class
"""
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from damop_model import *
import DefaultParams as D
date_form = DateFormatter("%m-%d")

class figures:
    '''
    class for running and plotting model flow, dam head level, relief flow, 
    runoff and inflow solutions.
    '''
    def __init__(self, ncols, figsize):
        '''
        setting up figures

        Parameters
        ----------
        ncols : int
            number of columns in each figure.
        figsize : (float, float)
            figure size

        '''
        self.fig, self.ax = plt.subplots(1, ncols, figsize = figsize,
                                         sharey = True)
        self.figro, self.axro = plt.subplots(1, ncols, figsize = figsize,
                                             sharey = True)
        if ncols == 2:
            self.axif = [0,0]
            for i in range(ncols):
                self.axif[i] = self.axro[i].twinx()
                if i > 0:
                    self.axif[i].sharey(self.axif[0])
                    
        else:
            self.axif = self.axro.twinx()
        
        self.figgen, self.axgen = plt.subplots(1, ncols, figsize = figsize, 
                                               sharey = True)
        self.figR, self.axR = plt.subplots(1, ncols, figsize = figsize, 
                                           sharey = True)
                
        ylabels = ['Water Head at Dam, $m$', 
                   'Runoff, $m^3$', 'Flowrate, $m^3s^{-1}$',
                   'Relief Flowrate, $m^3s^{-1}$', 'Inflow into Dam, $m^3$']
        self.ax_list = [self.ax, self.axro, self.axgen, self.axR, self.axif]
        
        self.fig_list = [self.fig, self.figro, self.figgen, self.figR]
        self.fignames = ['Waterhead', 'RO-IF', 'Flowrate', 'ReliefFlow']


        if ncols == 1:
            for i in range(len(self.ax_list)):
                a = self.ax_list[i]
                a.xaxis.set_major_formatter(date_form)
                a.set_ylabel(ylabels[i])
                a.set_xlabel('Time')
        else:
            for j in range(ncols):
                for i in range(len(self.ax_list)):
                    a = self.ax_list[i][j]
                    a.xaxis.set_major_formatter(date_form)
                    a.set_ylabel(ylabels[i])
                    a.set_xlabel('Time')
                                        
        for figi in self.fig_list:
            figi.tight_layout()
        self.ncols = ncols

    def run_and_plot(self, runoffarr, timarr, plotlabel,zorder = 1,day = D.day, 
                     Acatch = D.Acatch, kappa = D.kappa, Hmax = D.Hmax, 
                     Hmin = D.Hmin, Wmax = D.Wmax, Wmin = D.Wmin, 
                     Rmax = D.Rmax, sigma = D.sigma, plotRO = False, col =None,
                     color = None, alpha = 1):
        '''
        Runs and plots model solutions given runoff data and time array

        Parameters
        ----------
        runoffarr : array
            array of runoff values
        timarr : array
            array of time stamps for plotting
        plotlabel : str
            label to put in legend 
        zorder : int, optional
            zorder for plot. The default is 1.
        day : float, optional
            length of day in seconds. The default is D.day.
        Acatch : float , optional
            catchment area of the dam. The default is D.Acatch.
        kappa : float, optional
            Half the reservoir area. The default is D.kappa.
        Hmax : float, optional
            Maximum safe water head at the dam. The default is D.Hmax.
        Hmin : float, optional
            Mimmum allowed water head at the dam. The default is D.Hmin.
        Wmax : float, optional
            Maximum flow rate through the turbines. The default is D.Wmax.
        Wmin : float, optional
            Minimum flow rate through the turbines. The default is D.Wmin.
        Rmax : float, optional
            Maximum relief flow rate avoiding the turbines. The default is 
            D.Rmax.
        sigma : float, optional
            Enerfy conversion efficiency of the hydro plant. The default is
            D.sigma.
        plotRO : Bool, optional
            Whether or not to plot runoff and inflow. The default is False.
        col : int, optional
            Index of column to plot in. The default is None.
        color : float, optional
            color to plot in. The default is None.
        alpha : float, optional
            Opacity of plot. The default is 1.

        Returns
        -------
        None.

        '''
        #running damop model
        inflow, x, w, r, gout = damop_model(runoffarr, day, Acatch, 
                                            kappa, Hmax, Hmin, Wmax, Wmin, 
                                            Rmax, sigma)
        if np.sum(w<0) > 0:
            #if any flow is negative
            print('Unphysical Flow')
        if self.ncols == 1:
            ax1, ax2, ax3 = self.ax, self.axgen, self.axR
        else:
            ax1, ax2, ax3 = self.ax[col], self.axgen[col], self.axR[col]
            
        
        #plotting solutions
        ax1.plot(timarr, x, label = plotlabel, zorder = zorder, color = color,
                 alpha = alpha)
        ax2.plot(timarr, w, label = plotlabel, zorder = zorder, color = color,
                 alpha = alpha)
        ax3.plot(timarr, r, label = plotlabel, zorder = zorder, color = color,
                 alpha = alpha)
        
        #plots run off only if specified (runoff is the same when parameters 
        #are varied)
        if plotRO:
            if self.ncols == 1:
                ax4, ax5 = self.axif, self.axro
            else:
                ax4, ax5 = self.axif[col], self.axro[col]
            ax4.plot(timarr, inflow, label = 'Inflow', linestyle = '--', 
                     color = 'black')
            #just here for the legend 
            ax5.plot(timarr[0], runoffarr[0], label = 'Inflow', 
                     linestyle = '--', color = 'black')
            ax5.plot(timarr, runoffarr, label = 'Runoff')
            
    def plot_titles(self, titles):
        '''
        Function to plot subplot titles for figures with more than one column
        Parameters
        ----------
        titles : list
            list of titles, one for each column.

        Returns
        -------
        None.

        '''
        #only for 2 subfigures
        for i in self.ax_list:
            i[0].set_title(titles[0])
            i[1].set_title(titles[1])
    def plot_legend(self):
        '''
        Plot legend in all figures

        Returns
        -------
        None.

        '''
        if self.ncols == 1:
            for i in self.ax_list[:-1]:
                left, right = self.axgen.get_xlim()
                i.set_xlim(left, right)
                i.legend()
            ax4, ax5 = self.axif, self.axro


        else:
            for i in self.ax_list[:-1]:
                i[self.ncols-1].legend(loc='upper right',
                                       bbox_to_anchor=(1.4,1))        
                for j in range(self.ncols):
                    left, right = i[j].get_xlim()
                    i[j].set_xlim(left, right)
            ax4, ax5 = self.axif[self.ncols-1], self.axro[self.ncols-1]
            ax5.legend()
                
    def plot_min_max_head(self, minhead, maxhead, color = 'black',
                          labels = True):
        '''
        Plot minimum and maximum water head levels 

        Parameters
        ----------
        minhead : float
            minimum water head level.
        maxhead : float
            maximum water head level.
        color : str, optional
            color to plot in. The default is 'black'.
        labels : Bool, optional
            Whether or not to label minimum and maximum in legend. The default
            is True.

        Returns
        -------
        None.

        '''
        if labels:
            labels = ['Minimum', 'Maximum']
        else:
            labels = [None, None]
        if self.ncols == 1:
            left, right = self.ax.get_xlim()
            self.ax.hlines(minhead, left, right, linestyle = ':',
                           label = labels[0], color = color)
            self.ax.hlines(maxhead, left, right, linestyle = 'dashdot',
                           label = labels[1], color = color)
            self.ax.set_xlim(left, right)
        else:
            for i in range(self.ncols):
                left, right = self.axgen[i].get_xlim()
                self.ax[i].hlines(minhead, left, right, linestyle = ':',
                                  label = labels[0], color = color)
                self.ax[i].hlines(maxhead, left, right, linestyle = 'dashdot',
                                  label = labels[1], color = color)
                self.ax[i].set_xlim(left, right)
                
    def plot_min_max_flow(self, minflow, maxflow, color = 'black',labels=True):
        '''
        Plot minimum and maximum flow rates

        Parameters
        ----------
        minflow: float
            minimum flow rate.
        maxflow: float
            maximum flow rate.
        color : str, optional
            color to plot in. The default is 'black'.
        labels : Bool, optional
            Whether or not to label minimum and maximum in legend. The default
            is True.

        Returns
        -------
        None.

        '''
        if labels:
            labels = ['Minimum', 'Maximum']
        else:
            labels = [None, None]
        if self.ncols == 1:
            left, right = self.axgen.get_xlim()
            self.axgen.hlines(minflow, left, right, linestyle = ':',
                              label = labels[0], color = color)
            self.axgen.hlines(maxflow, left, right, linestyle = 'dashdot',
                              label = labels[1], color = color)
            self.axgen.hlines(0, left, right, color = 'black')
            self.axgen.set_xlim(left, right)
        else:
            for i in range(self.ncols):
                left, right = self.axgen[i].get_xlim()
                self.axgen[i].hlines(minflow, left, right, linestyle = ':',
                                     label = labels[0], color = color)
                self.axgen[i].hlines(maxflow, left, right,linestyle ='dashdot',
                                     label = labels[1], color = color)
                self.axgen[i].hlines(0, left, right, color = 'black')
                self.axgen[i].set_xlim(left, right)
    def min_max_head_labels(self):
        '''
        Plots singular label for minimum and maximum head in black. 
        Useful if there are multiple minimums and maximums on the same plot

        Returns
        -------
        None.

        '''
        labels = ['Minimum', 'Maximum']
        if self.ncols == 1:
            self.ax.plot([], [], color = 'black', linestyle = ':',
                         label = 'Minimum')
            self.ax.plot([], [], color = 'black', linestyle = 'dashdot',
                     label = 'Maximum')
        else:
            for i in range(self.ncols): 
                self.ax[i].plot([], [], color = 'black', linestyle = ':', 
                                label = 'Minimum')
                self.ax[i].plot([], [], color = 'black', linestyle = 'dashdot',
                                label = 'Maxmimum')
                
    def min_max_flow_labels(self):
        '''
        Plots singular label for minimum and maximum flow in black. 
        Useful if there are multiple minimums and maximums on the same plot

        Returns
        -------
        None.

        '''

        if self.ncols == 1:
            self.axgen.plot([], [], color = 'black', linestyle = ':', 
                            label = 'Minimum')
            self.axgen.plot([], [], color = 'black', linestyle = 'dashdot',
                     label = 'Maximum')
        else:
            for i in range(self.ncols): 
                self.axgen[i].plot([], [], color = 'black', linestyle = ':',
                                   label = 'Minimum')
                self.axgen[i].plot([], [], color = 'black', 
                                   linestyle = 'dashdot', label = 'Maxmimum')
                
                
    def savefigs(self, label):
        for i in range(len(self.fig_list)):
            self.fig_list[i].savefig(f'Plots/{label}_{self.fignames[i]}.pdf', 
                                      bbox_inches = 'tight')
        


