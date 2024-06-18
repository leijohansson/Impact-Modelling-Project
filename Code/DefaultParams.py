# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 14:23:36 2024

@author: Linne
"""

Hdam = 161 #dam height #m
A = 13e6 #reservoir area #m^2
kappa = A/2
Acatch = 254e6 #hydrological catchment area #m^2
Hmax = 0.5*Hdam #maximum safe head of water at dam #m
Hmin = 0.2*Hmax #minimum allowed head of water at dam #m
day = 24*60*60 #one day #s
tau = 180*day #Reservoir emptying timescale at max flow rate #s
Wmax = kappa/tau*Hdam #Maximum flow rate through turbines #m^3s^-1
Wmin = 0.1*Wmax #Minimum flow rate through turbines #m^3s^-1
Rmax = 0.2*Wmax #Maximum relief flow avoiding turbines #m^3s^-1
Gmax = 153e6 #Maximum power generation rate by turbines #Watts
sigma = 0.9 #Energy conversion efficiency of hydro plant
mu = Gmax/(sigma*Wmax*Hdam)

lonpick = 136.502
latpick = 35.667
