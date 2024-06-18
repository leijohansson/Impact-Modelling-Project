"""
Applying an impact model for hydroelectric dam management driven by
a time series of runoff data

Author: 2020, John Methven
Revised to include relief flow and smoother solution
        July 2021, John Methven
"""
#
# Library functions needed to run damop_model()
# Not including standard imports like numpy
#
from scipy import optimize
from scipy import signal
import numpy as np

def damop_model(runoffarr, dt, catcharea, kappa, hmax, hmin, wmax, wmin, rmax, sigma):
    '''
    Implementation of the dam operation model of Hirsch et al (2014).
    Input: 
    :runoffarr  - input time series for runoff data
    :dt         - runoff accumulation interval per record
    :catcharea  - catchment area for the dam
    :kappa      - parameter relating reservoir depth to volume
    :hmax       - maximum water head (constraint on optimization)
    :hmin       - minimum water head
    :wmax       - maximum flow rate through turbines
    :wmin       - minimum flow rate to maintain some power generation
    :rmax       - maximum relief flow rate, bypassing turbines in flood
    :sigma      - operational efficiency of power generation by dam
    Output: 
    :inflow     - input time series for inflow to reservoir  
    :x          - output time series for water head at dam
    :w          - output solution for optimum flow rate through turbines
    :r          - output solution for relief flow rate
    :gout       - value of time integrated generation for optimum solution (MW-days)
    '''       
    # print()
    # print('damop_model has been called with the constraints:')
    # print('wmax = ',wmax,'   wmin = ',wmin,'   hmax = ',hmax,'   hmin = ',hmin)
    #
    # Convert runoff data from units of m to an equivalent inflow in m^3 s^-1
    # Assume that the same runoff rate applies to the entire catchment area for dam
    #
    runoffave = np.mean(runoffarr)
    inflow = catcharea*runoffarr/dt
    n = len(inflow)
    inmax = max(inflow)
    #
    # Set parameter used to control computational mode using filter similar to Robert-Asselin
    # Recommend 0 because filter introduces an offset of W relative to I in optimization.
    #
    alpha = 0.0
    #
    # Apply running mean to the inflow data if required for smoother solution 
    # to the optimisation. Averaging window length = nwin.
    #
    nwin = 3
    inflow = running_mean(inflow, nwin)
    #
    # Scale mu so that the sum of generation over time points is approx one.
    # This gives a better numerical solution in the optimisation for max generation
    # by reducing numerical truncation error in the calculation.
    #
    mu = 1.0/(n*sigma*wmax*hmax)
    #
    # The dam management optimization model is set up in the mathematical form of a 
    # quadratic programming problem.
    # The only input time series is the inflow to the reservoir.
    # The model solves for the water head at the dam maximizing power generation.
    # This then gives the flow rate through the turbines.
    # However, contraints are applied on maximum and minimum water level 
    # and maximum/minimum flow rate through the turbines.
    #
    # The equation for generation can be written in the form
    # 
    # G = 0.5*H^T P H + q^T H
    #
    # where H is the head time series we are solving for (a 1-D array) and 
    # P is a matrix and q is also a 1-D time series (scaled inflow).
    # The notation ^T means the transpose of the matrix. 
    # Quadratic programming aims to minimize -G which is equivalent to max(G).
    #
    q = -mu*sigma*inflow
    umat = np.zeros((n, n))
    inmat = np.zeros((n, n))
    cmat = np.zeros((n, n))
    for i in range(n):
        umat[i, i] = 1
        inmat[i, i] = inflow[i]

    for j in range(n-2):
        i = j+1
        cmat[i, i-1] = -1 + 0.5*alpha
        cmat[i, i]   = -alpha
        cmat[i, i+1] = 1 + 0.5*alpha
    
    pscal = mu*sigma*(kappa/dt)*cmat
    wscal = -0.5*(kappa/dt)*cmat
    #
    # Set constraints on the water head at the dam: hmin <= h <= hmax
    # Optimization requires that constraints actually need to be applied in form:
    # Amat x <= b  (where in this problem the vector x is head time series, h).
    # For Amat x >= b it is necessary to re-arrange to -Amat x <= -b.
    # Therefore to apply hmin <= h <= hmax, the matrix Amat is the unit matrix.
    #
    hscal = umat
    hmaxcons = np.ones(n)*hmax
    hmincons = np.ones(n)*hmin    
    #
    # Set constraints on the flow rate 
    # based on the parameters Wmax, Rmax and Wmin.
    # The form of the contraints means that it must be applied to range of W*h:
    # Wmin*hmin <= W*h <= (Wmax+Rmax)*hmax
    #
    gscal = wscal + inmat
    gmaxcons = np.zeros(n)
    gmincons = np.zeros(n)
    for i in range(n):
        gmaxcons[i] = (wmax+rmax)*hmax
        gmincons[i] = wmin*0.5*(hmin+hmax)
    #
    # Construct a single matrix describing Amat and vector for constraint values b
    # in the form required by optimize.minimize
    #
    vmat = np.concatenate((gscal, -gscal, hscal, -hscal), axis=0)
    vcons = np.concatenate((gmaxcons, -gmincons, hmaxcons, -hmincons))
    
    # print('Now apply quadratic minimization technique')
    
    def gen(x, sign=1.):
        return sign * (0.5*np.dot(x.T, np.dot(pscal, x)) + np.dot(q.T, x))
    
    def jac(x, sign=1.):
        return sign * (np.dot(x.T, pscal) + q.T)
    
    cons = {'type':'ineq',
            'fun':lambda x: vcons - np.dot(vmat, x),
            'jac':lambda x: -vmat}
    
    opt = {'disp':True, 'maxiter':100, 'ftol':1e-08}

    #
    # Obtain solution by minimization nouter times. Smooth the input first guess 
    # and results for head, h, which removes noise and any numerical instability in 
    # optimal solution for the flow rate time series, W.
    # Note that the minimize method does not always find a solution consistent 
    # with the contraints imposed (depending on the first guess data) and these
    # failed attempts are not included in the average solution.
    #
    nouter = 3
    istsuccess = 1
    ic = -1
    afac = 0.5
    xinit = hmax*(afac + 0.1*np.random.randn(n))
    nwin = min([41, 2*round(0.2*n)+1])
    # print('running mean window length, nwin = ',nwin)
    xinit = running_mean(xinit, nwin)
    
    for io in range(nouter):
    #while istsuccess == 1:
        #
        # First guess values for x (water head).
        # Random variation on top of constant level.
        # Smooth to reduce 2-grid noise in input data.
        #
        ic = ic+1
        res_cons = optimize.minimize(gen, xinit, jac=jac, constraints=cons,
                                 method='SLSQP', options=opt)
        xup = res_cons['x']
        fup = res_cons['fun']  
        stexit = res_cons['status']
    
        if stexit != 4:
            if istsuccess == 1:
                x = xup
                x = running_mean(x, nwin)
                xinit = x
                f = fup
                # print('Constrained optimization')
                # print(res_cons)
                # print('iter ',ic,' f = ',f)
                istsuccess = 0
            else:
                if (fup/f) < 2:
                    afac = float(ic+1)/nouter
                    x = afac*x + (1-afac)*xup
                    x = running_mean(x, nwin)
                    xinit = x
                    f = afac*f + (1-afac)*fup
                    # print('iter ',ic,' f = ',f)
        if ic == nouter:
            # print(nouter,' outer iterations finished without reaching result')
            istsuccess = 1
    # end outer loop
    
    #
    # Optimisation returns the head in variable x
    # Total flow rate ft = W+R is calculated from head and known inflow rate
    # Total flow is diverted into relief flow when it exceeds Wmax (and Rmax > 0)
    #
    ft = np.dot(wscal, x) + inflow
    w = np.copy(ft)
    r = np.zeros(n)
    excessflow = np.where(ft > wmax)
    if rmax > 0:
        w[excessflow] = wmax
        r[excessflow] = ft[excessflow]-wmax
    
    gout = -f
    
    return inflow, x, w, r, gout


def running_mean(xarr, nwin):
    '''
    Apply running mean filter through array
    Inputs:
        xarr    - array to filter
        nwin    - number of points in the filter window (odd number expected)
    Output:
        xfilt   - same length as xarr after application of filter
    '''
    n = len(xarr)
    xfilt = np.copy(xarr)
    ist = int(nwin/2)
    xconv = np.convolve(xarr, np.ones(nwin),'valid')/nwin
    nconv = len(xconv)
    xfilt[ist:n-ist] = xconv[:]
    xfilt[0:ist] = xconv[0]
    xfilt[n-ist:n] = xconv[nconv-1]
    
    return xfilt
