import scipy as sp

import pylab as pl
from netCDF4 import Dataset, MFDataset, MFTime, date2num, num2date
import dataanalysis as da
import os
import datetime as dt
from scipy import stats


def select_var_from_models(varname):
    modeldata={}
    for j, gname in enumerate(nc.groups.keys()):
        if gname=='OBS_d22':
            continue
        try:
            var = G[gname].variables[varname][:]
            var[var.mask==True]=sp.nan
        except KeyError:
            continue
        if sp.isnan(var[0]).all():
            continue
        modeldata.update({gname: var})
    return modeldata



def finite(v):
    return v[sp.isfinite(v)]

def scqqplot(obs, var, label=' ',color='k', ax1=None, ax2=None, prob=sp.arange(0.01, 1., 0.01)):
    if (ax1 == None or ax2 == None):
        fig=pl.figure(figsize=[10,5])
        ax1=fig.add_subplot(121)
        ax2=fig.add_subplot(122)
    #obsc,varc = returnclean(obs,var) # only use times when both obs and model data are available
    if len(obs)==0: # exit function of no common data is available
        return

    maxval = max(list(obs)+list(var))
    minval = min(list(obs)+list(var))
    ax1.plot(obs, var, '.', color=color, label=label+' N='+str(len(obs)))
    ax1.plot([minval,maxval],[minval,maxval],'-r')
    ax1.set_ylabel('model')
    ax1.set_xlabel('observation')
#    ax1.text(0.05,0.9,'N='+str(len(obs)),transform=ax1.transAxes)
    ax1.axis('equal')
    try:
        ax1.legend(loc='lower right',fontsize='small')          
    except TypeError:
        ax1.legend(loc='lower right')          
    obsq = sp.stats.mstats.mquantiles(finite(obs),prob=prob)
    varq = sp.stats.mstats.mquantiles(finite(var),prob=prob)
    ax2.plot(obsq, varq, 'x',ms=5,mew=2, color=color, label=label)
    ax2.plot([minval,maxval],[minval,maxval],'-r')
    ax2.set_ylabel('model')
    ax2.set_xlabel('observation')
    ax2.set_title(str(prob[1]-prob[0])+'-step quantiles')
    try:
        ax2.legend(loc='lower right',fontsize='small')
    except TypeError:
        ax2.legend(loc='lower right')
    ax1.axis([minval,maxval,minval,maxval])
    ax1.axis([minval,maxval,minval,maxval])
    ax1.grid('on')
    ax2.grid('on')
    ax1.axis('scaled')
    ax1.axis('scaled')

#def tsplot(obs, var, label=' ', ax=None):
#    if ax=None:
#       fig=figure()
#    pl.plot(time, obs, label=label,lw=1)

def forecastskillplot(obs, var, reinitializationstep, statfunc, color='k', label=' ', ax=None):
    leadtime=sp.arange(var.shape[0]) * reinitializationstep
    stat = [ statfunc(*returnclean(obs,var[i])) for i in  range(var.shape[0]) ]
    ax.plot(leadtime, stat, color=color, label=label,lw=2)
    ax.set_xlabel('model lead time')
    return leadtime, stat
     

def rmsd(a,b):
    '''
    root mean square deviation
    '''
    a,b = sp.array(a),sp.array(b)
    n = len(a)
    diff2 = (a-b)**2
    return sp.sqrt(diff2.sum()/n)

rmse = rmsd # root mean square error

def msd(a,b):
    '''
    mean square deviation
    '''
    a,b = sp.array(a),sp.array(b)
    n = len(a)
    diff2 = (a-b)**2
    return diff2.sum()/n


def amerr(a,b):
    '''
    absolute mean error
    '''
    a,b = sp.array(a),sp.array(b)
    n = len(a)
    diff = abs(a-b)
    return diff.sum()/n

def bias(a,b):
    '''
    bias
    '''
    a,b = sp.array(a),sp.array(b)
    mask = sp.logical_and(sp.isfinite(a),sp.isfinite(b))
    a, b = a[mask], b[mask]
    return a.mean()-b.mean()

def pearsonr(a,b):
    a,b = sp.array(a),sp.array(b)
    mask = sp.logical_and(sp.isfinite(a),sp.isfinite(b))
    a, b = a[mask], b[mask]
    return stats.pearsonr(a,b)[0]

def returnclean(a,b):
    a,b = sp.array(a),sp.array(b)
    mask = sp.logical_and(sp.isfinite(a),sp.isfinite(b))
    a, b = a[mask], b[mask]
    return a,b

def correlate_vectors(u1,v1, u2,v2):
    ''' 
    u1,v1 : series of vector 1 components
    u2,v2 : series of vector 2 components
    
    returns:
     1 if vectors are always the same
     0 if they are uncorrelated
    -1 if they are opposite'''
    u1 = sp.array(u1)
    u2 = sp.array(u2)
    v1 = sp.array(v1)
    v2 = sp.array(v2)
    upper = u1*u2 + v1*v2
    lower = u1**2 + v1**2 + u2**2 + v2**2
    r =  2*upper.sum() / lower.sum()
    return r

def DD_FF(u,v):
    ''' calculates wind/current speed and direction from u and v components
    #
    if u and v are easterly and northerly components,
    DD is heading direction of the wind. 
    to get meteorological-standard, call DD_FF(-u, -v)
    '''
    DD = ma.arctan2(u, v)*180/sp.pi
    DD[DD < 0] = 360 + DD[DD < 0]
    FF = ma.sqrt(u**2 + v**2)
    return DD, FF

def UV(DD,FF,met=False):
    '''
    #requires testing
    DD is the heading direction of the wind
    to get met. standart, call UV(DD,FF,met=True)
    '''
    u = FF * sp.sin(DD*sp.pi/180)
    v = FF * sp.cos(DD*sp.pi/180)
    if met==True:
        u,v=-u,-v
    return u,v



