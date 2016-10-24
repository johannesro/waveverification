#!/usr/bin/env python
'''
 usage:
 ./plot_wavedata_timeseries

reads wave verification station data from /lustre/storeA/project/fou/hi/waveverification/data'
and plots time series.

''' 
import scipy as sp
import numpy.ma as ma
import pylab as pl
from netCDF4 import Dataset, date2num, num2date
import os
import datetime as dt

station='draugen' # select station
year = 2016
month = 9
varname = 'Hs' # 'Tm02','FF','Tp', 'DDM' # select variable name
models = ['MWAM4', 'MWAM8', 'ECWAM']
sensor = 0 # specify which wave sensor to use

# advanced user options:
setdates = False # set this swith to True if you want to plot specific days, as set below:w
if setdates:
    t1 = dt.datetime(2016,9,4) # set specific dates for time series plot
    t2 = dt.datetime(2016,9,10)

# set color table for models
ct = {'Subjective': 'b', 'WAM10': 'c', 'WAM4':'m', 'ECWAM':'k', 'LAWAM':'0.25', 'AROME': 'b', 'HIRLAM8': 'y', 'MWAM4':'r', 'EXP':'y', 'MWAM4exp':'w', 'MWAM10':'w', 'MWAM8':'g'}

# settings for time axis annotation
from matplotlib import dates
minorLocator=dates.DayLocator(range(33))
majorLocator=dates.DayLocator(range(5,31,5))
fmt=dates.DateFormatter('%d.%m.%Y') 

# open file
path = '/lustre/storeA/project/fou/hi/waveverification/data'
timep = r'%4.4d%2.2d' % (year, month)
filename = station+'_'+timep+'.nc'
nc       = Dataset(os.path.join(path,filename),mode='r')
time = num2date(nc.variables['time'][:],nc.variables['time'].units)
G = nc.groups # provide each model/observation group as dictionary
OBS  = G['OBS_d22'] # select observation group from dictionary
print ' available models:'
print G.keys()

obs = ma.array(OBS.variables[varname][sensor])
obs.data[obs.mask==True] = sp.nan # make sure all masked values are nan 
obs.mask = sp.logical_or(obs.mask, sp.isnan(obs.data))
units = OBS.variables[varname].units
print 'available variables:'
print OBS.variables.keys()


# read variable from models in netcdf file:
modeldata = {}
for model in models:
    var = G[model].variables[varname][:]
    var[var.mask==True]=sp.nan
    modeldata[model] = var

#
# make time series plot 
#
fig = pl.figure(figsize=[12,6])
ax = fig.add_subplot(111)

if setdates:
    pl.xlim([t1,t2])

ax.xaxis.set_minor_locator(minorLocator)
mask = sp.isfinite(obs)
ax.plot(sp.array(time)[mask], sp.array(obs)[mask], '.', label='obs',lw=1)

for gname, var in modeldata.iteritems():
    ax.plot(time, var[0],'-',color=ct[gname], label=gname, lw=1.5)

ax.legend(fontsize='small')
ax.grid('on',which='minor')
ax.grid('on',which='major',linestyle='--',linewidth=0.5)
ax.set_title(station+' '+varname+' ['+units+']')
fig.tight_layout(pad=0.2)
fig.savefig(station+'_'+varname+'_tseries.png')

pl.show()

