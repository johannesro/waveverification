#!/usr/bin/env python
'''
 usage:
 ./plot_wavedata_timeseries.py

reads wave verification station data from /lustre/storeA/project/fou/hi/waveverification/data'
and plots time series.

''' 
import scipy as sp
import numpy.ma as ma
import pylab as pl
from netCDF4 import Dataset, date2num, num2date
import os
import datetime as dt
import validationtools as vt

station='heidrun' # select station
year = 2016
month = 10
varname = 'Hs' # 'Tm02','FF','Tp', 'DDM' # select variable name
#models = ['MWAM4', 'MWAM8', 'ECWAM']
models = ['MWAM8']
sensor = 0 # specify which wave sensor to use

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
# make scatter and qq plot
#
fig=pl.figure(figsize=[10,5])
ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)
for gname, var in modeldata.iteritems():
   vt.scqqplot(obs, var[0],color=ct[gname],  label=gname, ax1=ax1, ax2=ax2)
ax1.legend(loc='lower right')
ax1.set_title(station+' '+varname+' ['+units+'] %4d-%2d' %(year,month))
fig.tight_layout(pad=0.6)
fig.savefig(station+'_'+varname+'_scatterqq_%4d%2d.png' %(year,month))

#
# plot statistics as function of forcast time
#
fig = pl.figure(figsize=[10,8])
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)
ax1.set_title(station+' '+varname+' forecast skill'+' %4d-%2d' %(year,month))
for gname, var in modeldata.iteritems():
    vt.forecastskillplot(obs, var, G[gname].getncattr('reinitialization_step'), vt.amerr, color=ct[gname],  label=gname, ax=ax1)
    vt.forecastskillplot(obs, var, G[gname].getncattr('reinitialization_step'), vt.rmse, color=ct[gname],  label=gname, ax=ax2)
    vt.forecastskillplot(-obs, -var, G[gname].getncattr('reinitialization_step'), vt.bias, color=ct[gname],  label=gname, ax=ax3)
    vt.forecastskillplot(obs, var, G[gname].getncattr('reinitialization_step'), vt.pearsonr, color=ct[gname],  label=gname, ax=ax4)
ax1.legend(loc='lower right')
ax1.set_ylabel('MAE ['+units+']')
ax2.set_ylabel('RMSE ['+units+']')
ax3.set_ylabel('model bias ['+units+']')
ax4.set_ylabel('cor. coef.')
ax4.set_xlabel('model lead time [h]')
fig.tight_layout(pad=0.6)
fig.savefig( station+'_'+varname+'_forecastskill_%4d%2d.png' %(year,month))



pl.show()

