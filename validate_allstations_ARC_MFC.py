#!/usr/bin/env python
# ./validate 201409
import scipy as sp # scientific function
import numpy as np
import numpy.ma as ma # dealing with masked arrays
import matplotlib
#matplotlib.use('Agg')
import pylab as pl # for plotting
from netCDF4 import Dataset, MFDataset, MFTime, date2num, num2date
import dataanalysis as da
import os # communicate with the OS
import datetime as dt # dealing with time
import validationtools as vt
from stationlist import ARCMFClocations as locations
from stationlist import WMsensors, bestWMsensor
import sys
#import calendar

print("The Python version is %s.%s.%s" % sys.version_info[:3])

interactive=True

timeplist=['201609','201610','201611','201612']

timestr='2016 Sep-Nov'

print('time: '+timestr)


# plotpath
#ppath = '/vol/hindcast3/waveverification/'+timep+'/'
#ppath = '/lustre/storeA/project/fou/hi/waveverification/'+timep+'/'
#ppath = '/disk1/anac/waveverifiction/'
ppath = '/home/johannesro/'

# set color table for models
ct = {'Subjective': 'b', 'WAM10': 'c', 'WAM4':'m', 'ECWAM':'k', 'LAWAM':'0.25', 'AROME': 'b', 'HIRLAM8': 'y', 'MWAM4':'r', 'EXP':'y', 'MWAM4exp':'w', 'MWAM10':'w', 'MWAM8':'g'}

def select_var_from_models(G,varname):
    modeldata={}
    for j, gname in enumerate(G.keys()):
        if gname=='OBS_d22':
            continue
        try:
            var = G[gname].variables[varname][:]
            var[var.mask==True]=sp.nan
            # check if we are dealing with directions and ensure meteorological convention
            if (G[gname].variables[varname].units[0:6] == 'degree'):
                try:
                    if (G[gname].variables[varname].Convention=='oceanographic'):
                        var=var+180
                        var[var>360.]=var[var>360.]-360.
                except AttributeError:
                    var=var
        except KeyError:
            continue
        if sp.isnan(var[0]).all():
            continue
        modeldata.update({gname: var})
    return modeldata


from matplotlib import dates
minorLocator=dates.DayLocator(range(33))
majorLocator=dates.DayLocator(range(5,31,5))
fmt=dates.DateFormatter('%d.%m.%Y') 

varname = 'Hs'
obs_all = []
mod_all = {'ECWAM':[],'MWAM8':[]}

os.system('mkdir -p '+ppath)

for station, parameters in locations.iteritems(): # loop over stations
    obs_long = []
    mod_long = {'ECWAM':[],'MWAM8':[]}


    for timep in timeplist: # loop over months
        print ' '
        print 'read data for station '+station+' for '+timep

        path = '/lustre/storeA/project/fou/hi/waveverification/data'
        filename = station+'_'+timep+'.nc'
        nc       = Dataset(os.path.join(path,filename),mode='r')
        time = num2date(nc.variables['time'][:],nc.variables['time'].units)
        G = nc.groups
        OBS  = G['OBS_d22']

    # Specify which WM sensor to use for validation
        try:
            sensor = bestWMsensor[station]
        except KeyError:
            sensor = 0

        obs = ma.array(OBS.variables[varname][sensor])
        obs.data[obs.mask==True] = sp.nan # make sure all masked values are nan 
        obs.mask = sp.logical_or(obs.mask, sp.isnan(obs.data))
        units = OBS.variables[varname].units

        if (all(sp.isnan(obs.data)) or all(obs.mask==True)):
            print 'no data for '+station+' during '+timestr
            continue

        # select variable from  each model:
        modeldata = select_var_from_models(G,varname)
        #        
        #nc.close()
        #
        obs_long.append(obs)
        for gname, var in modeldata.iteritems():
            if gname in mod_long.keys():
                mod_long[gname].append(var)

    # make arrays from lists
    obs_long = np.concatenate(obs_long)
    mod_longa={}
    for mod in mod_long.keys():
        mod_longa[mod] = np.concatenate(mod_long[mod],axis=1)

    print '%s, %s' % (station, str(mod_long['ECWAM'][1].shape))

    # append to list for all stations:
    obs_all.append(obs_long)
    for gname, var in modeldata.iteritems():
        if gname in mod_all.keys():
            print('append ' +gname+ ' for ' + station)
            mod_all[gname].append(mod_longa[gname])

# make arrays from list
obs = np.array(obs_all)
modeldata = {}
for gname in mod_all.keys():
    modeldata[gname] = np.dstack(mod_all[gname])




#
# make scatter and qq plot
#
fig=pl.figure(figsize=[10,5])
ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)
for gname, var in modeldata.iteritems():
   var2 = np.transpose(var,axes=[2,1,0])
   vt.scqqplot(obs[:,:-24], var2[:,:-24,0],color=ct[gname],  label=gname, ax1=ax1, ax2=ax2)# , prob=sp.arange(0.001,0.999,0.001))
ax1.legend(loc='lower right')
ax1.set_title('analysis '+varname+' ['+units+']'+' '+timestr)
pfilename = 'allstations_'+varname+'_scatterqq_analysis.png'
fig.tight_layout(pad=0.4)
fig.savefig(os.path.join(ppath,pfilename))

#
# make scatter and qq plot for 2-day forecast
#
day2index={'MWAM8':2,'ECWAM':4}
fig=pl.figure(figsize=[10,5])
ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)
for gname, var in modeldata.iteritems():
   var2 = np.transpose(var,axes=[2,1,0])
   vt.scqqplot(obs[:,:-24], var2[:,:-24,day2index[gname]],color=ct[gname],  label=gname, ax1=ax1, ax2=ax2)# , prob=sp.arange(0.001,0.999,0.001))
ax1.legend(loc='lower right')
ax1.set_title('2-day forecast '+varname+' ['+units+']'+' '+timestr)
pfilename = 'allstations_'+varname+'_scatterqq_forecast48h.png'
fig.tight_layout(pad=0.4)
fig.savefig(os.path.join(ppath,pfilename))



#
# plot statistics as function of forcast time
#

fig = pl.figure(figsize=[10,8])
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)
ax1.set_title('allstations '+varname+' forecast skill'+' '+timestr)
for gname, var in modeldata.iteritems():
    print gname
    #vart = np.transpose(var,axes=[0,1,2])
    amerr = vt.forecastskillplot(obs.T, var, G[gname].getncattr('reinitialization_step'), vt.amerr, color=ct[gname],  label=gname, ax=ax1)
    rms = vt.forecastskillplot(obs.T, var, G[gname].getncattr('reinitialization_step'), vt.rmse, color=ct[gname],  label=gname, ax=ax2)
    print 'RMS'
    print rms
    bias = vt.forecastskillplot(-obs.T, -var, G[gname].getncattr('reinitialization_step'), vt.bias, color=ct[gname],  label=gname, ax=ax3)
    corr = vt.forecastskillplot(obs.T, var, G[gname].getncattr('reinitialization_step'), vt.pearsonr, color=ct[gname],  label=gname, ax=ax4)
    print 'correlation'
    print corr
ax1.legend(loc='lower right')
ax1.set_ylabel('MAE ['+units+']')
ax2.set_ylabel('RMSE ['+units+']')
ax3.set_ylabel('model bias ['+units+']')
ax4.set_ylabel('cor. coef.')
ax4.set_xlabel('model lead time [h]')
ax1.set_xlim([0,5*24])
ax1.set_xticks([24,48,72,96,120])
ax2.set_xlim([0,5*24])
ax2.set_xticks([24,48,72,96,120])
ax3.set_xlim([0,5*24])
ax3.set_xticks([24,48,72,96,120])
ax4.set_xlim([0,5*24])
ax4.set_xticks([24,48,72,96,120])

fig.tight_layout(pad=0.2)
pfilename = 'allstations_'+varname+'_forecastskill.png'
fig.savefig(os.path.join(ppath,pfilename))



pl.show()

#
# compute statistics
#


