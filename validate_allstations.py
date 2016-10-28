#!/usr/bin/env python
# ./validate 201409
import scipy as sp
import numpy.ma as ma
import matplotlib
matplotlib.use('Agg')
import pylab as pl
from netCDF4 import Dataset, MFDataset, MFTime, date2num, num2date
import dataanalysis as da
import os
import datetime as dt
import validationtools as vt
from stationlist import bestlocations as locations
from stationlist import WMsensors, bestWMsensor
import sys
#import calendar

print("The Python version is %s.%s.%s" % sys.version_info[:3])

interactive=True

if interactive:
    timep='201609'
    #timep='2013-2014'
else:
    if len(sys.argv) > 1:
        timep = sys.argv[1]
    else:
        now = dt.datetime.now()
        #timep = str(now.year)+str(now.month)
        timep = now.strftime('%Y%m')

if len(timep)==6:
    month = timep[4:6]
    year = timep[0:4]
    timestr = month+'-'+year
    t1 = dt.datetime(int(year), int(month), 1)
    if int(month) < 12:
        t2 = dt.datetime(int(year), int(month)+1, 1)
    else:
        t2 = dt.datetime(int(year)+1, 1, 1) 
else:
    t1 = dt.datetime(int(year), 1, 1)
    t2 = dt.datetime(int(year)+1, 1, 1)
    timestr=timep
print('time: '+timestr)
print(t1, t2)

# plotpath
#ppath = '/vol/hindcast3/waveverification/'+timep+'/'
ppath = '/lustre/storeA/project/fou/hi/waveverification/'+timep+'/'

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
mod_all = {'ECWAM':[], 'MWAM4':[], 'MWAM8':[], 'EXP':[], 'Subjective':[]}

for station, parameters in locations.iteritems():
    print ' '
    print 'read data for station '+station+' for '+timep
#
# open file
    #path = '/vol/hindcast3/waveverification/data'
    path = '/lustre/storeA/project/fou/hi/waveverification/data'
#filename = starttime.strftime(station+'_201312.nc')
    filename = station+'_'+timep+'.nc'
    nc       = Dataset(os.path.join(path,filename),mode='r')
    time = num2date(nc.variables['time'][:],nc.variables['time'].units)
    G = nc.groups
    OBS  = G['OBS_d22']

    os.system('mkdir -p '+ppath+varname)

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

    # append to list for all stations:
    obs_all.append(obs)
    for gname, var in modeldata.iteritems():
        mod_all[gname].append(var)

    nc.close()


# make arrays from list
obs = sp.array(obs_all)
modeldata = {}
for gname in mod_all.keys():
    modeldata[gname] = sp.array(mod_all[gname])

#
# make scatter and qq plot
#
fig=pl.figure(figsize=[10,5])
ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)
for gname, var in modeldata.iteritems():
   vt.scqqplot(obs, var[:,0,:],color=ct[gname],  label=gname, ax1=ax1, ax2=ax2)
ax1.legend(loc='lower right',fontsize='small')
ax1.set_title(station+' '+varname+' ['+units+']'+' obs#'+str(sensor+1)+' '+timestr)
pfilename = 'allstations_'+varname+'_scatterqq.png'
fig.tight_layout(pad=0.2)
fig.savefig(os.path.join(ppath+varname,pfilename))


#
# plot statistics as function of forcast time
#
'''
fig = pl.figure(figsize=[10,8])
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)
ax1.set_title(station+' '+varname+' forecast skill'+' '+timestr)
for gname, var in modeldata.iteritems():
    vt.forecastskillplot(obs, var, G[gname].getncattr('reinitialization_step'), vt.amerr, color=ct[gname],  label=gname, ax=ax1)
    vt.forecastskillplot(obs, var, G[gname].getncattr('reinitialization_step'), vt.rmse, color=ct[gname],  label=gname, ax=ax2)
    vt.forecastskillplot(-obs, -var, G[gname].getncattr('reinitialization_step'), vt.bias, color=ct[gname],  label=gname, ax=ax3)
    vt.forecastskillplot(obs, var, G[gname].getncattr('reinitialization_step'), vt.pearsonr, color=ct[gname],  label=gname, ax=ax4)
ax1.legend(loc='lower right',fontsize='small')
ax1.set_ylabel('MAE ['+units+']')
ax2.set_ylabel('RMSE ['+units+']')
ax3.set_ylabel('model bias ['+units+']')
ax4.set_ylabel('cor. coef.')
ax4.set_xlabel('model lead time [h]')
fig.tight_layout(pad=0.2)
pfilename = station+'_'+varname+'_forecastskill.png'
fig.savefig(os.path.join(ppath+varname,pfilename))
'''

if interactive:
    pl.show()

#
# compute statistics
#


