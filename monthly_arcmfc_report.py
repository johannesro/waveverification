#!/usr/bin/env python
# ./validate 201409
import scipy as sp
import numpy as np
import pylab as pl
import netCDF4
import dataanalysis as da
import os
import datetime as dt
import validationtools as vt
from stationlist import ARCMFClocations as locations
from stationlist import WMsensors, bestWMsensor
import sys
from collectdatatools import validationfile
from dateutil.relativedelta import relativedelta

"""
This script establishes the arcmfc-report file. 
For securing results two harddisks store{A,B} are mirrored and mirrored 
with the script ... located in my home folder. ... is executed each time
the arcmfc-report is created.
"""
print("The Python version is %s.%s.%s" % sys.version_info[:3])

#timeplist=['201707','201708','201709'] # Usage still possible with limitations, see in the following few lines the issues with copying files and naming of folders. This will soon be depricated
now=dt.datetime.now()
pm=now - relativedelta(months=1) # converts to previous month
if pm.month<10:
    timeplist=[str(pm.year)+ '0' +str(pm.month)]
else:
    timeplist=[str(pm.year) +str(pm.month)]

timestr='2018 ' + pm.strftime("%b")
timestrt='2018_' + pm.strftime("%b")
print('time: '+timestr)

# copy necessary files to my folder (only interim solution)
# cp *_201711.nc /lustre/storeA/users/patrikb/waveverification/data/.
cmd = 'cp /lustre/storeA/project/fou/hi/waveverification/data/*_' + timeplist[0] + '.nc /lustre/storeA/users/patrikb/waveverification/data/.'
os.system(cmd)
# plotpath
#ppath = '/lustre/storeA/project/fou/hi/waveverification/Arc-MFC/'+timestrt+'/'
ppath = '/lustre/storeA/project/fou/om/waveverification/Arc-MFC/monthly/'+timestrt+'/'

# set color table for models
ct = {'Subjective': 'b', 'WAM10': 'c', 'WAM4':'m', 'ECWAM':'k', 'LAWAM':'0.25', 'AROME': 'b', 'HIRLAM8': 'y', 'MWAM4':'r', 'EXP':'y', 'MWAM4exp':'w', 'MWAM10':'w', 'MWAM8':'g'}

def select_var_from_models(vf,varname):
    modeldata={}
    for gname in vf.models:
        try:
            mod = vf.get_modelgroup(gname,create=False)
            varraw = mod[varname][:]
            try:
                var = varraw.data
                var[varraw.mask==True]=sp.nan
            except AttributeError:
                var = varraw
            #
            # fill data gap in second forecast range
            #
            var[1] = 0.5*(var[0]+var[2])
            #
            # check if we are dealing with directions and ensure meteorological convention
            #if (G[gname].variables[varname].units[0:6] == 'degree'):
            #    try:
            #        if (G[gname].variables[varname].Convention=='oceanographic'):
            #            var=var+180
            #            var[var>360.]=var[var>360.]-360.
            #    except AttributeError:
            #        var=var
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
mod_all = {'MWAM8':[]}
time_all = []

os.system('mkdir -p '+ppath)

for station, parameters in locations.iteritems():
    print(' ')
    print('verification of station '+station+' for '+timestr)
   

    obs_long = []
    mod_long = {'MWAM8':[]}


    for timep in timeplist: # loop over months
        print ' '
        print 'read data for station '+station+' for '+timep

        # open file
        #path = '/lustre/storeA/project/fou/hi/waveverification/data'
        path = '/lustre/storeA/users/patrikb/waveverification/data'
        
        year,month = int(timep[0:4]),int(timep[4:6])
        vf = validationfile(path,station,year,month)
        time = vf.time
        OBS  = vf.get_obs()
        if station=='draugen':
            time_all = time_all + list(time)



    # Specify which WM sensor to use for validation
        try:
            sensor = bestWMsensor[station]
        except KeyError:
            sensor = 0

        obsraw = OBS[varname][sensor]
        try: 
            obs = obsraw.data 
            obs[obsraw.mask==True] = sp.nan # make sure all masked values are nan 
        except AttributeError:
            obs = obsraw
        units = vf.nc.variables[varname+'_OBS'].units
     
        if all(sp.isnan(obs.data)):
            print('no data for '+station+' during '+timestr)
            continue

        # select variable from  each model:
        modeldata = select_var_from_models(vf,varname)
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

    #print '%s, %s' % (station, str(mod_long['MWAM8'][1].shape))

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


var = modeldata['MWAM8'].transpose(2,1,0)

print('var.shape ', var.shape)
print('obs.shape ', obs.shape)

#
# compute statistics
#

# reshape time axis to join hours
nhours = 6

obsS = obs.reshape(obs.shape[0],nhours,obs.shape[1]/nhours)
varS = var.reshape(var.shape[0],nhours,var.shape[1]/nhours,var.shape[2])

time = time_all[::6]
timeunit = "days since 2001-01-01 12:00:00 UTC"
time_start = time[0].strftime('%Y%m%d')
time_end = time[-1].strftime('%Y%m%d')
print(time_start, time_end)

# produce netcdf file:
nc = netCDF4.Dataset(os.path.join(ppath,'product_quality_stats_ARCTIC_ANALYSIS_FORECAST_WAV_002_006_'+time_start+'-'+time_end+'.nc'),'w')
nc.contact = 'patrikb@met.no'
nc.product = 'Arctic wave model WAM'
nc.production_centre = 'Arctic MFC'
nc.production_unit = 'Norwegian Meteorological Institute'
nc.creation_date = str(dt.datetime.now())
nc.thredds_web_site = 'http://thredds.met.no/thredds/myocean/ARC-MFC/mywave-arctic.html'

ncdims = {'string_length':28, 'areas':1, 'metrics':4, 'surface':1, 'forecasts':10} # and time, unlim
metric_names = [name.ljust(28) for name in ["mean of product", "mean of reference", "mean square difference", "number of data values"]]


nc.createDimension('time',size=None)
for name,dim in ncdims.iteritems():
    nc.createDimension(name,size=dim)

nc_time = nc.createVariable('time','f8',dimensions=('time'))
nc_time[:] = netCDF4.date2num(time,units=timeunit)
nc_time.units = timeunit
nc_time.long_name = 'validity time'

nc_metricnames = nc.createVariable('metric_names','S1', dimensions=(u'metrics',u'string_length')) 
nc_metricnames[:] = netCDF4.stringtochar(np.array(metric_names))

nc_areanames = nc.createVariable('area_names','S1', dimensions=(u'areas',u'string_length')) 
nc_areanames[0] = netCDF4.stringtochar(np.array('North Sea and Norwegian Sea'.ljust(28)))
nc_areanames.long_name = 'area names'
nc_areanames.description = 'region over which statistics are aggregated'

nc_leadtime = nc.createVariable('forecasts','f4',dimensions=('forecasts'))
nc_leadtime.long_name = 'forecast lead time'
nc_leadtime.units = 'hours'
nc_leadtime[:] = sp.arange(12,229,24)


varArcMFCname = {'Hs': 'stats_VHM0'}
varstandardname = {'Hs':'sea_surface_wave_significant_height'}

ncvar = nc.createVariable(varArcMFCname[varname],'f4', dimensions=('time', 'forecasts', 'surface', 'metrics', 'areas'), fill_value=9999.)
ncvar.standard_name = varstandardname[varname]
ncvar.parameter = varArcMFCname[varname]
ncvar.units = 'm'
ncvar.reference_source = 'wave data from offshore platforms available from d22 files at the Norwegian Meteorological Institute'


# calculate statistics (bias and root-mean-square difference)

for leadtime in range(10):
    # mean of product
    ncvar[:,leadtime,0,0,0] = sp.array([sp.mean( vt.returnclean(obsS[:,:,i],varS[:,:,i,leadtime])[1]) for i in range(obsS.shape[2])])
    # mean of reference
    ncvar[:,leadtime,0,1,0] = sp.array([sp.mean( vt.returnclean(obsS[:,:,i],varS[:,:,i,leadtime])[0]) for i in range(obsS.shape[2])])
    # mean square difference
    ncvar[:,leadtime,0,2,0] = sp.array([vt.msd( *vt.returnclean(obsS[:,:,i],varS[:,:,i,leadtime])) for i in range(obsS.shape[2])])
    # number of data values
    ncvar[:,leadtime,0,3,0] = sp.array([ len(vt.returnclean(obsS[:,:,i],varS[:,:,i,leadtime])[1]) for i in range(obsS.shape[2])])

nc.close()







