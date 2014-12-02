#!/usr/bin/env python
# 
# collect model data, observations and forecasts at offshore stations and put time series into station files
#
# usage:
# ./collectdata.py [year month [day [numdays]]]
#
# if day is given, the script starts here;
# if day is not given, the script starts checks for which day the verification file (if it already exists) has been last updated
# and continues to fill the file from that day;
# either way the script check of the verification files (for each station) already exists or not; if it doesn't, a new file is created.
# if no argument is given, the script starts or continues to collect data for the current month

import scipy as sp
import numpy as np
import netCDF4 as nc4 #import Dataset, num2date, date2num 
import datetime as dt
import calendar
#import METread
import METread_rpy as METread
import pylab as pl
import os
import sys
#from sys import argv
from stationlist import locations as locations
import collectdatatools as cdt

#
# user parameters
#

collectsubjective = True
collectobservation = True
models = ['WAM10', 'WAM4', 'HIRLAM8', 'AROME']

# 
# get parameters from calling script/user:
#
if len(sys.argv) < 2:
    now = dt.datetime.now()
    startyear, startmonth = now.year, now.month
else:
    startyear, startmonth = int(sys.argv[1]), int(sys.argv[2])

if len(sys.argv) > 3:
    startday = int(sys.argv[3]) # start date for data collection
    updatemode = False
else:
    startday = 1 # This value will again be overwritten from the value in the .nc file, if it exists 
    updatemode = True # Check if there is an nc file that can be updated

starttime = dt.datetime(startyear,startmonth,startday,0)

if len(sys.argv) > 4:
    numdays = int(sys.argv[4])
else:
    numdays = calendar.monthrange(startyear,startmonth)[1] # number of days to be processed
    diff = dt.datetime.today()- dt.datetime(startyear,startmonth,numdays)
    if diff.days < 0:
        if dt.datetime.today().hour < 23: # check how late it is
            numdays = numdays + diff.days - 1 # collect data up until yesterday's forecast
        else:
            numdays = numdays + diff.days # include today's forecasts as well.

print numdays
print updatemode
modlength = 67 # hours in each model run
outpath = '/disk4/waveverification/data'

# errors to be catched during processing:
errlist = IOError, EOFError, KeyError, IndexError


for station, coordinate in locations.iteritems():
    month = startmonth
    filename = starttime.strftime(station+'_%Y%m.nc')

    print(' ')
    print('station '+station)

# Initialize station file object. This creates the file if it doesn't already exists
    stafile = cdt.validationfile(outpath, station, startyear, startmonth)

    stationstartday, stationstarttime, stationnumdays = startday, starttime, numdays # these parameters may be different for each station!
    if updatemode:
        print(int(stafile.nc.last_update_day))
        stationstartday = int(stafile.nc.last_update_day) + 1
        stationstarttime = dt.datetime(startyear,startmonth,stationstartday,0)
        stationnumdays = numdays - stationstartday + 1 # number of days to be processed
    
    print(stationstartday, stationstarttime, stationnumdays)
    print('Collect station data starting from '+str(stationstarttime)+' for '+str(stationnumdays)+' days.')
    print(' ')

    gOBS    = stafile.nc.groups['OBS_d22'].variables # could be stafile.OBS_d22 ??
    if collectsubjective and station in cdt.subjlist:
        #gSubj = stafile.nc.groups['Subjective'].variables
        gSubj = stafile.get_subjective()


#
# collect observations
#
    if collectobservation:
        try:
            print('read d22 files ')
            obs_d22 = METread.obs_d22(station, stationstarttime, numdays=stationnumdays)
            for i,time in enumerate(obs_d22['time']): # do each time step seperately because some steps might be missing
                nci = METread.find_pos1d1(pl.date2num(stafile.time),pl.date2num(time))
                for varname, ncid in gOBS.iteritems():
                    try:
                        starti = (stationstartday-1) * 24
                        ncid[0:obs_d22[varname].shape[0],nci] = obs_d22[varname][:,i]
                    except IndexError:
                        print(varname+' in d22 has wrong dimension')
        except errlist:
            print('error reading and saving D22 tile')
        finally:
            stafile.nc.sync()


#
# Collect subjective analysis
#
    if collectsubjective and station in cdt.subjlist:
        try:
            cdt.collectsubjective(startyear, startmonth, stationstartday, gSubj, station, fstep=6, modlength=67, numdays=stationnumdays)
        except errlist:
            print('subjective analysis not read')
        finally:
            stafile.nc.sync()

#
# collect model data
#
    for mod in models:
        gmod = stafile.get_modelgroup(mod)
        fstep = cdt.reini_dict[mod]
        try:
            cdt.collect(startyear, startmonth, stationstartday, gmod, getattr(METread,mod+'_modrun'), locations[station], fstep=fstep, modlength=67, numdays=stationnumdays)
        except errlist:
            print('model data '+mod+' not processed')
        finally:
            stafile.nc.sync()

    stafile.nc.last_update_day = stationstartday + stationnumdays - 1
    stafile.nc.close()



