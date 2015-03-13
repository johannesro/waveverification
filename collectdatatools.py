import scipy as sp

import netCDF4 as nc4 #import Dataset, MFDataset, MFTime, date2num, num2date
import datetime as dt
from scipy import stats
import numpy as np
import calendar
import METread
import pylab as pl
import os
from stationlist import testlocations as locations
import collectdatatools as cdt
import variables

modlength = 67 # hours in each model run
nleadtimes6 = int(sp.ceil(1.0*modlength/6))
nleadtimes12 = int(sp.ceil(1.0*modlength/12))
nleadtimes24 = int(sp.ceil(1.0*modlength/24))
subjlist = ['ekofiskL', 'ekofisk', 'draugen', 'valhall', 'oseberg', 'osebergc']
subj_reinitime = {'ekofiskL':6, 'ekofisk':6, 'draugen':12, 'valhall':12, 'oseberg':12, 'osebergc':12}
reini_dict = {'WAM4': 12, 'WAM10':12, 'AROME':6, 'HIRLAM8':6, 'LAWAM':12, 'ECWAM':12, 'MWAM4':12, 'MWAM10':12, 'WAMAROME2W':24, 'WAMAROME1W':24}

class validationfile():
    '''

    '''
    def __init__(self, path, station, year, month):
        self.starttime = dt.datetime(year,month,1)
        filename  = self.starttime.strftime(station+'_%Y%m.nc')
        self.filename = os.path.join(path,filename)
        self.station = station
        self.year = year
        self.month = month
        os.system('mkdir -p '+path)
        if not os.path.exists(self.filename):
            print('create new netcdf station file: '+self.filename)
            self.create_file()
        else:
            print('open existing netcdf station file: '+os.path.join(path,filename))
            self.nc = nc4.Dataset(self.filename,mode='a',format='NETCDF4')
        self.time = nc4.num2date(self.nc.variables['time'], self.nc.variables['time'].units)
#
    def create_file(self):
        nc      = nc4.Dataset(self.filename,mode='w',format='NETCDF4')
        self.nc = nc
        self.nc.title = 'wave verification data for '+ self.station
        self.nc.station_name = self.station
        self.nc.last_update_day = 0

        ndays = calendar.monthrange(self.year,self.month)[1]
        dimtime =     nc.createDimension('time', size=24*ndays) #size = None?
        dimleadtime6  = nc.createDimension('lead_time6h', size=nleadtimes6)
        dimleadtime12 = nc.createDimension('lead_time12h', size=nleadtimes12)
        dimleadtime24 = nc.createDimension('lead_time24h', size=nleadtimes24)
        dimwaveobs =  nc.createDimension('n_waveobs', size=3)
        dimwindobs =  nc.createDimension('n_windobs', size=6)

        nctime     = nc.createVariable('time',sp.float64, dimensions=('time',))
        ncleadtime6 = nc.createVariable('lead_time6h',sp.int32, dimensions=('lead_time6h',))
        ncleadtime12 = nc.createVariable('lead_time12h',sp.int32, dimensions=('lead_time12h',))
        ncleadtime24 = nc.createVariable('lead_time24h',sp.int32, dimensions=('lead_time24h',))
# add time to nc file
        ncleadtime6[:] = range(0,modlength+1,6)
        ncleadtime12[:] = range(0,modlength+1,12)
        ncleadtime24[:] = range(0,modlength+1,24)

# generate time for netcdf file
        nctime.units = 'seconds since 1970-01-01 00:00:00'
        ncleadtime6.units = 'hours'
        ncleadtime12.units = 'hours'
        ncleadtime24.units = 'hours'
        t0 = self.starttime - dt.datetime(1970,1,1)
        nctime[:] = t0.total_seconds() + sp.array(range(ndays*24))*3600.

# initialize d22 and subjective variable group
        self.init_d22group()
        if self.station in subjlist:
            self.init_subjective()
        self.nc.sync()
#        
    def get_modelgroup(self, gname):
        '''
        Check if the requested model group exists, if not create it. 
        Return the model group as dictionary of netcdf variables
        '''
        try:
            group = self.nc.groups[gname].variables
        except KeyError:
            self.init_modelgroup(gname)
            group = self.nc.groups[gname].variables
        return group
#
    def init_modelgroup(self,gname):
        setattr(self, gname, {})
        ncgroup = self.nc.createGroup(gname)
        ncgroup.reinitialization_step = reini_dict[gname]
        varlist = getattr(variables, gname)
        for varname,specs in varlist.iteritems():
            lt = str(reini_dict[gname])
            ncvar = ncgroup.createVariable(varname,sp.float32, dimensions=('lead_time'+lt+'h','time'))
            getattr(self, gname)[varname]=ncvar
            ncvar.long_name = specs[0]
            ncvar.units = specs[1]
            if specs[1]=='degree':
                ncvar.Convention = specs[3]
                ncvar.NorthDegree = '0.'
                ncvar.WestDegree = '90.'

    def init_d22group(self):
#
# create variable group for Observations from d22 files
#
        ncOBS      = self.nc.createGroup('OBS_d22')
        self.gOBS = {}
        for varname,specs in variables.d22.iteritems():
            if (varname=='FF') or (varname=='DD'):
                self.gOBS.update({  varname: ncOBS.createVariable(varname,sp.float32, dimensions=('n_windobs','time'))})
            else:
                self.gOBS.update({  varname: ncOBS.createVariable(varname,sp.float32, dimensions=('n_waveobs','time'))})
# variable attributes
            self.gOBS[varname].long_name = specs[0]
            self.gOBS[varname].units = specs[1]
        for varname in ['DD','DDP','DDM']:
            self.gOBS[varname].Convention = 'meteorological'
            self.gOBS[varname].NorthDegree = '0.'
            self.gOBS[varname].WestDegree = '90.'

# 
# create variable group for subjective analysis
#
    def init_subjective(self):
            ncObj = self.nc.createGroup('Subjective')
            ncObj.reinitialization_step = 6
            gObj = {'Hs': ncObj.createVariable('Hs', sp.float32, dimensions=('lead_time6h','time')),
                   'Tp': ncObj.createVariable('Tp', sp.float32, dimensions=('lead_time6h','time')),
                   'FF': ncObj.createVariable('FF', sp.float32, dimensions=('lead_time6h','time'))}
            gObj['Hs'].units = 'm'
            gObj['Tp'].units = 's'
            gObj['FF'].units = 'm/s'
    def get_subjective(self):
        '''
        Check if the subjective  group exists, if not create it. 
        Return the subjective group as dictionary of netcdf variables
        '''
        try:
            group = self.nc.groups['Subjective'].variables
        except KeyError:
            self.init_subjective()
            group = self.nc.groups['Subjective'].variables
        return group
#
# netcdf attributes
# nc standard: variables: units, long_name
# more: mod
# global: Conventions = 'CF', _FillValue, missing_value, title, history



# each variable in the group dictionary is attempted to be read from the specified model reader:
def collect(year, month, sday, ncgroup, modelreader, location, fstep=6, modlength=67, numdays=None):
    '''
    go through each day of the month, including the three previous days and write the data from the model file into the new netcdf file
    the data will be sorted by lead time of the forcast to yield time series of quasi-constant lead times

    year, month: collect during this month
    sday: day of month where to start collection
    ''' 
    forecasttimes = range(0,24,fstep)
    nleadtimes = int(sp.ceil(1.0*modlength/fstep))
    daysofmonth = calendar.monthrange(year,month)[1]
    if numdays == None:
        numdays = daysofmonth
    firstdaytoread = sday-1 if sday>1 else sday-4 # collect data from previous month on the first of each month
    #sday starts at 1, firstdaytoread at 0 for the 1rst of each month
    for day in range(firstdaytoread,sday+numdays-1): # day also start at 0 for the 1rst of each month (python indexing)
        print('day of month: '+str(day+1))
        cyear, cmonth, cday = year, month, day
        if day < 0: 
            cmonth = month-1 
            if month == 1 and day < 0: # go to previous year if we are at the beginning of January
                cyear = year - 1
                cmonth = 12
                cday = 31 + day 
            cday = calendar.monthrange(year,cmonth)[1] + day # get number of days of the previous month and substract current day (add -1 or -2)
        for hour in forecasttimes: # [0, 12] or [0,6,12,18]
            modrun = dt.datetime(cyear,cmonth,cday+1,hour) # get initialization time of model run
            # time series from the model with the given initialization time
            data = modelreader(location, modrun, ncgroup.keys()) # This reader should be universal for each model! (variable list not yet used as parameter)
            # cut the time series into 6 (or 12) hour long pieces and put each piece into the netcdf file according to forecast lead time
            for modstart in range(0, modlength, fstep): # loop over fstep-hour long pieces of model data (Lead Time)
                ncoffset = day*24+hour # time offset of the data in the new netcdf file
                modend   = modstart + fstep
                ncstart   = modstart + ncoffset
                ncend     = modend   + ncoffset
                if not (ncstart < 0 or ncend > daysofmonth*24-1) :
                    for varname, ncid in ncgroup.iteritems():
                        try:
                            ncid[modstart/fstep,ncstart:ncend] = data[varname][modstart:modend]
                        except ValueError:
                            ncend = ncstart + data[varname][modstart:modend].shape[0]
                            ncid[modstart/fstep,ncstart:ncend] = data[varname][modstart:modend]
    return 0


def collectsubjective(year, month, sday, ncgroup, station, fstep=6, modlength=67, numdays=None):
    '''
    go through each day of the month, including the three previous days and write the data from the model file into the new netcdf file
    the data will be sorted by lead time of the forcast to yield time series of quasi-constant lead times

    year, month: collect during this month
    sday: day of month where to start collection
    ''' 
    forecasttimes = range(0,24,fstep)
    nleadtimes = int(sp.ceil(1.0*modlength/fstep))
    daysofmonth = calendar.monthrange(year,month)[1]
    if numdays == None:
        numdays = daysofmonth
    firstdaytoread = sday-1 if sday>3 else sday-4 #sday starts at 1, firstdaytoread at 0 for the 1rst of each month
    for day in range(firstdaytoread,sday+numdays-1): # day also start at 0 for the 1rst of each month (python indexing)
        print('read subjective analysis for day of month: '+str(day+1))
        cyear, cmonth, cday = year, month, day
        if day < 0: 
            cmonth = month-1 
            if month == 1 and day < 0: # go to previous year if we are at the beginning of January
                cyear = year - 1
                cmonth = 12
                cday = 31 + day 
            cday = calendar.monthrange(year,cmonth)[1] + day # get number of days of the previous month and substract current day (add -1 or -2)
        for hour in forecasttimes: # [0, 12] or [0,6,12,18]
            modrun = dt.datetime(cyear,cmonth,cday+1,hour) # get initialization time of model run
            # time series from the model with the given initialization time
#            data = modelreader(location, modrun, ncgroup.keys()) # This reader should be universal for each model! (variable list not yet used as parameter)
            data = METread.subjectiveforecast(modrun, station, ncgroup.keys())
            for varname, ncid in ncgroup.iteritems():
                time,var = data[varname]
                for i in range(len(time)):
                    #print(' ')
                    #print(time[i], modrun)
                    lead_time = time[i] - modrun
                    lead_time = lead_time.total_seconds() /3600.
                    lead_index = int(sp.floor(lead_time /fstep) )
                    timeofmonth = time[i] - dt.datetime(year,month,1)
                    #print timeofmonth
                    timeofmonth = timeofmonth.total_seconds() /3600.
                    #print((timeofmonth<0))
                    #print((timeofmonth> daysofmonth)) 
                    #print(lead_index > int(modlength/fstep))
                    if not ((timeofmonth<0) or (timeofmonth>(daysofmonth*24-1)) or (lead_index > int(modlength/fstep))):
                        #print var[i]
                        #print(lead_index, int(timeofmonth))
                        ncid[lead_index, int(timeofmonth)] = var[i]

    return 0



