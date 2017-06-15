import scipy as sp

import netCDF4 as netCDF4 #import Dataset, MFDataset, MFTime, date2num, num2date
import datetime as dt
from scipy import stats
import numpy as np
import calendar
import METread
import pylab as pl
import os
from stationlist import testlocations as locations
import collectdatatools as cdt
import config

#modlength = 67 # hours in each model run
# these list should be in config.py
subjlist = ['ekofiskL', 'ekofisk', 'draugen', 'valhall', 'oseberg', 'osebergc']
subj_reinitime = {'ekofiskL':6, 'ekofisk':6, 'draugen':12, 'valhall':12, 'oseberg':12, 'osebergc':12}
reini_dict = {'WAM4': 12, 'WAM10':12, 'AROME':6, 'HIRLAM8':6, 'LAWAM':12, 'ECWAM':12, 'MWAM4':6, 'MWAM8':24, 'EXP':6, 'WAMAROME2W':24, 'WAMAROME1W':24}

nleadtimes6 = int(sp.ceil(1.0*67/6))
nleadtimes12 = int(sp.ceil(1.0*241/12))
nleadtimes24 = int(sp.ceil(1.0*241/24))

def nleadtimes(ml, reini):
    return 

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
            self.nc = netCDF4.Dataset(self.filename,mode='a',format='NETCDF3')
            mstr = self.nc.models
            self.models = mstr.strip().split(' ')
        self.time = netCDF4.num2date(self.nc.variables['time'][:], self.nc.variables['time'].units)

    def create_file(self):
        nc      = netCDF4.Dataset(self.filename,mode='w',format='NETCDF4')
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
        ncleadtime6[:] = range(0,67,6)
        ncleadtime12[:] = range(0,241,12)
        ncleadtime24[:] = range(0,241,24)

# generate time for netcdf file
        nctime.units = 'seconds since 1970-01-01 00:00:00'
        ncleadtime6.units = 'hours'
        ncleadtime12.units = 'hours'
        ncleadtime24.units = 'hours'
        t0 = self.starttime - dt.datetime(1970,1,1)
        nctime[:] = t0.total_seconds() + sp.array(range(ndays*24))*3600.
        #self.time = netCDF4.num2date(self.nc.variables['time'][:], self.nc.variables['time'].units)

# list for model names
        nc.models = ''
        self.models = []

# initialize d22 and subjective variable group
        self.init_obs()
        if self.station in subjlist:
            self.init_subjective()
        self.nc.sync()
#        
    def get_obs(self):
        '''
        Return dictionary of observation nc variables
        '''
        gOBS = {}
        variables = self.nc.getncattr('observations').strip().split(' ')
        for var in variables:
            gOBS[var] = self.nc.variables[var + '_OBS']
        return gOBS
#
    def get_modelgroup(self, gname, create=True):
        '''
        Check if the requested model group exists, if not create it. 
        Return the model group as dictionary of netcdf variables
        '''
        if gname in self.models:
            variables = self.nc.getncattr(gname).strip().split(' ')
            group = {}
            for var in variables:
                group[var] = self.nc.variables[var+'_'+gname]
        else:
            if create == True:
                group = self.init_modelgroup(gname)
        setattr(self, gname, group) # set dictionary with nc variable handles as attribute to validationfile object
        return group
#
    def init_modelgroup(self,gname):
        group = {}
        self.nc.setncattr(gname+'_reinitialization_step' , reini_dict[gname])
        self.nc.setncattr(gname , '')
        self.nc.models = self.nc.models+gname+' '
        self.models.append(gname)
        varlist = getattr(config, gname)
        for varname,specs in varlist.iteritems():
            lt = str(reini_dict[gname])
            ncvar = self.nc.createVariable(varname+'_'+gname,sp.float32, dimensions=('lead_time'+lt+'h','time'))
            self.nc.setncattr(gname , self.nc.getncattr(gname)+varname+' ') # append variable to nc attribute 
            group[varname] = ncvar # append nc object to dictionary
            ncvar.standard_name = specs['standard_name']
            ncvar.short_name = specs['short_name']
            ncvar.units = specs['units']
            if specs['units']=='degree':
                ncvar.Convention = specs['convention']
                ncvar.NorthDegree = '0.'
                ncvar.WestDegree = '90.'
        setattr(self, gname, group) # set dictionary with nc variable handles as attribute to validationfile object
        return group

    def init_obs(self):
#
# create variable group for Observations
#
        self.gOBS = {}
        self.nc.observations = ''
        for varname,specs in config.d22.iteritems():
            self.nc.observations = self.nc.observations + varname + ' '  
            if (varname=='FF') or (varname=='DD'):
                self.gOBS.update({  varname: self.nc.createVariable(varname+'_OBS',sp.float32, dimensions=('n_windobs','time'))})
            else:
                self.gOBS.update({  varname: self.nc.createVariable(varname+'_OBS',sp.float32, dimensions=('n_waveobs','time'))})
# variable attributes
            self.gOBS[varname].standard_name = specs['standard_name']
            self.gOBS[varname].short_name = specs['short_name']
            self.gOBS[varname].units = specs['units']
        for varname in ['DD','DDP','DDM']:
            self.gOBS[varname].Convention = 'meteorological'
            self.gOBS[varname].NorthDegree = '0.'
            self.gOBS[varname].WestDegree = '90.'
        return self.gOBS

# 
# create variable group for subjective analysis
#
    def init_subjective(self):
            self.models.append('Subjective')
            self.nc.models = self.nc.models+'Subjective '
            self.nc.setncattr('Subjective_reinitialization_step', 6)
            gObj = {'Hs': self.nc.createVariable('Hs_Subjective', sp.float32, dimensions=('lead_time6h','time')),
                   'Tp': self.nc.createVariable('Tp_Subjective', sp.float32, dimensions=('lead_time6h','time')),
                   'Tm02': self.nc.createVariable('Tm02_Subjective', sp.float32, dimensions=('lead_time6h','time')),
                   'FF': self.nc.createVariable('FF_Subjective', sp.float32, dimensions=('lead_time6h','time')),
                   'DD': self.nc.createVariable('DD_Subjective', sp.float32, dimensions=('lead_time6h','time'))}
            gObj['Hs'].units = 'm'
            gObj['Tp'].units = 's'
            gObj['Tm02'].units = 's'
            gObj['FF'].units = 'm/s'
            gObj['DD'].units = 'degree'
            gObj['DD'].Convention = 'meteorological'
            self.nc.Subjective = 'Hs Tp Tm02 FF DD'
            self.Subjective = gObj
            return gObj

    def get_subjective(self, create=True):
        '''
        Check if the subjective  group exists, if not create it. 
        Return the subjective group as dictionary of netcdf variables
        '''
        if 'Subjective' in self.models:
            group = {}
            for var in self.nc.Subjective.strip().split(' '):
                group[var] = self.nc.variables[var+'_Subjective']
        else:
            if create==True:
                group = self.init_subjective()
        self.Subjective = group
        return group
#
# netcdf attributes
# nc standard: variables: units, short_name
# more: mod
# global: Conventions = 'CF', _FillValue, missing_value, title, history



# each variable in the group dictionary is attempted to be read from the specified model reader:
def collect(year, month, sday, modelgroup, modelreader, location, fstep=6, modlength=67, numdays=None, numdays_previousmonth=3):
    '''
    go through each day of the month, including the three previous days and write the data from the model file into the new netcdf file
    the data will be sorted by lead time of the forcast to yield time series of quasi-constant lead times

    year, month: collect during this month
    sday: day of month where to start collection
    modelgroup: dictionary with nc variable instances
    fstep: number of hours between each model reinitialization
    ''' 
    forecasttimes = range(0,24,fstep)
    nleadtimes = int(sp.ceil(1.0*modlength/fstep))
    daysofmonth = calendar.monthrange(year,month)[1]
    if numdays == None:
        numdays = daysofmonth
    # collect data from the last 5 days of previous month on the first of each month
    firstdaytoread = sday-1 if sday>1 else sday-(numdays_previousmonth+1) 
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
            variables = modelgroup.keys()
            data = modelreader(location, modrun, variables) # This reader should be universal for each model! (variable list not yet used as parameter)
            ensure_hourly(data)
            #print data
    	    #print len(data['time']), len(data['Hs'])
            try:
                modoffset = modrun-data['time'][0] #check at which hour the model output starts. Assume hourly output thereafter!
            except TypeError:
                continue
            mo = int(modoffset.total_seconds()/3600.) # get number of hours
            # cut the time series into 6 (or 12) hour long pieces and put each piece into the netcdf file according to forecast lead time
            for modstart in range(0, modlength, fstep): # loop over fstep-hour long pieces of model data (Lead Time)
                ncoffset = day*24+hour # time offset of the data in the new netcdf file
                modend   = modstart + fstep
                ncstart   = modstart + ncoffset
                ncend     = modend   + ncoffset
                if not (ncstart < 0 or ncend > daysofmonth*24) :
                    for varname, ncid in modelgroup.iteritems():
                        try:
                            ncid[modstart/fstep,ncstart:ncend] = data[varname][modstart+mo:modend+mo]
                        except ValueError:
                            ncend = ncstart + data[varname][modstart+mo:modend+mo].shape[0]
                            ncid[modstart/fstep,ncstart:ncend] = data[varname][modstart+mo:modend+mo]
                        except TypeError:
                            print 'variable %s not available from model' %varname
                            continue
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
                    lead_time = time[i] - modrun
                    lead_time = lead_time.total_seconds() /3600.
                    lead_index = int(sp.floor(lead_time /fstep) )
                    timeofmonth = time[i] - dt.datetime(year,month,1)
                    timeofmonth = timeofmonth.total_seconds() /3600.
                    if not ((timeofmonth<0) or (timeofmonth>(daysofmonth*24-1)) or (lead_index > int(modlength/fstep))):
                       ncid[lead_index, int(timeofmonth)] = var[i]

    return 0

def ensure_hourly(data):
	'''
	'''
	time = data['time']
	try:
		if not ((len(time)-1) == (time[-1]-time[0]).total_seconds()/3600):
			make_hourly(data)
	except AttributeError:
		return 1
	return 0

def make_hourly(data):
	'''
	'''
	time = data['time']
	ntimesnew = int((time[-1]-time[0]).total_seconds() / 3600)
	timenew = [time[0]+dt.timedelta(i)/24 for i in range(ntimesnew)]
	for varname, var in data.iteritems():
		if (varname == 'time'):
			continue
		var = sp.interp(pl.date2num(timenew), pl.date2num(time), var)
                data[varname] = var
	data['time'] = timenew
	return 0

