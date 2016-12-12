#
# File reading routines for time series of MET Norway data
#


# common syntax for all readers:
#
# def source(location, run, varnamelist):
#    return vardict
#
# location:    (lat, lon) touple or list
# run:   datetime object that identifies the initialization time of the model to be read
# varnamelist: list of variable names to read, i.e. ['FF', 'DD'] (currently only implemented for netcdf files)
#
# returns:
# datadict:     dictionary with time and timeseries


import numpy as np
import numpy.ma as ma
import scipy as sp
import scipy.spatial as spatial
import scipy.interpolate as interpolate
import datetime as dt
import scipy.ndimage as nd
import os
import d22 # modlue for reading .d22 files from offshore stations
import netCDF4 as nc4 #import Dataset, num2date, date2num
import MySQLdb
import dataanalysis as da
import pylab as pl
from mpl_toolkits.basemap import Basemap #package for map projection
#import sphrot


# variable list for ECMWF netcdf files
# should better list these names in variables.py!
ECvardict ={ 'Hs':  'significant_wave_height', 
        'Tp':  'peak_wave_period', 
	'Tm02': 'mean_wave_period',
	'DDM': 'wave_direction',
	'Hs_s': 'significant_swell_wave_height',
	'Tm02_s': 'sea_surface_swell_wave_period', # or is this Tp_s?
	'DDM_s': 'sea_surface_swell_wave_to_direction'}

MWAMvardict = { 'Hs':'hs', 'Tp':'tp', 'Tp_s':'tp_swell', 'Tm02':'tm2', 'DDM':'thq', 'Hs_s':'hs_swell', 'Tm02_s':'tm2_swell', 'DDM_s':'thq_swell', 'FF':'ff', 'DD':'dd'}

MWAM8vardict = { 'Hs':'VHM0', 'Tp':'VTPK', 'Tm02':'VTM02', 'DDM':'VMDR', 'Hs_s':'VHM0_SW', 'Tm02_s':'VTM02_SW', 'DDM_s':'VMDR_SW', 'DDP':'VPED', 'FF':'ff', 'DD':'dd'}


WAMAROMEvardict = {'Hs': 'significant_wave_height', 'FF':'wind_speed_10m'}

def obs_d22(station, day, numdays=1, varlist=[]):
    '''
    station: offshore station name as string
    day: datetime object (should be with 00 hours)
    returns dictionary with hourly time, Hs, Tp, Tm, FF, DD
    each of the returned variable is a 2d array with time and the number of available sensors
    '''
    time, WMlist, WIlist = d22.read_d22(station,start=day, end=day+dt.timedelta(numdays-1))
    datadict = {'time':time['1hr']}

    # organize each variable as 2d-arrays with time and sensor number as dimensions
    Hs, Tp, Tm02 ,Tm01, DDP, DDM, WMnames, FF, DD, WInames = [],[],[],[], [], [], [], [], [], []
    for WM in WMlist:
        if sp.isfinite(WM['Hs']).sum()>1:
            Hs.append(WM['Hs'])
            Tp.append(WM['Tp'])
            Tm02.append(WM['Tm02'])
            Tm01.append(WM['Tm01'])
            DDP.append(WM['DDP'])
            DDM.append(WM['DDM'])
            WMnames.append(WM['name'])

    for WI in WIlist:
        if sp.isfinite(WI['FF']).sum()>1:
            FF.append(WI['FF'])
            DD.append(WI['DD'])
            WInames.append(WI['name'])

    datadict.update({'Hs':sp.array(Hs), 'Tp':sp.array(Tp), 'Tm02':sp.array(Tm02) , 'Tm01':sp.array(Tm01) , 
                     'FF':sp.array(FF), 'DD':sp.array(DD),'DDM':sp.array(DDM), 'DDP':sp.array(DDP)})
    datadict.update({'WMnames':WMnames, 'WInames':WInames}) # save sensor names

    return datadict


def AROME_modrun(location, run, varnamelist, step=1):
    '''
    location: (lat, lon) touple or list
    run: datetime object of model run initialization time (should be 00, 03, 06,... hours)
    '''
    # ensure correct formats
    location = list(location)
    varnamelist = varnamelist+['x_wind_10m','y_wind_10m']
    # check file, preferably from opdata -> /vol/data from ppi!!!
    filename=run.strftime("/vol/data/arome2_5_main/arome_metcoop_default2_5km_%Y%m%d_%H.nc")
    if not os.path.exists(filename):
        filename=run.strftime("/lustre/storeB/immutable/short-term-archive/DNMI_AROME_METCOOP/%Y/%m/%d/AROME_MetCoOp_%H_DEF.nc_%Y%m%d")
        #filename=run.strftime("/starc/DNMI_AROME_METCOOP/%Y/%m/%d/AROME_MetCoOp_%H_DEF.nc_%Y%m%d")
    print(' ')
    print('reading '+filename)
    if os.path.isfile(filename):
        datadict = nctimeseries(filename, varnamelist, location)
        datadict['DD'], datadict['FF'] = DD_FF(datadict['x_wind_10m'],datadict['y_wind_10m'])
    else:
       print('file '+filename+' does not exist; returning missing values')
       nan = sp.ones(67)*sp.nan
       datadict = {'time':nan,'FF':nan, 'DD':nan, 'x_wind_10m':nan, 'y_wind_10m':nan}
    return datadict

def WAMAROME2W_modrun(location, run, varnamelist, step=1):
    '''
    location: (lat, lon) touple or list
    run: datetime object of model run initialization time (should be 00, 03, 06,... hours)
    '''
    # ensure correct formats
    location = list(location)
    #varnamelist = varnamelist+['x_wind_10m','y_wind_10m']
    # check file, preferably from opdata
    filename=run.strftime("/vol/fou/atmos2/jakobks/MyWave/WAM_2WCoupled/WAM_2WCoupled_%Y%m%d00.nc")
    print(' ')
    print('reading '+filename)
    if os.path.isfile(filename):
        #datadict = nctimeseries(filename, WAMAROMEvardict.values(), location)
        WAdatadict = nctimeseries(filename, WAMAROMEvardict.values(), location)
        datadict = {'time':WAdatadict['time']}
        for varname, WAname in WAMAROMEvardict.iteritems():
            datadict.update({varname: WAdatadict[WAname]})
        #datadict['DD'], datadict['FF'] = DD_FF(datadict['x_wind_10m'],datadict['y_wind_10m'])
    else:
       print('file '+filename+' does not exist; returning missing values')
       nan = sp.ones(49)*sp.nan
       datadict = {'time':nan,'FF':nan, 'Hs':nan}#, 'x_wind_10m':nan, 'y_wind_10m':nan}
    return datadict

def WAMAROME1W_modrun(location, run, varnamelist, step=1):
    '''
    location: (lat, lon) touple or list
    run: datetime object of model run initialization time (should be 00, 03, 06,... hours)
    '''
    # ensure correct formats
    location = list(location)
    #varnamelist = varnamelist+['x_wind_10m','y_wind_10m']
    # check file, preferably from opdata
    filename=run.strftime("/vol/fou/atmos2/jakobks/MyWave/WAM_1WCoupled/WAM_1WCoupled_%Y%m%d00.nc")
    print(' ')
    print('reading '+filename)
    if os.path.isfile(filename):
        #datadict = nctimeseries(filename, varnamelist, location)
        WAdatadict = nctimeseries(filename, WAMAROMEvardict.values(), location)
        datadict = {'time':WAdatadict['time']}
        for varname, WAname in WAMAROMEvardict.iteritems():
            datadict.update({varname: WAdatadict[WAname]})
        #datadict['DD'], datadict['FF'] = DD_FF(datadict['x_wind_10m'],datadict['y_wind_10m'])
    else:
       print('file '+filename+' does not exist; returning missing values')
       nan = sp.ones(49)*sp.nan
       datadict = {'time':nan,'FF':nan, 'Hs':nan}#, 'x_wind_10m':nan, 'y_wind_10m':nan}
    return datadict


def LAWAM_modrun(location, run, varnamelist, step=1):
    '''
    location: (lat, lon) touple or list
    run: datetime object of model run initialization time (should be 00, 03, 06,... hours)
    '''
    # ensure correct formats
    location = list(location)
    # check file, preferably from opdata
    filename=run.strftime("/disk4/ECMWF_LAWAM/LAW_wave_%Y%m%d_%H.nc")
    print(' ')
    print('reading '+filename)
    if os.path.isfile(filename):
        ecdatadict = nctimeseries(filename, ECvardict.values(), location)
        datadict = {'time':ecdatadict['time']}
        for varname, ecname in ECvardict.iteritems():
            datadict.update({varname: ecdatadict[ecname]})
    else:
       print('file '+filename+' does not exist; returning missing values')
       nan = sp.ones(241)*sp.nan
       datadict = {'time':nan}
       for varname in ECvardict.keys():
           datadict.update({varname: nan})
    return datadict



def ECWAM_modrun(location, run, varnamelist, step=1):
    '''
    location: (lat, lon) touple or list
    run: datetime object of model run initialization time (should be 00, 03, 06,... hours)
    '''
    # ensure correct formats
    location = list(location)
    # check file, preferably from opdata
    filename=run.strftime("/vol/data/ec/ecwam_%Y%m%dT%HZ.nc")
    print(' ')
    print('reading '+filename)
    if os.path.isfile(filename):
        ecdatadict = nctimeseries(filename, ECvardict.values(), location)
        datadict = {'time':ecdatadict['time']}
        for varname, ecname in ECvardict.iteritems():
            datadict.update({varname: ecdatadict[ecname]})
    else:
       print('file '+filename+' does not exist; returning missing values')
       nan = sp.ones(241)*sp.nan
       datadict = {'time':nan}
       for varname in ECvardict.keys():
           datadict.update({varname: nan})
    return datadict

def MWAM10_modrun(location, run, varnamelist, step=1):
    '''
    location: (lat, lon) touple or list
    run: datetime object of model run initialization time (should be 00, 03, 06,... hours)
    '''
    # ensure correct formats
    location = list(location)
    # check file, preferably from opdata
    filename = run.strftime("/vol/data/wave/MyWave_wam10_WAVE_%Y%m%dT%HZ.nc")
    if not os.path.isfile(filename):
        filename = run.strftime("/vol/hindcast3/waveverification/mywaveWAM_archive/MyWave_wam4_WAVE_%Y%m%dT%HZ.nc")
    gfile = nc4.Dataset('/lustre/storeA/project/fou/hi/waveverification/WAM10grid.nc')  #this step could later be droppet, if lat/lon info are available in operational files
    lat, lon = gfile.variables['latitude'][:], gfile.variables['longitude'][:]
    gfile.close()
    print(' ')
    print('reading '+filename)
    if os.path.isfile(filename):
        MWdatadict = nctimeseries(filename, MWAMvardict.values(), location, grid=(lat,lon))
        #print MWdatadict
        datadict = {'time': MWdatadict['time']}
        for varname, MWAMname in MWAMvardict.iteritems():
            datadict.update({varname: MWdatadict[MWAMname]})
    else:
       print('file '+filename+' does not exist; returning missing values')
       nan = sp.ones(67)*sp.nan
       datadict = {'time':nan}
       for varname in MWAMvardict.keys():
           datadict.update({varname: nan})
    return datadict



def MWAM4_modrun(location, run, varnamelist, step=1):
    '''
    location: (lat, lon) touple or list
    run: datetime object of model run initialization time (should be 00, 03, 06,... hours)
    '''
    # ensure correct formats
    location = list(location)
    # check file, preferably from opdata
    filename = run.strftime("/vol/data/wave/MyWave_wam4_WAVE_%Y%m%dT%HZ.nc")
    if not os.path.isfile(filename):
        filename = run.strftime("/vol/hindcast3/waveverification/mywaveWAM_archive/MyWave_wam4_WAVE_%Y%m%dT%HZ.nc")
    if not os.path.isfile(filename):
        #filename = run.strftime("/starc/DNMI_WAVE/%Y/%m/%d/MyWave_wam4_WAVE_%Y%m%dT%HZ.nc")
        filename = run.strftime("/lustre/storeB/immutable/short-term-archive/DNMI_WAVE/%Y/%m/%d/MyWave_wam4_WAVE_%Y%m%dT%HZ.nc")
    print(' ')
    print('reading '+filename)
    if os.path.isfile(filename):
        MWdatadict = nctimeseries(filename, MWAMvardict.values(), location)#, grid=(lat,lon))
        datadict = {'time': MWdatadict['time']}
        for varname, MWAMname in MWAMvardict.iteritems():
            datadict.update({varname: MWdatadict[MWAMname]})
    else:
       print('file '+filename+' does not exist; returning missing values')
       nan = sp.ones(67)*sp.nan
       datadict = {'time':nan}
       for varname in MWAMvardict.keys():
           datadict.update({varname: nan})
    return datadict

def EXP_modrun(location, run, varnamelist, step=1):
    '''
    location: (lat, lon) touple or list
    run: datetime object of model run initialization time (should be 00, 03, 06,... hours)
    '''
    # ensure correct formats
    location = list(location)
    # check file, preferably from opdata
    filename = run.strftime("/vol/data/wave_exp/MyWave_wam4_WAVE_%Y%m%dT%HZ.nc")
    print(' ')
    print('reading '+filename)
    if os.path.isfile(filename):
        MWdatadict = nctimeseries(filename, MWAMvardict.values(), location) #, grid=(lat,lon))
        #print MWdatadict
        datadict = {'time': MWdatadict['time']}
        for varname, MWAMname in MWAMvardict.iteritems():
            datadict.update({varname: MWdatadict[MWAMname]})
    else:
       print('file '+filename+' does not exist; returning missing values')
       nan = sp.ones(67)*sp.nan
       datadict = {'time':nan}
       for varname in MWAMvardict.keys():
           datadict.update({varname: nan})
    return datadict


def MWAM8_modrun(location, run, varnamelist, step=1):
    '''
    location: (lat, lon) touple or list
    run: datetime object of model run initialization time (should be 00, 03, 06,... hours)
    Note that the ARC-MFC WAM always starts at 06 hours only!
    '''
    # ensure correct formats
    location = list(location)
    # check file, preferably from opdata (from ppi: /vol/data/)
    filename = run.strftime("/vol/data/wave/MyWave_wam8_WAVE_%Y%m%dT06Z.nc")
    if not os.path.isfile(filename):
        #filename = run.strftime("/starc/DNMI_WAVE/%Y/%m/%d/MyWave_wam8_WAVE_%Y%m%dT%HZ.nc")
        filename = run.strftime("/lustre/storeB/immutable/short-term-archive/DNMI_WAVE/%Y/%m/%d/MyWave_wam8_WAVE_%Y%m%dT06Z.nc")
    print(' ')
    print('reading '+filename)
    if os.path.isfile(filename):
        MWdatadict = nctimeseries(filename, MWAM8vardict.values(), location)#, grid=(lat,lon))
        datadict = {'time': MWdatadict['time']}
        for varname, MWAMname in MWAM8vardict.iteritems():
            datadict.update({varname: MWdatadict[MWAMname]})
    else:
       print('file '+filename+' does not exist; returning missing values')
       nan = sp.ones(67)*sp.nan
       datadict = {'time':nan}
       for varname in MWAMvardict.keys():
           datadict.update({varname: nan})
    return datadict




def Nordic4_stormsurge_modrun(location, run, varnamelist, step=1):
    '''
    location: (lat, lon) touple or list
    run: datetime object of model run initialization time (should be 00, 03, 06,... hours)
    varnamelist: list of variables, i.e. h, ubar, vbar, zeta_detided_raw, zeta_detided, h
    '''
    # ensure correct formats
    location = list(location)
    # check file, preferably from opdata
    #filename=run.strftime("/starc/DNMI_ROMS/%Y/%m/%d/Nordic-4km_stormsurge_%H.nc_%Y%m%d")
    filename=run.strftime("/lustre/storeB/immutable/short-term-archive/DNMI_ROMS/%Y/%m/%d/Nordic-4km_stormsurge_%H.nc_%Y%m%d")
    print(' ')
    print('reading '+filename)
    if os.path.isfile(filename):
        datadict = nctimeseries(filename, varnamelist, location)
    else:
       print('file '+filename+' does not exist; returning missing values')
       nan = sp.ones(144)*sp.nan
       datadict = {'time':nan}
       for varname in varnamelist:
           datadict.update({varname: nan})
    return datadict



def nctimeseries(filename, varnamelist, coordinate, grid=None):
    '''
    Extract time series from the nearest model point at given coordinate from netCDF file
    attempts to read all variables in varnamelist. Returns None if a variable doesnt exist
    coordinate: tuple of (lat, lon)
    varnamelist: list of strings, i.e. ['x_wind_10m', 'y_wind_10m']
    returns a dictionary containing the time series of each variable as well as time list
    '''
    lat0, lon0 = coordinate[0], coordinate[1]
    try:
        nc = nc4.Dataset(filename, 'r') # do only if not already opened
    except TypeError:
        nc = filename
    nctime = nc.variables['time'] # use this if we're dealing with single netcdf file
    time = nc4.num2date(nctime[:], units=nctime.units)
    if grid==None:
        try: 
            lon = nc.variables['longitude'][:]
            lat = nc.variables['latitude'][:]
        except KeyError:
            lon = nc.variables['lon'][:]
            lat = nc.variables['lat'][:]
    else:
        lat,lon = grid[0], grid[1]
    try: # use this for generic grids that provide full lon and lat arrays
        y,x = find_pos(lon, lat, lon0, lat0)
    except ValueError:  # use this if we have a lonlat grid with lon and lat vector
        x = find_pos1d1(lon, lon0)
        y = find_pos1d1(lat, lat0)
    data = {'time':time}
    for var in varnamelist:
        try:
            if len(nc.variables[var].shape) == 3:
                data[var] = nc.variables[var][:,y,x]
            if len(nc.variables[var].shape) == 4:
                data[var] = nc.variables[var][:,0,y,x]
        except KeyError:
            #print 'no variable '+var+' in file '+filename
            data[var] = None
    # check for alternative varnames 
    varnames = {'hs':'Hs', 'tp':'Tp', 'ff':'FF', 'dd':'DD', 'tm2':'Tm02', 'thq':'DDM', 'significant_wave_height':'Hs', 'wind_speed_10m':'FF'}
    for altname,varname in varnames.iteritems():
        try:
            data[varname]=data[altname]
        except KeyError:
            continue
    return data

def find_pos1d1(x,x0):
    x=sp.array(x)
    xdist=(x-x0)**2
    flist=list(xdist[sp.isfinite(xdist)])
    flist.sort()
    return list(xdist).index(flist[0])

def find_pos1d(x,y,x0,y0):
    '''find index in x,y that suits the point x0,y0 most
    x, y are one-dimensional arrays or sequences of x and y coordinates
    returns the best fitting index of x,y as integer'''
    x,y = sp.array(x), sp.array(y)
    xdist=x-x0
    ydist=y-y0
    f=xdist**2+ydist**2 #calculate distance**2
    flist=list(f[sp.isfinite(f)]) # make a list without NaNs
    flist.sort() # sort the list, so that the smallest distance**2 is at the top
    return list(f).index(flist[0]) # return index of closest distance**2



def find_pos(lon,lat,lon_0,lat_0):
    '''find the best fitting position of a coordinate pair in an image
    lon: longitude image
    lat: latitude image
    lon_0, lat_0: coordinates of desired position; 
    f: scale factor to speed up the calculation, but reduces the accurance of the output position; 
    returns a tuple y,x which gives the best position in lon,lat of lon_0,lat_0'''
    best_lat_i=lat.copy()
    best_lat_i[best_lat_i>lat_0]=2*lat_0-best_lat_i[best_lat_i>lat_0]
    best_lon_i=lon.copy()
    best_lon_i[best_lon_i>lon_0]=2*lon_0-best_lon_i[best_lon_i>lon_0]
    best_position=best_lat_i/lat_0.__abs__()+best_lon_i/lon_0.__abs__()
    best_y,best_x=nd.measurements.maximum_position(best_position)
    best_y,best_x=best_y,best_x
    return best_y,best_x

def DD_FF(u,v,met=True):
    ''' calculates wind/current speed and direction from u and v components
    #
    if u and v are easterly and northerly components,
    returns:
    FF : wind speed
    DD : wind direction in meteorological standard (directions the wind is coming from)
    call with met=False for oceanographic standard
    '''
    if met==False:
        u,v = -u, -v
    DD = ma.arctan2(-u, -v)*180/sp.pi
    DD[DD < 0] = 360 + DD[DD < 0]
    FF = ma.sqrt(u**2 + v**2)
    return DD, FF


def subjectiveforecast(runtime, station, varlist):
    varid={'Hs': 200, 'FF':298, 'DD':299, 'Tp':201, 'Tm02':203}
    runstr = runtime.isoformat(sep=' ')
    if station=='ekofisk' or station=='ekofiskL':
        #print('ekofisk=ekofisk_os')
        station='Ekofisk_OS'

    quba = MySQLdb.connect(host="kommersiellqubadb",user="qbload", passwd="load", db="kommersiellquba")
    #qhis = MySQLdb.connect(host="qubadb",user="qbload", passwd="load", db="qhist",port=19100)
    cquba = quba.cursor()
    #cqhis = qhis.cursor()
    cquba.execute("select stationid from station where name like '"+station+"'")
    staid = str(cquba.fetchall()[0][0])
    #print('Station ID for '+station+': '+staid)

# loop over variables
    vardict  = {}
    for var in varlist:
# fetch all forcasts of the variable at given valid time as function of forecast time
        #print('read from quba')
        cquba.execute("select valid,value from subjective where levelid=0 and run = \""+runstr+"\" and stationid="+staid+" and pindexid = "+str(varid[var]))
        x1 = cquba.fetchall()
        ftime1 = [ x1[i][0] for i in range(len(x1))] # 
        value1 = [ x1[i][1] for i in range(len(x1))] # 

        #print('read from qhis')
        #cqhis.execute("select valid,value from subjective where levelid=0 and run = \""+runstr+"\" and stationid="+staid+" and pindexid = "+str(varid[var]))
        #x2 = cqhis.fetchall()
        #ftime2 = [ x2[i][0] for i in range(len(x2))] # 
        #value2 = [ x2[i][1] for i in range(len(x2))] # 

# put data from quba and qhis together:
        time = ftime1
        value = value1
        vardict.update({var:(time,value)})
    return vardict

class stationfiles():
    '''
    Class to read a list of subsequent validation station files
    example:
    sfiles = stationfils(glob.glob('ekofisk_2014??.nc'))

    '''
    def __init__(self, files):
        self.time, self.nclist = [], []
        files.sort()
        for filename in files:
            nc = nc4.Dataset(filename,'r')
            self.nclist.append(nc)
            self.time=self.time+list(nc4.num2date(nc.variables['time'],nc.variables['time'].units))

    def get_var(self, varname, groupname,treshhold=None):
        var = []
        for nc in self.nclist:
            ncg = nc.groups[groupname]
            var.append(ncg.variables[varname][:])
        var = np.ma.concatenate(var, axis=1)
        var[var.mask==True]=sp.nan
        if not treshhold==None:
            var[var > treshhold]=sp.nan
        return var



'''
def nctimeseries(filename, varnamelist, coordinate):
   # Extract time series from the nearest model point at given coordinate from netCDF file
   # attempts to read all variables in varnamelist. Returns None if a variable doesnt exist
   # coordinate: tuple of (lat, lon)
   # varnamelist: list of strings, i.e. ['x_wind_10m', 'y_wind_10m']
   # returns a dictionary containing the time series of each variable as well as time list
 
    lat0, lon0 = coordinate[0], coordinate[1]
    nc = nc4.Dataset(filename, 'r')
    nctime = nc.variables['time']
    time = nc4.num2date(nctime[:], units=nctime.units)
    lon = nc.variables['longitude'][:] # this is probably a waist of time. should only read x and y positions on first use of model
    lat = nc.variables['latitude'][:]
    y,x = find_pos(lon, lat, lon0, lat0)
    data = {'time':time}
    for var in varnamelist:
        try:
            if len(nc.variables[var].shape) == 3:
                data[var] = nc.variables[var][:,y,x]
        except KeyError:
            data[var] = None
    return data

'''


def north_direction(lat):
    '''get the north direction relative to image positive y coordinate'''
    dlatdx = nd.filters.sobel(lat,axis=1,mode='constant',cval=sp.nan) #gradient in x-direction
    dlatdy = nd.filters.sobel(lat,axis=0,mode='constant',cval=sp.nan)
    ydir = lat[-1,0] -lat[0,0] # check if latitude is ascending or descending in y axis
    # same step might have to be done with x direction.
    return sp.arctan2(dlatdx,dlatdy*sp.sign(ydir) )*180/sp.pi

def rotate_wind(U,V,alpha):
    alpha = sp.array(alpha)*sp.pi/180
    alpha = alpha.flatten()
    R = sp.array([[sp.cos(alpha), -sp.sin(alpha)], [sp.sin(alpha), sp.cos(alpha)] ])
    shpe = U.shape
    origwind = sp.array((U.flatten(), V.flatten()))
    if len(R.shape)==2:
        rotwind = dot(R, origwind) # for constant rotation angle
    else:
        # for rotation angle given as array with same dimensions as U and V:
        # k-loop with rotwind(k) = dot(R(i,j,k), origwind(j,k)) (einstein summation indices)
        rotwind = sp.einsum("ijk,ik -> jk", R, origwind) 
    Urot ,Vrot = rotwind[0,:], rotwind[1,:]
    Urot = Urot.reshape(shpe)
    Vrot = Vrot.reshape(shpe)
    return Urot, Vrot


class ncfile():
    '''
    Class to handle netcdf files with geophysical fields
    '''
    def __init__(self,filename):
        try:
            self.nc = nc4.Dataset(filename, 'r') # do only if not already opened
        except RuntimeError:
            self.nc = nc4.MFDataset(filename, 'r')
        except TypeError:
            self.nc = filename
        nctime = self.nc.variables['time'] # use this if we're dealing with single netcdf file
        self.time = nc4.num2date(nctime[:], units=nctime.units)
        try: 
            self.lon = self.nc.variables['longitude'][:]
            self.lat = self.nc.variables['latitude'][:]
        except KeyError:
            self.lon = self.nc.variables['lon'][:]
            self.lat = self.nc.variables['lat'][:]
        self.spatial=False
    def __repr__(self):
        print 'advanced nc file interface by johannesro@met.no'
        #print self.nc
        ostr = 'lon: '+ str(self.lon.min()) + ' ' + str(self.lon.max()) +'\n'
        ostr = ostr + 'lat: '+ str(self.lat.min()) + ' ' + str(self.lat.max()) +'\n'
        ostr = ostr + 'time: ' + str(self.time[0]) +' - '+ str(self.time[-1])
        return ostr
    def close(self):
        self.nc.close()
    def field(self, timestep, varname):
        '''
        extract field at timestep
        '''
        timei = da.find_pos1d1(pl.date2num(self.time),pl.date2num(timestep))
        return self.nc.variables[varname][timei]
    def uv(self, timestep):
        '''
        extract u and v components as eastern/northern components; as fields at timestep
        '''
        northdir = north_direction(self.lat)
        self.northdir= northdir
        try:
            u_grid = self.field(timestep, 'Uwind')
            v_grid = self.field(timestep, 'Vwind')
        except KeyError:
            try:
                u_grid = self.field(timestep, 'x_wind')
                v_grid = self.field(timestep, 'y_wind')
            except KeyError:
                u_grid = self.field(timestep, '10u')
                v_grid = self.field(timestep, '10v')
        u,v = rotate_wind(u_grid,v_grid,northdir)
        return u,v 
    def trajectory(self,lon, lat, time, varnamelist,interpolate=True):
        '''
        pick selected variables (provided as list of strings) along a trajectory (given by lon, lat, time)
        '''
        result = {}
        for varname in varnamelist:
            result.update({varname:sp.zeros_like(lon)})
        if self.spatial==False:
            self.initiatespatial()
        x0,y0 = self.pMap(lon, lat)
        points = zip(x0,y0) # points in trajectory           #### xy order!
        d,p = self.Tree.query(points,k=1) #nearest point, gives distance and point index
        for i in range(len(time)):
            timei = findtime(self.time, time[i])
            xi, yi = self.xip[p[i]], self.etap[p[i]] # get grid indices 
            #select variables
            for varname in varnamelist:
                var = self.nc.variables[varname]
                if interpolate==True:
                    try:
                        yis, xis, timeis = slice(yi-1,yi+2), slice(xi-1,xi+2), slice(timei-1,timei+2)
                        values = var[timeis,yis,xis] # select 3x3x3 box around closest value  #####xy order!
                        # space-interplote each time slice
                        valuet = []
                        for ii in range(3):
                            valuet.append(int_point(self.y[yis,xis], self.x[yis,xis], values[ii], y0[i], x0[i]))
                        # time interpolation 
                        # value = int_time(self.time[timeis], sp.array(valuet), time[i])
                        value = np.interp(pl.date2num(time[i]), pl.date2num(self.time[timeis]), sp.array(valuet))
                    except IndexError:
                        value = var[timei,yi,xi]
                        print('Warning: closest point is located on edge of model domain!')
                else:
                    value = var[timei,yi,xi]
                result[varname][i]=value
        return result
    def initiatespatial(self):
        self.spatial=True
        self.pMap = Basemap(projection='stere',lat_0=68.0,lon_0=15.0,llcrnrlon=3.,llcrnrlat=60.3,urcrnrlon=47.0,urcrnrlat=71.,resolution='c')
        self.x, self.y = self.pMap(self.lon, self.lat)
        XI,ETA = sp.meshgrid(range(self.y.shape[1]),range(self.y.shape[0]))
        self.xip, self.etap = XI.ravel(), ETA.ravel()
        self.Tree = spatial.KDTree(zip(self.x.ravel(), self.y.ravel()))
    def timeseries(filename, varnamelist, coordinate):
        '''
        Extract field at time step from netCDF file
        attempts to read all variables in varnamelist. Returns None if a variable doesnt exist
        coordinate: tuple of (lat, lon)
        varnamelist: list of strings, i.e. ['x_wind_10m', 'y_wind_10m']
        returns a dictionary containing the time series of each variable as well as time list
        '''
        time = lat0, lon0 = coordinate[0], coordinate[1]
        try: # use this for generic grids that provide full lon and lat arrays
            y,x = find_pos(lon, lat, lon0, lat0)
        except ValueError:  # use this if we have a lonlat grid with lon and lat vector
            x = find_pos1d1(lon, lon0)
            y = find_pos1d1(lat, lat0)
        data = {'time':time}
        for var in varnamelist:
            try:
                if len(nc.variables[var].shape) == 3:
                    data[var] = nc.variables[var][:,y,x]
                if len(nc.variables[var].shape) == 4:
                    data[var] = nc.variables[var][:,0,y,x]
            except KeyError:
                data[var] = None
        # check for alternative varnames 
        varnames = {'hs':'Hs', 'tp':'Tp', 'ff':'FF', 'dd':'DD', 'tm2':'Tm02', 'thq':'DDM', 'significant_wave_height':'Hs', 'wind_speed_10m':'FF'}
        for altname,varname in varnames.iteritems():
            try:
                data[varname]=data[altname]
            except KeyError:
                continue
        return data

class Nora10file(ncfile):
    '''
    Class to handle Nora10 netcdf files that don't include longitute/latitde fields
    '''
    def __init__(self,filename):
        try:
            self.nc = nc4.Dataset(filename, 'r') # do only if not already opened
        except RuntimeError:
            self.nc = nc4.MFDataset(filename, 'r')
        except TypeError:
            self.nc = filename
        nctime = self.nc.variables['time'] # use this if we're dealing with single netcdf file
        self.time = nc4.num2date(nctime[:], units=nctime.units)
        gridf = nc4.Dataset('/home/johannesro/Nora10/grid/Nora10grid.nc','r')
        self.lon = gridf.variables['longitude'][:]
        self.lat = gridf.variables['latitude'][:]
        self.spatial = False
    def uv(self, timestep):
        '''
        extract u and v components as eastern/northern components; as fields at timestep
        '''
        #northdir = north_direction(self.lat)
        #self.northdir= northdir
        try:
            u = self.field(timestep, 'Uwind')
            v = self.field(timestep, 'Vwind')
        except KeyError:
            try:
                u = self.field(timestep, 'x_wind')
                v = self.field(timestep, 'y_wind')
            except KeyError:
                u = self.field(timestep, '10u')
                v = self.field(timestep, '10v')
        #print 'wind not rotated'
        return u,v 


def int_point_m(lonfield,latfield,field,lon,lat,Imap):
    '''use the LinearNDInterpolater to interpolate one point in a field

    lonfield, latfield: coordinates of field
    field: 2d-field of values
    lon, lat: point to which field is interpolated
    Imap: basemap instance
    '''
    # get coordinates to cartesian coordinate system
    xfield,yfield = Imap(lonfield, latfield) 
    x,y = Imap(lon,lat)
    # do the interpolation (first define points for interpolation
    points = sp.array( [xfield.flatten(), yfield.flatten()] ).transpose()
    I = interpolate.LinearNDInterpolator(points,field.flatten())
    return I(x, y)

def int_point(xfield,yfield,field,x,y):
    '''use the LinearNDInterpolater to interpolate one point in a field
    '''
    points = sp.array( [xfield.flatten(), yfield.flatten()] ).transpose()
    I = interpolate.LinearNDInterpolator(points,field.flatten())
    return I(x, y)

def int2_time(time, y, time_0):
    f = interpolate.interp1d(pl.date2num(time), y, kind='linear')
#    print time, time_0
    return f(pl.date2num(time_0))

def int_time(time, y, time_0):
    #print time, time_0
    return np.interp(pl.date2num(time_0), pl.date2num(time), y)


def findtime(time,time0):
    i = da.find_pos1d1(pl.date2num(time),pl.date2num(time0))
    return i

def gsplitread(station):
    '''read Nora10 data from gsplit files for a given station name'''
    f = '/home/johannesro/WAM_NORA10/gsplit_'+station+'.txt'
    data = pl.loadtxt(f, skiprows=4).transpose()
    dd={}
    dd['time'] = [dt.datetime(int(data[0,i]), int(data[1,i]), int(data[2,i]), int(data[3,i]), 0, 0) for i in range(data.shape[1])]
    dd['Hs'] = data[6]
    dd['Tp'] = data[7]
    dd['FF'] = data[4]
    dd['DDM'] = data[10]
    dd['Tm02'] = data[8] 
    return dd


