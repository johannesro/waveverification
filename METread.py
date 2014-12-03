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
import datetime as dt
import scipy.ndimage as nd
import os
import d22 # modlue for reading .d22 files from offshore stations
import netCDF4 as nc4 #import Dataset, num2date, date2num
import MySQLdb


# variable list for ECMWF netcdf files
ECvardict ={ 'Hs':  'swh', 'Tp':  'pp1d', 'Tm02': 'mp2', 'DDM': 'mwd', 'Hs_s': 'shts', 'Tm02_s': 'p2ps', 'DDM_s': 'mdts', 'FF': 'wind'}


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
    # check file, preferably from opdata
    filename=run.strftime("/opdata/arome2_5_main/arome_metcoop_default2_5km_%Y%m%d_%H.nc")
    if not os.path.exists(filename):
        filename=run.strftime("/starc/DNMI_AROME_METCOOP/%Y/%m/%d/AROME_MetCoOp_%H_DEF.nc_%Y%m%d")
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
       nan = sp.ones(67)*sp.nan
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
    filename=run.strftime("/disk4/ECMWF_ECWAVE/ECWAM_wave_%Y%m%d_%H.nc")
    print(' ')
    print('reading '+filename)
    if os.path.isfile(filename):
        ecdatadict = nctimeseries(filename, ECvardict.values(), location)
        datadict = {'time':ecdatadict['time']}
        for varname, ecname in ECvardict.iteritems():
            datadict.update({varname: ecdatadict[ecname]})
    else:
       print('file '+filename+' does not exist; returning missing values')
       nan = sp.ones(67)*sp.nan
       datadict = {'time':nan}
       for varname in ECvardict.keys():
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
    filename=run.strftime("/starc/DNMI_ROMS/%Y/%m/%d/Nordic-4km_stormsurge_%H.nc_%Y%m%d")
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



def nctimeseries(filename, varnamelist, coordinate):
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
    try: 
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
    except KeyError:
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
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
        except KeyError:
            data[var] = None
    # check varnames 
    varnames = {'Hs':'hs', 'Tp':'tp', 'FF':'ff','DD':'dd','Tm02':'tm2', 'DDM':'thq'}
    for varname,altname in varnames.iteritems():
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

    quba = MySQLdb.connect(host="qubadb",user="qbload", passwd="load", db="quba",port=19100)
    qhis = MySQLdb.connect(host="qubadb",user="qbload", passwd="load", db="qhist",port=19100)
    cquba = quba.cursor()
    cqhis = qhis.cursor()
    cquba.execute("select stationid from station where name like '"+station+"'")
    staid = str(cquba.fetchall()[0][0])
    #print('Station ID for '+station+': '+staid)

# loop over variables
    vardict  = {}
    for var in varlist:
# fetch all forcasts of the variable at given valid time as function of forecast time
        #print('read from quba')
        cquba.execute("select valid,value from subjective where run = \""+runstr+"\" and stationid="+staid+" and pindexid = "+str(varid[var]))
        x1 = cquba.fetchall()
        ftime1 = [ x1[i][0] for i in range(len(x1))] # 
        value1 = [ x1[i][1] for i in range(len(x1))] # 

        #print('read from qhis')
        cqhis.execute("select valid,value from subjective where run = \""+runstr+"\" and stationid="+staid+" and pindexid = "+str(varid[var]))
        x2 = cqhis.fetchall()
        ftime2 = [ x2[i][0] for i in range(len(x2))] # 
        value2 = [ x2[i][1] for i in range(len(x2))] # 

# put data from quba and qhis together:
        time = ftime2+ftime1
        value = value2+value1
        vardict.update({var:(time,value)})
    return vardict


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


