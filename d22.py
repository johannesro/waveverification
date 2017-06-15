#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Read d22 files from /opdata
# use d22plot.py for plotting
#
# ole.aarnes@met.no
# johannes.rohrs@met.no

######################################################################
import numpy as np
import scipy as sp
import scipy.ndimage as ndimage
import sys
import datetime as dt
import pylab as pl
from numpy import ma


# define dictionary data structures
def newWM(name):
    return {'name':name,'Hs_10min':[],'Hs':[],'Tp_10min':[],'Tp':[],'Tm02_10min':[],'Tm02':[],
            'DDM_10min':[],'DDM':[],'Tm01_10min':[],'Tm01':[],'DDP_10min':[],'DDP':[]}

def newWL(name):
    return {'name':name, 'Hlat':[]}
  
def newWI(name):
    return {'name':name,'FF_10min':[],'DD_10min':[],'FF':[],'DD':[]}

def tryfloat(string):
    try:
        fl = sp.float32(string.strip())
    except:
        fl = sp.nan
    return fl

def read_d22(station,start=None,end=None):
#if True:
    rig = station
    if start==None:
        start = dt.datetime.now() - dt.timedelta(days=1)
        start = start.strftime("%Y%m%d")
    if end==None:
        end = dt.datetime.now()
        end = end.strftime("%Y%m%d")
    
    # convert to datetime object if neccessary:
    st=start
    en=end

    try:
        st=dt.datetime(int(start[0:4]),int(start[4:6]),int(start[6:8]))
        en=dt.datetime(int(end[0:4]),int(end[4:6]),int(end[6:8]))
    except TypeError:
        print(' ')

####################################
# Importing and processing data
####################################

    searchlines=[]

# Read all lines in fil and append to searchlines
    for d in range(int(pl.date2num(st)),int(pl.date2num(en))+1): 
        try:
            dy=pl.num2date(d).strftime("%Y%m%d")
            #print('try to read d22 from /opdata...')
            f = open("/vol/data/offshore/"+rig+"/d22/"+dy+".d22", "r")
        except IOError:
            dy=pl.num2date(d).strftime("%Y/%Y%m%d")
            try:
                 #print('try to read d22 from hindcast...')
                 f = open("/lustre/storeB/immutable/short-term-archive/DNMI_OFFSHORE/"+rig+"/d22/"+dy+".d22", "r")       
            except IOError:
                #print('try to read d22 from starc')
                try:
                    f = open("/lustre/WATZMANN/storeB/immutable/short-term-archive/DNMI_OFFSHORE/"+rig+"/d22/"+dy+".d22", "r")
                except IOError:
                    print('no d22 file for station '+station+' at '+dy)
                    continue
        searchlines = searchlines + f.readlines()
        f.close()

#Creat dictionaries and variables
    tseries=[]
    dat = {'10min':[],'1hr':[]}
    WMlist = [newWM(name) for name in ['WM1', 'WM2', 'WM3']]
    #WLlist = [newWL(name) for name in ['WL1', 'WL2', 'WL3']]
    WIlist = [newWI(name) for name in ['WIA', 'WIB', 'WIC', 'WID', 'WIE']]

#Extract data of choice - reading searchlines
#    new = True
    for i, line in enumerate(searchlines):
        if "!!!!" in line: 
            tseriesl = []
            for l in searchlines[i+3:i+5]:
                tseriesl.append(l.strip())
            tseries.append(' '.join(tseriesl))
            date_object = dt.datetime.strptime(' '.join(tseriesl),'%d-%m-%Y %H:%M')
            dat['10min'].append(date_object)
#            for W in WMlist+WIlist:
#                W['flag']=False
            for WM in WMlist:
                Hs,Tp,Tm02,Tm01,DDP,DDM = sp.nan, sp.nan,sp.nan, sp.nan, sp.nan,sp.nan
                WM['Hs_10min'].append(Hs)
                WM['Tp_10min'].append(Tp)
                WM['Tm02_10min'].append(Tm02)
                WM['Tm01_10min'].append(Tm01)
                WM['DDP_10min'].append(DDP)
                WM['DDM_10min'].append(DDM)
#                WM['flag']=True
            for WI in WIlist:
                FF,DD = sp.nan, sp.nan
                WI['FF_10min'].append(FF)
                WI['DD_10min'].append(DD)
        # look for values in the subsequent lines:       
        #for j, linej in enumerate(searchlines[i+1:])
        for WM in WMlist:
            if str(WM['name']) in line: # if values are found, replace the NaNs with the values
                try:
                    Hs   = tryfloat(searchlines[i+3-1])
                    Tp   = tryfloat(searchlines[i+6-1])
                    Tm02   = tryfloat(searchlines[i+12-1]) #Tm02
                    Tm01 = tryfloat(searchlines[i+13-1])
                    DDP  = tryfloat(searchlines[i+20-1])
                    DDM  = tryfloat(searchlines[i+21-1])
                    WM['Hs_10min'][-1] = Hs
                    WM['Tp_10min'][-1] = Tp
                    WM['Tm02_10min'][-1] = Tm02
                    WM['Tm01_10min'][-1] = Tm01 
                    WM['DDP_10min'][-1] = DDP
                    WM['DDM_10min'][-1] = DDM
                except IndexError:
                    continue

        for WI in WIlist:
            if str(WI['name']) in line:
                try:
                    FF=tryfloat(searchlines[i+10])
                    DD=tryfloat(searchlines[i+13])
                    WI['FF_10min'][-1] = FF
                    WI['DD_10min'][-1] = DD
                except IndexError:
                    continue
#                WI['flag']=True
#                if "!!!!" in linej: #continue
#                    break

#Convert data to arrays
    dat['10min']=sp.array(dat['10min'])
    for WM in WMlist + WIlist: #+ WLlist:
      for var in WM.keys():
        if var != 'name':
            WM[var] = sp.array(WM[var])
            #Set all negative values to nan
            WM[var][WM[var]<0] = sp.nan
            #If variable does not exist
            if len(sp.array(WM[var]))==0:
                WM[var] = np.empty(len(dat['10min']))
                WM[var][:] = sp.nan  
#Create 1 hour averages (max for Tp)
            if var == 'Hs_10min':
                hs = WM['Hs_10min']
                hs[hs>30.] = sp.nan # set maximum hs to 30m
                hs[hs <= 0.] = sp.nan
                WM['Hs'] = sp.sqrt(hourmean(sp.power(hs,2))) # use the wave energy for averaging!
            if var == 'Tm02_10min':
                tm02 = WM['Tm02_10min']
                tm02[tm02>20] = sp.nan
                tm02[tm02<4] = sp.nan
                WM['Tm02'] = sp.sqrt(hourmean(sp.power(tm02,2)))
            if var == 'Tm01_10min':
                WM['Tm01'] = sp.sqrt(hourmean(sp.power(WM['Tm01_10min'],2)))
            if var == 'DDP_10min': # use every 6th direction measurement
                WM['DDP'] = WM['DDP_10min'][::6]
            if var == 'DDM_10min':
                WM['DDM'] = WM['DDM_10min'][::6]
            if var == 'FF_10min':
                u,v = UV(WM['DD_10min'],WM['FF_10min'],met=True)
                um,vm = hourmean(u), hourmean(v)
                WM['DD'],WM['FF'] = DD_FF(-um,-vm)
            if var == 'Tp_10min':
                tp = WM['Tp_10min']
                tp[tp>30] = sp.nan
                WM['Tp'] = hourmax(tp)

    dat['1hr']=dat['10min'][::6]
    return dat, WMlist, WIlist

 
def hourmean(varin):
    ''' average function that uses an average defined by d22mean
    and some despiking before the averaging '''
    var = despike(varin)
    var = ndimage.generic_filter(var, d22mean, size=6, origin=0,mode='constant')
    return var[::6]

def d22mean(sample):
    ''' return the mean if at least 3 of 6  10min values are finite
     Returns NaN only if more than 3 values missing, or if all are the same (std=0)
     '''
    if (sum(sp.isfinite(sample)) > 2 and sp.std(sample) > 0.):
        average = sp.mean(sample[sp.isfinite(sample)])
    else:
        average = sp.nan
    return average

def hourmax(varin):
    ''' average function that uses an average defined by d22mean
    and some despiking before the averaging '''
    #var = despike(varin)
    var = ndimage.generic_filter(varin, max, size=6, origin=0,mode='constant')
    return var[::6]


def despike(samplein):
    '''
    replace spikes with the median of sourounding 3 values if they are more than 50% (or only a third) of the median
    '''
    sample=samplein.copy()
    med = ndimage.median_filter(sample,size=3,mode='mirror')
    mask = sp.logical_or(sample/med > 1.5, sample/med < 0.66)
    sample[mask] = med[mask]
    return sample

def UV(DD,FF,met=False):
    '''
    #requires testing
    DD is the heading direction of the wind
    to get met. standart, call UV(DD,FF,met=True)
    '''
    u = FF * sp.sin(DD*sp.pi/180)
    v = FF * sp.cos(DD*sp.pi/180)
    if met==True:
        u,v=-u,-v
    return u,v


def DD_FF(u,v):
    ''' calculates wind/current speed and direction from u and v components
    #
    if u and v are easterly and northerly components,
    DD is heading direction of the wind. 
    to get meteorological-standard, call DD_FF(-u, -v)
    '''
    DD = ma.arctan2(u, v)*180/sp.pi
    DD[DD < 0] = 360 + DD[DD < 0]
    FF = ma.sqrt(u**2 + v**2)
    return DD, FF


