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
from stationlist import locations, WMsensors, bestWMsensor
import sys
#import calendar

interactive=False

if interactive:
    timep='201411'
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
ppath = '/disk4/waveverification/'+timep+'/'
# set color table for models
ct = {'Subjective': 'r', 'WAM10': 'c', 'WAM4':'m', 'ECWAM':'k', 'LAWAM':'0.25', 'AROME': 'g', 'HIRLAM8': 'y', 'MWAM10':'b' }

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
                if (G[gname].variables[varname].Convention=='oceanographic'):
                    var=var+180
                    var[var>360.]=var[var>360.]-360.
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

if interactive:
    s='draugen'
    locations = {s: locations[s]}

for station, parameters in locations.iteritems():
    print ' '
    print 'verification of station '+station+' for '+timep
#
# open file
    path = '/disk4/waveverification/data'
#filename = starttime.strftime(station+'_201312.nc')
    filename = station+'_'+timep+'.nc'
    nc       = Dataset(os.path.join(path,filename),mode='r')
    time = num2date(nc.variables['time'],nc.variables['time'].units)
    G = nc.groups
    OBS  = G['OBS_d22']

    os.system('mkdir -p '+ppath+'Hs')
    os.system('mkdir -p '+ppath+'Hs_directional')
    os.system('mkdir -p '+ppath+'Tp')
    os.system('mkdir -p '+ppath+'Tp_directional')
    os.system('mkdir -p '+ppath+'Tm02')
    os.system('mkdir -p '+ppath+'Tm02_directional')
    os.system('mkdir -p '+ppath+'FF')
    os.system('mkdir -p '+ppath+'FF_directional')
    os.system('mkdir -p '+ppath+'DD')
    os.system('mkdir -p '+ppath+'DDM')
    os.system('mkdir -p '+ppath+'DDP')
    os.system('mkdir -p '+ppath+'directions')



# Specify which WM sensor to use for validation
    try:
        sensor = bestWMsensor[station]
    except KeyError:
        sensor = 0

#
# Perform validation for each variable
#
    for varname in ['Hs','Tm02','FF','Tp']:
        print 'make plots for parameter '+varname
        obs = ma.array(OBS.variables[varname][sensor])
        obs.data[obs.mask==True] = sp.nan # make sure all masked values are nan 
        obs.mask = sp.logical_or(obs.mask, sp.isnan(obs.data))
        units = OBS.variables[varname].units
# some quick outlier detection:
#        obs[obs>500] = sp.nan
#        if varname=='Hs':
#            obs[obs>30] = sp.nan

        

        if (all(sp.isnan(obs.data)) or all(obs.mask==True)):
            print 'no data for '+station+' during '+timestr
            continue

        # select variable from  each model:
        modeldata = select_var_from_models(G,varname)
        #        

        # make scatter and qq plot
        #
        fig=pl.figure(figsize=[10,5])
        ax1=fig.add_subplot(121)
        ax2=fig.add_subplot(122)
        for gname, var in modeldata.iteritems():
           if sp.isnan(var[0]).all() or sp.isnan(obs).all():
               continue
           vt.scqqplot(obs, var[0],color=ct[gname],  label=gname, ax1=ax1, ax2=ax2)
        ax1.legend(loc='lower right',fontsize='small')
        ax1.set_title(station+' '+varname+' ['+units+']'+' obs#'+str(sensor+1)+' '+timestr)
        pfilename = station+'_'+varname+'_scatterqq.png'
        fig.tight_layout(pad=0.2)
        fig.savefig(os.path.join(ppath+varname,pfilename))
        #pl.close()
        #        
        # make scatter and qq plots for direction intervals
        #
        for interval in [(0,90), (90,180), (180,270), (270,360)]:
            fig=pl.figure(figsize=[10,5])
            ax1=fig.add_subplot(121)
            ax2=fig.add_subplot(122)
            direction = OBS.variables['DD'][0]
            mask = sp.logical_and(direction >= interval[0], direction < interval[1])
            for gname, var in modeldata.iteritems():
                if sp.isnan(var[0][mask]).all() or sp.isnan(obs[mask]).all():
                    continue
                vt.scqqplot(obs[mask], var[0][mask],color=ct[gname],  label=gname, ax1=ax1, ax2=ax2, prob=sp.arange(0.02,1.,0.02))
            ax1.set_title(station+' '+varname+' ['+units+'] for wind direction interval '+str(interval))
            pfilename = station+'_'+varname+'_scatterqq_dir'+str(interval[0])+'.png'
            fig.tight_layout(pad=0.2)
            fig.savefig(os.path.join(ppath+varname+'_directional',pfilename))
            #pl.close()

        #
        # make time series plot (with all available observations)
        #
        fig = pl.figure(figsize=[12,8])
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212,sharex=ax1)
        pl.xlim([t1,t2])
        ax2.xaxis.set_minor_locator(minorLocator)
        OBSarray = OBS.variables[varname][:]
        OBSarray[OBSarray>1000]=sp.nan
        for i, Nobs in enumerate(OBSarray.tolist()):
            obsnum = i+1
            if all(sp.isnan(Nobs)):
                continue
            ax1.plot(time, Nobs, label='observation #'+str(obsnum),lw=1)
        for gname, var in modeldata.iteritems():
            ax1.plot(time, var[0],'--',color=ct[gname], label=gname, lw=2)
        ax1.legend(fontsize='small')
        ax1.grid('on',which='minor');ax1.grid('on',which='major',linestyle='--',linewidth=0.5)
        ax1.set_title(station+' '+varname+' ['+units+']')
        # put wind direction into the same panel!
        ax2.plot(time, OBS.variables['DDP'][0],'.',label='peak wave direction (obs)')
        ax2.plot(time, OBS.variables['DD'][0],'.',label='wind direction (obs)')
        ax2.legend(fontsize='small')
        ax2.grid('on',which='minor');ax2.grid('on',which='major',linestyle='--',linewidth=0.5)
        ax1.xaxis.set_major_locator(majorLocator)
        ax1.xaxis.set_minor_locator(minorLocator)
        ax1.xaxis.set_major_formatter(fmt)
        pl.ylim([0,360])
        pl.yticks([90,180,270,360],['E','S','W','N'])
        pfilename = station+'_'+varname+'_tseries.png'
        fig.tight_layout(pad=0.2)

        fig.savefig(os.path.join(ppath+varname,pfilename))
        #pl.close()
#
# plot statistics as function of forcast time
#
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
        
#
# compare directions
#
    #
    # make time series plot of directions
    #
    fig = pl.figure(figsize=[12,8])
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312,sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(313,sharex=ax1, sharey=ax1)
    pl.xlim([t1,t2])
    ax3.xaxis.set_minor_locator(minorLocator)
    for varname, ax in {'DD':ax1, 'DDM':ax2, 'DDP':ax3}.iteritems():
        modeldata = select_var_from_models(G,varname)
        OBSarray = OBS.variables[varname][:]
        OBSarray[OBSarray>1000]=sp.nan
        for i, Nobs in enumerate(OBSarray.tolist()):
            obsnum = i+1
            if all(sp.isnan(Nobs)):
                continue
            ax.plot(time, Nobs,'.', label='observation #'+str(obsnum),lw=1)
        for gname, var in modeldata.iteritems():
            ax.plot(time, var[0],'--', color=ct[gname], label=gname, lw=2)
    ax1.legend(fontsize='x-small')
    ax2.legend(fontsize='x-small')
    ax3.legend(fontsize='x-small')
    ax1.grid('on',which='minor');ax1.grid('on',which='major',linestyle='--',linewidth=0.5)
    ax2.grid('on',which='minor');ax2.grid('on',which='major',linestyle='--',linewidth=0.5)
    ax3.grid('on',which='minor');ax2.grid('on',which='major',linestyle='--',linewidth=0.5)
    ax1.set_title(station+'\n wind direction DD [degree, met.]')
    ax2.set_title('mean wave direction DDM [degree, met.]')
    ax3.set_title('peak wave direction DDP [degree, met.]')
    ax1.xaxis.set_major_locator(majorLocator)
    ax1.xaxis.set_minor_locator(minorLocator)
    ax1.xaxis.set_major_formatter(fmt)
    pl.ylim([0,360])
    pl.yticks([90,180,270,360],['E','S','W','N'])
    pfilename = station+'_directions_tseries.png'
    fig.tight_layout(pad=0.2)
    fig.savefig(os.path.join(ppath+'directions',pfilename))


    if not interactive:
        pl.close('all')
        nc.close()



if interactive:
    pl.show()
#
# compute statistics
#



#osolete:
'''

# 
# Direction vs. wind speed scatter plot with color code for model error
#
    Hs_modeldata = select_var_from_models(G,'Hs')
    DDM_modeldata = select_var_from_models(G,'DDM')
    DD = OBS.variables['DD'][0]
    DDM = OBS.variables['DDM'][0]
    FF = OBS.variables['FF'][0]
    Hs = OBS.variables['Hs'][0]
    for gname, Hs_mod in Hs_modeldata.iteritems():
        DDM_mod = DDM_modeldata[gname]
        dHs = Hs_mod[0] - Hs
        dDDM = DDM_mod[0] - DDM
        fig = pl.figure(figsize=[10,5])
        #
        ax1 = fig.add_subplot(121)
        ax1.set_title('wave height (Hs) error in '+gname)
        pl.yticks([90,180,270,360],['E','S','W','N'])
        sc = ax1.scatter(FF, DD ,c=dHs ,cmap=pl.cm.seismic,vmin=-3, vmax=3)
        CB = pl.colorbar(sc)
        CB.set_label('model-obs: Hs [m]')
        #
        ax2 = fig.add_subplot(122)
        ax2.set_title('mean wave direction (DDM) error in '+gname)
        pl.yticks([90,180,270,360],['E','S','W','N'])
        #pl.xticks([90,180,270,360],['E','S','W','N'])
        sc = ax2.scatter(Hs, DD, c=dDDM , cmap=pl.cm.BrBG,vmin=-45,vmax=45)
        CB = pl.colorbar(sc)
        CB.set_label('model-obs: DDM [degree]')
        #
        ax1.set_ylim([0,360])
        ax1.set_ylabel('wind direction DD [degree]')
        ax1.set_xlabel('FF [m/s]')
        ax2.legend()
        ax2.set_ylim([0,360])
        #ax2.set_xlim([0,360])
        ax2.set_ylabel('wind directino DD [degree]')
        ax2.set_xlabel('wave height Hs [m]')
        pfilename = station+'_'+gname+'_FFvsDD_error.png'
        fig.tight_layout(pad=0.5)
        fig.savefig(os.path.join(ppath+'directions',pfilename))

'''

