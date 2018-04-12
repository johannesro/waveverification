#!/usr/bin/env python
# ./validate 201409
import scipy as sp
#import numpy.ma as ma
import matplotlib
matplotlib.use('Agg')
import pylab as pl
from netCDF4 import Dataset, MFDataset, MFTime, date2num, num2date
import dataanalysis as da
import os
import datetime as dt
import validationtools as vt
from stationlist import locations as locations
from stationlist import WMsensors, bestWMsensor
import sys
from collectdatatools import validationfile
#import calendar

print("The Python version is %s.%s.%s" % sys.version_info[:3])

interactive=False

if interactive:
    year=2017
    month=5
else:
    if len(sys.argv) > 1:
        year = int(sys.argv[1])
        month = int(sys.argv[2])
    else:
        now = dt.datetime.now()
        year = now.year
        month = now.month

t1 = dt.datetime(year,month,1)
if month==12:
    t2 = dt.datetime(year+1,1,1)
else:
    t2 = dt.datetime(year,month+1,1)

timestr = t1.strftime('%Y%m')

# plotpath
#ppath = '/vol/hindcast3/waveverification/'+timep+'/'
ppath = '/lustre/storeA/project/fou/om/waveverification/'+timestr+'/'

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

if interactive:
    s='ekofiskL'
    locations = {s: locations[s]}

for station, parameters in locations.iteritems():
    print(' ')
    print('verification of station '+station+' for '+timestr)

    # open file
    path = '/lustre/storeA/project/fou/om/waveverification/data'
    
    vf = validationfile(path,station,year,month)
    time = vf.time
    OBS  = vf.get_obs()

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
        print('make plots for parameter '+varname)
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

        #        
        # make scatter and qq plots for direction intervals
        #
        #for interval in [(0,90), (90,180), (180,270), (270,360)]:
        #    fig=pl.figure(figsize=[10,5])
        #    ax1=fig.add_subplot(121)
        #    ax2=fig.add_subplot(122)
        #    direction = OBS['DD'][0]
        #    mask = sp.logical_and(direction >= interval[0], direction < interval[1])
        #    for gname, var in modeldata.iteritems():
        #        if sp.isnan(var[0][mask]).all() or sp.isnan(obs[mask]).all():
        #            continue
        #        vt.scqqplot(obs[mask], var[0][mask],color=ct[gname],  label=gname, ax1=ax1, ax2=ax2, prob=sp.arange(0.02,1.,0.02))
        #    ax1.set_title(station+' '+varname+' ['+units+'] for wind direction interval '+str(interval))
        #    pfilename = station+'_'+varname+'_scatterqq_dir'+str(interval[0])+'.png'
        #    fig.tight_layout(pad=0.2)
        #    fig.savefig(os.path.join(ppath+varname+'_directional',pfilename))
        #    #pl.close()

        #
        # make time series plot (with all available observations)
        #
        fig = pl.figure(figsize=[12,8])
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212,sharex=ax1)
        pl.xlim([t1,t2])
        ax2.xaxis.set_minor_locator(minorLocator)
        OBSarray = OBS[varname][:]
        OBSarray[OBSarray>1000]=sp.nan
        for i, Nobs in enumerate(OBSarray.tolist()):
            obsnum = i+1
            if all(sp.isnan(Nobs)):
                continue
            mask = sp.isfinite(Nobs)
            ax1.plot(sp.array(time)[mask], sp.array(Nobs)[mask], '.', label='observation #'+str(obsnum),lw=1)
            #ax1.plot(time, Nobs, label='observation #'+str(obsnum),lw=1)
        for gname, var in modeldata.iteritems():
            ax1.plot(time, var[0],'-',color=ct[gname], label=gname, lw=1.5)
        ax1.legend(fontsize='small')
        ax1.grid('on',which='minor');ax1.grid('on',which='major',linestyle='--',linewidth=0.5)
        ax1.set_title(station+' '+varname+' ['+units+']')
        # put wind direction into the same panel!
        ax2.plot(time, OBS['DDP'][0],'.',label='peak wave direction (obs)')
        ax2.plot(time, OBS['DD'][0],'.',label='wind direction (obs)')
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
            reini = vf.nc.getncattr(gname+'_reinitialization_step')
            vt.forecastskillplot(obs, var[:], reini, vt.amerr, color=ct[gname],  label=gname, ax=ax1)
            vt.forecastskillplot(obs, var[:], reini, vt.rmse, color=ct[gname],  label=gname, ax=ax2)
            vt.forecastskillplot(-obs, -var[:], reini, vt.bias, color=ct[gname],  label=gname, ax=ax3)
            vt.forecastskillplot(obs, var[:], reini, vt.pearsonr, color=ct[gname],  label=gname, ax=ax4)
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
            modeldata = select_var_from_models(vf,varname)
            OBSarray = OBS[varname][:]
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
        
    vf.nc.close()
   
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
    DD = OBS['DD'][0]
    DDM = OBS['DDM'][0]
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

