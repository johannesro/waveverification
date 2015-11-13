
from METread import *
import rpy2.robjects as robjects # open an R session in python
import rpy2.robjects.numpy2ri as rpyn # convert R data to numpy objects
import sys
import os

os.system('export R_LIBS=/disk1/local/R-packages')
robjects.r('''library(miIO)''') # load the miIO package in the R session
miReadFelt = robjects.r('miReadFelt') # make the miReadFelt function available in python

feltnumbers = {'Hs': 200, 'Tp': 201, 'Tm02':202, 'DDP':203, 'DDM':204,  'Hs_s': 220, 'Tp_s': 221, 'Tm02_s':223, 'DDP_s':222,
               'DDM_s':224,'FF':298,'DD':299} # felt parameter number
rfeltname = {'Hs': "HS.L0", 'Tp': "TP.L0", 'Tm02':"TS.L0", 'DDP':"DPP.L0", 'DDM':"DDM.L0",'Hs_s':"HSSW.L0",'Tp_s':"TPSW.L0",'Tm02_s':"TMSW.L0",
             'DDP_s':"DDPSW.L0", 'DDM_s':"DDMSW.L0",'FF':"FF.L0",'DD':"DD.L0"} # name of variables in r2py interface


class Blackhole(object):
    def write(self, string):
        pass

stdout = sys.stdout

def rx(data, var):
    return rpyn.ri2numpy(data.rx2(var))

def WAM4_modrun(location, run, varnamelist, step=1):
    '''
    location: (lat, lon) touple or list
    run: datetime object of model run initialization time (should be 00 or 12 hours)
    '''
    # ensure correct formats
    location = list(location)
    # check file
    filename=run.strftime("/vol/hindcast3/waveverification/DNMI_WAVE/%Y/%m/%d/g4kmwave%H.dat_%Y%m%d")
    if not os.path.exists(filename):
        filename=run.strftime("/starc/DNMI_WAVE/%Y/%m/%d/g4kmwave%H.dat_%Y%m%d")
    print('reading '+filename)
    nan = sp.ones(66)*sp.nan
    if os.path.isfile(filename):
        # save location as r matrix
        pos = robjects.r.matrix(robjects.FloatVector(location),nrow=1,byrow=True) #matrix with first column: latitude, second column: longitude
        # save variable list as r matrix 
        # rows for each variable are [vertical coordinate, varnumber, level1, level2] given by column 8-11 of 'rfinh <felt_file>'
        varlist=[]
        for varname in varnamelist:
            varlist=varlist+[3,feltnumbers[varname],0,0]
        prm = robjects.r.matrix(robjects.IntVector(varlist),nrow=len(varnamelist),byrow=True) 
        # make r sequence for model lead times
        robjects.r('''prg <- seq(0, 66, by=1) ''')
        prg = robjects.r['prg']
        # use R function to read from felt file
        sys.stdout = Blackhole()
        x = miReadFelt(files=filename,sites=pos, prm=prm, prg=prg, collapse_time=False, acc=6, df=True);
        sys.stdout = stdout
        # get time as datetime object
        y,m,d,h = rx(x,'YEAR'),rx(x,'MONTH'),rx(x,'DAY'),rx(x,'HOUR')
        time = []
        for i in range(len(h)):
            try:
                time.append(dt.datetime(int(y[i]), int(m[i]), int(d[i]), int(h[i])))
            except ValueError: # some felt file have missing time steps. If so, fill time list with next hour
                try:
                    time.append(time[-1]+dt.timedelta(hours=1))
                except IndexError:
                    time.append(run) # if this is the first time step
        datadict = {'time':time}
        for varname in varnamelist:
            var = rx(x,rfeltname[varname])
            if not sp.array(var==robjects.rinterface.NULL).all():
                datadict.update({varname:var})
            else:
                datadict.update({varname:nan})
                print('variable '+varname+' not found in '+filename)
    else:
       print('file '+filename+' does not exist; returning missing values')
       datadict = {'time':nan}
       for name,rname in rfeltname.iteritems():
           datadict.update({name:nan})
    return datadict


def WAM10_modrun(location, run, varnamelist, step=1):
    '''
    location: (lat, lon) touple or list
    run: datetime object of model run initialization time (should be 00 or 12 hours)
    '''
    # ensure correct formats
    location = list(location)
    # check file
    filename=run.strftime("/vol/hindcast3/waveverification/DNMI_WAVE/%Y/%m/%d/g10kmwave%H.dat_%Y%m%d")
    if not os.path.exists(filename):
        filename=run.strftime("/starc/DNMI_WAVE/%Y/%m/%d/g10kmwave%H.dat_%Y%m%d")
    print('reading '+filename)
    nan = sp.ones(66)*sp.nan
    if os.path.isfile(filename):
        # save location as r matrix
        pos = robjects.r.matrix(robjects.FloatVector(location),nrow=1,byrow=True) #matrix with first column: latitude, second column: longitude
        # save variable list as r matrix 
        # rows for each variable are [vertical coordinate, varnumber, level1, level2] given by column 8-11 of 'rfinh <felt_file>'
        varlist=[]
        for varname in varnamelist:
            varlist=varlist+[3,feltnumbers[varname],0,0]
        
        prm = robjects.r.matrix(robjects.IntVector(varlist),nrow=len(varnamelist),byrow=True) 
        # make r sequence for model lead times
        robjects.r('''prg <- seq(0, 66, by=1) ''')
        prg = robjects.r['prg']
        # use R function to read from felt file
        sys.stdout = Blackhole()
        x = miReadFelt(files=filename,sites=pos, prm=prm, prg=prg, collapse_time=False, acc=6, df=True);
        sys.stdout = stdout
        # get time as datetime object
        y,m,d,h = rx(x,'YEAR'),rx(x,'MONTH'),rx(x,'DAY'),rx(x,'HOUR')
        time = []
        for i in range(len(h)):
            try:
                time.append(dt.datetime(int(y[i]), int(m[i]), int(d[i]), int(h[i])))
            except ValueError: # some felt file have missing time steps. If so, fill time list with next hour
                time.append(time[-1]+dt.timedelta(hours=1))
        datadict = {'time':time}
        for varname in varnamelist:
            var = rx(x,rfeltname[varname])
            if not sp.array(var==robjects.rinterface.NULL).all():
                datadict.update({varname:var})
            else:
                datadict.update({varname:nan})
                print('variable '+varname+' not found in '+filename)
    else:
       print('file '+filename+' does not exist; returning missing values')
       datadict = {'time':nan}
       for name,rname in rfeltname.iteritems():
           datadict.update({name:nan})
    return datadict



def HIRLAM8_modrun(location, run, varnamelist, step=1):
    '''
    location: (lat, lon) touple or list
    run: datetime object of model run initialization time (should be 00, 06, 12 or 18 hours)
    '''
    # ensure correct formats
    location = list(location)
    # check file
    filename=run.strftime("/starc/DNMI_HIRLAM8/%Y/%m/%d/h8km%H.dat_%Y%m%d")
    print('reading '+filename)
    if os.path.isfile(filename):
        # save location as r matrix
        pos = robjects.r.matrix(robjects.FloatVector(location),nrow=1,byrow=True) #matrix with first column: latitude, second column: longitude
        # save variable list as r matrix 
        # rows for each variable are [vertical coordinate, varnumber, level1, level2] given by column 8-11 of 'rfinh <felt_file>'
        prm = robjects.r.matrix(robjects.IntVector([2,33,1000,0]+[2,34,1000,0]),nrow=2,byrow=True) # more variables:,[2,31,1000,0],[2,17,1000,0])
        # make r sequence for model lead times
        robjects.r('''prg <- seq(0, 66, by=1) ''')
        prg = robjects.r['prg']
        # use R function to read from felt file
        sys.stdout = Blackhole()
        x = miReadFelt(files=filename,sites=pos, prm=prm, prg=prg, collapse_time=False, acc=6, df=True);
        sys.stdout = stdout
        FF = rx(x,"FF") # Put the data into a python array using ri2numpy
        DD = rx(x,"DD")
        # get time as datetime object
        y,m,d,h = rx(x,'YEAR'),rx(x,'MONTH'),rx(x,'DAY'),rx(x,'HOUR')
        #time = [dt.datetime(int(y[i]), int(m[i]), int(d[i]), int(h[i])) for i in range(len(h))]
        time = []
        for i in range(len(h)):
            try:
                time.append(dt.datetime(int(y[i]), int(m[i]), int(d[i]), int(h[i])))
            except ValueError: # some felt file have missing time steps. If so, fill time list with next hour
                time.append(time[-1]+dt.timedelta(hours=1))
        datadict = {'time':time}
        datadict.update({'FF':FF, 'DD':DD})
    else:
       print('file '+filename+' does not exist; returning missing values')
       nan = sp.ones(66)*sp.nan
       datadict = {'time':nan,'FF':nan, 'DD':nan}
    return datadict




