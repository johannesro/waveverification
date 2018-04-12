#!/usr/bin/env python2.7

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
from datetime import datetime, timedelta
from dateutil import tz
from matplotlib.dates import DayLocator, HourLocator, DateFormatter
import os
import METread

datapath = '/lustre/storeA/project/fou/om/waveverification/'
#os.chdir(datapath)
print(os.getcwd())

#### Read Hs and Tz 
#y, m, d, H, M, Hs, Tz = np.loadtxt('/home/kaihc/DATAWELL/ASCII/waveriderdata.txt', skiprows=1, unpack=True)
today = datetime.today()
url = '157.249.178.22/datawell/waved/Tennholmen/%i/%2.2i/' % (today.year, today.month)
filename = 'Tennholmen{0x324}%i-%2.2i.csv' % (today.year, today.month)
#location = (68.21 , 14.48 )
location = (67.3406, 13.5618)
os.system('rm ' + filename)
os.system('wget ' + url + filename)
timeseconds, Hs, Tz = np.loadtxt(filename, skiprows=1, usecols=(2,3,4), unpack=True)

# UTC/local time conversion
HERE = tz.tzlocal()
UTC = tz.gettz('UTC')

# Make datetime objects from data
#dte = [datetime(yi[i], mi[i], di[i], Hi[i], Mi[i]) for i in range(0,len(yi))]
reftime = datetime(1970,1,1)
dte = [reftime + timedelta(seconds=timeseconds[i]) for i in range(len(timeseconds))]

# Convert from UTC to local time
dteUTC = [dte[i].replace(tzinfo=UTC) for i in range(0,len(dte))]
dteLOCAL = [dteUTC[i].astimezone(HERE) for i in range(0,len(dteUTC))]

zoneName = time.tzname[time.daylight]

# preperations for plotting
halfdays = HourLocator(byhour=[0, 12])
hours = HourLocator(byhour=[3, 6, 9, 15, 18, 21])
dayfmt = DateFormatter('%H\n%b %d')
hourfmt = DateFormatter('%H')

#### Choose only the last N entries
N = 12
HsPlot = Hs[-N:]
TzPlot = Tz[-N:]
dteLOCALplot = dteLOCAL[-N:]


### fetch model data
modelinterval = 6
modelrun = today - timedelta(hours=today.hour % modelinterval) 
wam = METread.MWAM4_modrun(location, modelrun, ['Hs', 'Tp'])
if np.isnan(wam['Hs']).all():
    modelrun = today - (timedelta(hours=today.hour % modelinterval)) - timedelta(hours=modelinterval) 
    wam = METread.MWAM4_modrun(location, modelrun, ['Hs', 'Tp'])

print('reading model data complete')

wamtimeUTC  = [wam['time'][i].replace(tzinfo=UTC) for i in range(len(wam['time']))]
wamtimeLOCAL = [wamtimeUTC[i].astimezone(HERE) for i in range(len(wam['time']))]

#### Plot Hs and Ts in same figure
fig = plt.figure()
print('initiated figure')

ax1 = fig.add_subplot(111)
ax1.plot(dteLOCALplot, HsPlot, 'bd-', linewidth = 3, label='waverider')
ax1.plot(wamtimeLOCAL, wam['Hs'], 'b--', linewidth = 3, label='forecast')
plt.legend(loc='lower right')

xlabelStr = "Local time"
xlab = ax1.set_xlabel(xlabelStr, fontsize = 14)

# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel(r'Significant wave height $H_s$ [m]', fontsize = 14, color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')

ax2 = ax1.twinx()
ax2.plot(dteLOCALplot, TzPlot, 'gd-', linewidth = 3)
ax2.plot(wamtimeLOCAL, wam['Tm02'], 'g--', linewidth = 3)

ax2.set_ylabel(r'Mean wave period $T_{m02}$ [s]', fontsize = 14, color='g')

# set axes limits
ax1.set_ylim([0, 1.3*np.max([np.max(wam['Hs'])  ,np.max(HsPlot)]) ])
ax2.set_ylim([0, 1.1*np.max([np.max(wam['Tm02']),np.max(TzPlot)]) ])
ax1.set_xlim([today-timedelta(hours=12),today+timedelta(hours=12)])

for tl in ax2.get_yticklabels():
    tl.set_color('g')


titlestr = "Vestfjorden waverider"

# Plot every third hour with minor ticks, but remove at midnight

ax1.xaxis.set_minor_locator(hours)
ax1.xaxis.set_minor_formatter(hourfmt)

ax1.xaxis.set_major_locator(halfdays)
ax1.xaxis.set_major_formatter(dayfmt)

plt.title(titlestr, fontsize = 20)

#fig.autofmt_xdate()
ax1.grid('on', axis='x', which='both')

print('save figure Tennholmen_ts.png')
plt.savefig('Tennholmen_ts.png', bbox_extra_artists=(xlab,), bbox_inches='tight')

#os.system('rsync -u Tennholmen_ts.png projects.met.no:/var/www/waverider/.')
os.system('scp Tennholmen_ts.png projects.met.no:/var/www/waverider/.')

#plt.show()

