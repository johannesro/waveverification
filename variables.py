d22 = {'Hs':  ['significant wave height','m',200],
               'Tp':  ['peak period','s',201],
               'Tm02':['mean period','s',202],
               'DDP': ['peak direction','degree',203, 'meteorological'],
               'DDM': ['mean direction','degree',204, 'meteorological'],
               'FF':['wind speed','m/s',298],
               'DD':['wind direction','degree',299, 'meteorological']} 

WAM4 = {'Hs':  ['significant wave height','m',200],
               'Tp':  ['peak period','s',201],
               'Tm02':['mean period','s',202],
               'DDP': ['peak direction','degree',203,'oceanographic'],
               'DDM': ['mean direction','degree',204,'oceanographic'],
               'Hs_s':['swell significant wave height','m',220],
               'Tp_s':['swell peak period','s', 221],
               'Tm02_s':['swell mean period','s',223],
               'DDP_s':['swell peak direction','degree',222,'oceanographic'],
               'DDM_s':['swell mean direction','degree',224,'oceanographic'],
               'FF':['wind speed','m/s',298],
               'DD':['wind direction','degree',299,'oceanographic']} 

WAM10 = WAM4

ECWAM =    {'Hs':  ['significant wave height','m','swh'],
               'Tp':  ['peak period','s','pp1d'],
               'Tm02':['mean period','s','mp2'],
               'DDM': ['mean direction','degree','mwd','oceanographic'],
               'Hs_s':['swell significant wave height','m','shts'],
               'Tm02_s':['swell mean period','s','p2ps'],
               'DDM_s':['swell mean direction','degree','mdts','oceanographic'],
               'FF':['wind speed','m/s','wind']}

LAWAM = ECWAM


AROME = {'FF': ['wind speed', 'm/s', ''],
         'DD': ['wind direction', 'degree', '', 'meteorological']}

HIRLAM8=AROME
