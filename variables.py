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

MWAM4 = {'Hs':  ['significant wave height','m','hs'],
               'Tp':  ['peak period','s','tp'],
               'Tm02':['mean period','s','tm2'],
#               'DDP': ['peak direction','degree',203,'oceanographic'],
               'DDM': ['mean direction','degree','thq','oceanographic'],
               'Hs_s':['swell significant wave height','m','hs_swell'],
               'Tp_s':['swell peak period','s', 'tp_swell'],
               'Tm02_s':['swell mean period','s','tm2_swell'],
#               'DDP_s':['swell peak direction','degree',222,'oceanographic'],
               'DDM_s':['swell mean direction','degree','thq_swell','oceanographic'],
               'FF':['wind speed','m/s','ff'],
               'DD':['wind direction','degree','dd','oceanographic']} 

#MWAM10 = MWAM4
#MWAM8 = MWAM4
MWAM8 = {'Hs':  ['significant wave height','m','hs'],
               'Tp':  ['peak period','s','tp'],
               'Tm02':['mean period','s','tm2'],
#               'DDP': ['peak direction','degree',203,'oceanographic'],
               'DDM': ['mean direction','degree','thq','oceanographic'],
               'Hs_s':['swell significant wave height','m','hs_swell'],
#               'Tp_s':['swell peak period','s', 'tp_swell'],
#               'Tm02_s':['swell mean period','s','tm2_swell'],
#               'DDP_s':['swell peak direction','degree',222,'oceanographic'],
               'DDM_s':['swell mean direction','degree','thq_swell','oceanographic'],
               'FF':['wind speed','m/s','ff'],
               'DD':['wind direction','degree','dd','oceanographic']} 


EXP = {'Hs':  ['significant wave height','m','hs'],
               'Tp':  ['peak period','s','tp'],
               'Tm02':['mean period','s','tm2'],
#               'DDP': ['peak direction','degree',203,'oceanographic'],
               'DDM': ['mean direction','degree','thq','oceanographic'],
               'Hs_s':['swell significant wave height','m','hs_swell'],
               'Tp_s':['swell peak period','s', 'tp_swell'],
               'Tm02_s':['swell mean period','s','tm2_swell'],
#               'DDP_s':['swell peak direction','degree',222,'oceanographic'],
               'DDM_s':['swell mean direction','degree','thq_swell','oceanographic'],
               'FF':['wind speed','m/s','ff'],
               'DD':['wind direction','degree','dd','oceanographic']} 

ECWAM =    {'Hs':  ['significant wave height','m','significant_wave_height'],
               'Tp':  ['peak period','s','peak_wave_period'],
               'Tm02':['mean period','s','mean_wave_period'],
               'DDM': ['mean direction','degree','wave_direction','oceanographic'],
               'Hs_s':['swell significant wave height','m','significant_swell_wave_height'],
               'Tm02_s':['swell mean period','s','sea_surface_swell_wave_period'], # or is this Tp_s ?
               'DDM_s':['swell mean direction','degree','sea_surface_swell_wave_to_direction','oceanographic']}


LAWAM = ECWAM


AROME = {'FF': ['wind speed', 'm/s', ''],
         'DD': ['wind direction', 'degree', '', 'meteorological']}

WAMAROME2W = {'FF': ['wind speed', 'm/s', ''],
         'Hs': ['significant wave height', 'm', '']}

WAMAROME1W = {'FF': ['wind speed', 'm/s', ''],
         'Hs': ['significant wave height', 'm', '']}


HIRLAM8=AROME



