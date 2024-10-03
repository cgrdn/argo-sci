
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(context='paper', style='ticks', palette='colorblind')

# import argopy
# argopy.set_options(mode='expert')

# fetcher = argopy.DataFetcher()
# index = argopy.ArgoIndex().load()

import argopandas as argo

# select lab sea region, doxy floats, last 5 years, ascending
lab_sea = [-67, -43, 55, 62.5]
bx = argo.synthetic_prof.subset_parameter('DOXY')\
    .subset_rect(lab_sea[2], lab_sea[0], lab_sea[3], lab_sea[1])\
    .subset_date('2016-01')\
    .subset_direction('asc')

# load in the data - this will take a minute
df = bx.levels[[
    'PRES', 'PRES_QC', 
    'TEMP', 'TEMP_QC', 
    'PSAL', 'PSAL_QC', 
    'DOXY', 'DOXY_QC'
]]

prof = bx.prof[['JULD', 'LATITUDE', 'LONGITUDE']]
df = df.join(prof, on=['file', 'N_PROF'])
df['DATE'] = [pd.Timestamp(j -  365.25*20, unit='D') for j in df.JULD]
df['YEAR'] = [t.year for t in df.DATE]
df['YEAR_DAY'] = [t.day_of_year for t in df.DATE]

allowed_flags = [b'1', b'2', b'5']
df = df.loc[df.TEMP_QC.isin(allowed_flags)]
df = df.loc[df.PSAL_QC.isin(allowed_flags)]