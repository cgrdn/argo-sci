
from pathlib import Path

import numpy as np
import pandas as pd
import argopandas as argo

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks', palette='colorblind')
import cmocean.cm as cmo

# define function to make cartopy box work for argopandas
def reshape_extent(ex):
    return [ex[2], ex[0], ex[3], ex[1]]

# oceanic reagions, north atlantic encapsulates the other two
na_box = [-68.551034, 3.02453, 40.409517, 80.0]
lab_sea = [-67.366619, -41.982122, 53.179934, 67.562026]
scotian_shelf = [-68.551034, -58.248323, 40.409517, 48.092161]

# variables to get, add QC and ADJUSTED
varnames = ['PRES', 'PSAL', 'TEMP', 'DOXY']
varnames = [[v, f'{v}_QC', f'{v}_ADJUSTED', f'{v}_ADJUSTED_QC'] for v in varnames]
varnames = [item for sublist in varnames for item in sublist]

# profiles, data
ix = argo.synthetic_prof.subset_rect(*reshape_extent(na_box)).subset_parameter('DOXY').subset_date('2022-01', '2023-01')
prof = ix.prof
ix = ix.set_index('file')
ix = ix.join(prof[['PLATFORM_NUMBER', 'CYCLE_NUMBER']])
df = ix.levels[varnames]
ix['LAST_REPORTED'] = [ix.loc[ix.PLATFORM_NUMBER == p].date.iloc[-1] for p in ix.PLATFORM_NUMBER]
ix['deployment_date'] = [argo.float(f).synthetic_prof.date.min() for f in ix.PLATFORM_NUMBER]
ix['deployment_year'] = [d.year for d in ix.deployment_date]
last_cycle = ix.loc[[True if c == ix.loc[ix.PLATFORM_NUMBER == p].CYCLE_NUMBER.max() else False for p in ix.PLATFORM_NUMBER.unique() for c in ix.loc[ix.PLATFORM_NUMBER == p].CYCLE_NUMBER]]
