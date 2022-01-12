#!/usr/bin/python

import numpy as np
import pandas as pd
import argopandas as argo

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
sns.set(style='ticks', palette='colorblind')
import cmocean.cm as cmo

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# select lab sea region, doxy floats, last 5 years, ascending
# note - just do 2021 to start
lab_sea = [-67, -43, 55, 62.5]
bx = argo.synthetic_prof.subset_parameter('DOXY')\
    .subset_rect(lab_sea[2], lab_sea[0], lab_sea[3], lab_sea[1])\
    .subset_date('2021-01')\
    .subset_direction('asc')

# this float is bad
bx = bx[~bx.file.str.contains('4901779')]

# load in the data - this will take a minute
df = bx.levels[[
    'PRES', 'PRES_QC', 
    'TEMP', 'TEMP_QC', 
    'PSAL', 'PSAL_QC', 
    'DOXY', 'DOXY_QC'
]]

# get surface (z <= 10m) data
surf = df[df.PRES <= 10]

# average surface value for each profile
bx['surf_doxy'] = [surf.loc[f].DOXY.mean() for f in bx.file]

# define grid edges
boxsize = 0.5
xgrid = np.arange(lab_sea[0], lab_sea[1]+boxsize, boxsize)
ygrid = np.arange(lab_sea[2], lab_sea[3]+boxsize, boxsize)

X, Y = np.meshgrid(xgrid, ygrid)
N, M = X.shape

O2 = np.nan*np.ones(X.shape)

for i in range(N-1):
    for j in range(M-1):
        y1, y2 = ygrid[i], ygrid[i+1]
        x1, x2 = xgrid[j], xgrid[j+1]
        O2[i,j] = bx[
            (bx.longitude > x1) & (bx.longitude < x2) &\
            (bx.latitude > y1) & (bx.latitude < y2)
        ].surf_doxy.mean()
