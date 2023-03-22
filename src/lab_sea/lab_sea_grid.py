#!/usr/bin/python

from pathlib import Path

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
    .subset_date('2016-01', '2022-01')\
    .subset_direction('asc')
# bx = argo.prof.subset_rect(lab_sea[2], lab_sea[0], lab_sea[3], lab_sea[1])\
#     .subset_date('2016-01', '2022-01')\
#     .subset_direction('asc')

# this float is bad
bx = bx[~bx.file.str.contains('4901779')]

# load in the data - this will take a minute
df = bx.levels[[
    'PRES', 'PRES_QC', 
    'TEMP', 'TEMP_QC', 
    'PSAL', 'PSAL_QC', 
    'DOXY', 'DOXY_QC'
]]

# average surface value for each profile
# bx['surf_doxy'] = [df.loc[f].TEMP[df.loc[f].PRES <= 10].mean() for f in bx.file]
bx['surf_doxy'] = [df.loc[f].DOXY[df.loc[f].PRES <= 10].mean() for f in bx.file]
# bx['juld'] = bx.prof[['JULD']].values
# bx['pandas_date'] = [pd.Timestamp(j -  365.25*20, unit='D') for j in bx.juld]
bx['year'] = [t.year for t in bx.date]
allowed_flags = [b'1', b'2', b'3', b'5']

# define grid edges
boxsize = 2
xgrid = np.arange(lab_sea[0], lab_sea[1]+boxsize, boxsize)
ygrid = np.arange(lab_sea[2], lab_sea[3]+boxsize, boxsize)
X, Y = np.meshgrid(xgrid, ygrid)
N, M = X.shape

# set up map, make figure
projection = ccrs.LambertConformal(central_latitude=55, central_longitude=-55)
transform = ccrs.PlateCarree()
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
extent  = [-67, -42, 54.5, 63]

'''
-------------------------------------------------------------------------------
INDIVIDUAL YEAR FIGURES
-------------------------------------------------------------------------------
'''
years = np.sort(bx.year.unique())
nyear = years.shape[0]
maps = []
O2 = np.nan*np.ones((nyear, N-1, M-1))
HG = np.nan*np.ones((nyear, N-1, M-1))
for year in years:
    fig = plt.figure()
    axes = [fig.add_subplot(211, projection=projection), fig.add_subplot(212, projection=projection)]
    year_bx = bx.subset_date(f'{year}-01', f'{year+1}-01')
    axes[0].set_title(f'{year}', loc='left', fontweight='bold')
    for i in range(N-1):
        for j in range(M-1):
            y1, y2 = ygrid[i], ygrid[i+1]
            x1, x2 = xgrid[j], xgrid[j+1]
            area = (year_bx.longitude > x1) & (year_bx.longitude < x2) &\
                (year_bx.latitude > y1) & (year_bx.latitude < y2)
            O2[year - years[0],i,j] = year_bx[area].surf_doxy.mean()
    ct = np.histogram2d(year_bx.longitude, year_bx.latitude, bins=[xgrid, ygrid])
    ct[0][ct[0] == 0] = np.nan
    HG[year - years[0], :, :] = ct[0].T
    # for ax, data, cm, vmin, vmax, lab in zip(axes, [O2[year - years[0],:,:], ct[0].T], [cmo.thermal, cmo.amp], [-1, 0], [8, 60], [f'T ({chr(176)}C)', '$N_{obs}$']):
    for ax, data, cm, vmin, vmax, lab in zip(axes, [O2[year - years[0],:,:], ct[0].T], [cmo.dense, cmo.amp], [250, 0], [350, 60], ['', '$N_{obs}$']):
        # plot map
        ax.set_extent(extent)
        ax.add_feature(cfeature.GSHHSFeature('low', 
            edgecolor='black', facecolor=cfeature.COLORS['land']))
        ax.patch.set_facecolor('lightgrey')
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='white')
        gl.xlocator = mticker.FixedLocator([-65, -45])
        gl.ylocator = mticker.FixedLocator([55, 60])
        gl.right_labels = False
        gl.left_labels = True
        gl.top_labels = False
        gl.xformatter = lon_formatter
        gl.yformatter = lat_formatter

        m = ax.pcolormesh(X, Y, data, cmap=cm, transform=transform, vmin=vmin, vmax=vmax)
        cb = plt.colorbar(m, ax=ax)
        cb.set_label(lab)
    fig.savefig(Path(f'../../figures/2022/gridded_surface_{boxsize}deg_{year}.png'), dpi=350, bbox_inches='tight')

'''
-------------------------------------------------------------------------------
CLIMTOLOGY FIGURE
-------------------------------------------------------------------------------
'''
clim = bx.subset_date('2016-01', '2021-01')
clim_O2 = np.nan*np.ones((N-1, M-1))
for i in range(N-1):
    for j in range(M-1):
        y1, y2 = ygrid[i], ygrid[i+1]
        x1, x2 = xgrid[j], xgrid[j+1]
        area = (clim.longitude > x1) & (clim.longitude < x2) &\
            (clim.latitude > y1) & (clim.latitude < y2)
        clim_O2[i,j] = clim[area].surf_doxy.mean()

fig = plt.figure(constrained_layout=True)
axes = [
    fig.add_subplot(311, projection=projection), 
    fig.add_subplot(312, projection=projection),
    fig.add_subplot(313, projection=projection)
]

for ax, data, cm, vmin, vmax, lab, title in zip(axes, [clim_O2, O2[-1,:,:], O2[-1,:,:]-clim_O2], [cmo.dense, cmo.dense, plt.cm.bwr], [250, 250, -60], [350, 350, 60], ['O$_2$ ($\mathregular{\mu}$mol kg$^{-1}$)', 'O$_2$ ($\mathregular{\mu}$mol kg$^{-1}$)', '$\Delta$O$_2$ ($\mathregular{\mu}$mol kg$^{-1}$)'], ['2016-2020', '2021', 'Difference: (2021) - (2016-2020)']):
    # plot map
    ax.set_extent(extent)
    ax.add_feature(cfeature.GSHHSFeature('low', 
        edgecolor='black', facecolor=cfeature.COLORS['land']))
    ax.patch.set_facecolor('lightgrey')
    ax.set_title(title, loc='left', fontweight='bold')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, color='white')
    gl.xlocator = mticker.FixedLocator([-65, -45])
    gl.ylocator = mticker.FixedLocator([55, 60])
    gl.right_labels = False
    gl.left_labels = True
    gl.top_labels = False
    gl.xformatter = lon_formatter
    gl.yformatter = lat_formatter

    m = ax.pcolormesh(X, Y, data, cmap=cm, transform=transform, vmin=vmin, vmax=vmax)
    cb = plt.colorbar(m, ax=ax)
    cb.set_label(lab)
fig.set_size_inches(1.5*fig.get_figwidth(), 1.5*fig.get_figheight())
fig.savefig(Path(f'../../figures/2022/gridded_surface_{boxsize}deg_delta_60bwr.png'), dpi=350, bbox_inches='tight')

'''
-------------------------------------------------------------------------------
SCATTER PLOTS
-------------------------------------------------------------------------------
'''