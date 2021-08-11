#!/usr/bin/python

from pathlib import Path
from netCDF4 import Dataset

import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks', palette='colorblind')
import cmocean.cm as cmo
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import bgcArgoDMQC as bgc
import argopy
argopy.set_options(mode='expert')

na_box   = [-62, -36, 52, 65]
extent   = [-70, -32, 40, 70]

bgc_index = bgc.get_index()
bgc_index = bgc_index[~bgc_index.parameters.isna()]
bgc_index = bgc_index[['DOXY' in p for p in bgc_index.parameters]]

na_index  = bgc_index[bgc_index.longitude > na_box[0]]
na_index  = na_index[na_index.longitude < na_box[1]]
na_index  = na_index[na_index.latitude > na_box[2]]
na_index  = na_index[na_index.latitude < na_box[3]]

na_index  = na_index[na_index.date > 2.02e13]

projection = ccrs.LambertConformal()

bath = Dataset(Path('/Users/gordonc/Documents/data/GEBCO/GEBCO_2020.nc'))

lons = bath['lon'][:]
lats = bath['lat'][:]
elev = bath['elevation'][:]

ix = np.logical_and(lons > na_box[0]-30, lons < na_box[1]+30)
iy = np.logical_and(lats > na_box[2]-30, lats < na_box[3]+30)

can_lons = lons[ix]
can_lats = lats[iy]
can_elev = elev[iy,:]
can_elev = can_elev[:,ix]

dx = np.arange(0,can_lons.shape[0],40)
dy = np.arange(0,can_lats.shape[0],40)

can_lons = can_lons[dx]
can_lats = can_lats[dy]
can_elev = can_elev[dy,:]
can_elev = can_elev[:,dx]
can_elev = -np.ma.masked_array(can_elev.data, can_elev > 0)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=projection)

ax.plot(na_index.longitude, na_index.latitude, 'o', markerfacecolor='yellow', markersize=2, markeredgecolor='k', markeredgewidth=0.2, transform=ccrs.PlateCarree())
ax.set_extent(extent)
ax.coastlines(resolution='50m')
ax.gridlines(draw_labels=True)
im = ax.contourf(can_lons, can_lats, can_elev,
                transform=ccrs.PlateCarree(),
                cmap=cmo.deep,
                vmin=0, extend='max')
# ax.plot([na_box[0], na_box[1], na_box[1], na_box[0], na_box[0]], [na_box[2], na_box[2], na_box[3], na_box[3], na_box[2]], linewidth=1, color='red', transform=ccrs.PlateCarree())
plt.colorbar(im, ax=ax, orientation='vertical', label='Depth (m)')
ax.add_feature(cfeature.LAND.with_scale('50m'))

fig.savefig(Path('../figures/na_map.png'), bbox_inches='tight', dpi=350)

fetcher = argopy.DataFetcher()
ds = fetcher.float([int(w) for w in na_index.wmo.unique()]).to_xarray().to_dataframe()

