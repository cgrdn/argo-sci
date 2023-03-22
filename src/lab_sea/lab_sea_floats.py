#!/usr/bin/python

from pathlib import Path

import numpy as np
from netCDF4 import Dataset

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cmocean.cm as cmo
import seaborn as sns
sns.set(style='ticks', context='paper', palette='colorblind')

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import argopandas as argo

# define rectangle over lab sea
lab_sea = [-67, -43, 48.5, 66.25]
extent  = [-68, -40, 40, 70]
# choose projection/transform
projection = ccrs.LambertConformal(central_latitude=55, central_longitude=-55)
transform = ccrs.PlateCarree()
# get bathymetry
bath_file = Path('/Users/GordonC/Documents/data/GEBCO/GEBCO_2020.nc')
bath = Dataset(bath_file)
blat = bath['lat'][:]
blon = bath['lon'][:]
elev = bath['elevation'][:]

ix = np.logical_and(blon > extent[0] - 10, blon < extent[1] + 10)
iy = np.logical_and(blat > extent[2] - 10, blat < extent[3] + 10)

blon = blon[ix]
blat = blat[iy]
elev = elev[iy,:]
elev = elev[:,ix]
elev = -np.ma.masked_array(elev.data, elev > 0)

# create figure, add coastlines
fig = plt.figure()
ax = fig.add_subplot(projection=projection)
ax.set_extent(extent)
ax.add_feature(cfeature.GSHHSFeature('intermediate', edgecolor='black', facecolor='lightgray'))

# add bathymetry
im = ax.contourf(
    blon, blat, elev, list(range(0, 5000, 200)),
    transform=transform,
    cmap=cmo.deep,
    vmin=0, extend='max'
)
plt.colorbar(im, ax=ax, orientation='vertical', label='Depth (m)')

# get floats
bx = argo.bio_prof.subset_rect(lab_sea[2], lab_sea[0], lab_sea[3], lab_sea[1])
bx = bx.subset_date('2020', '2021')

# plot the track of each float
bx['wmo'] = [f.split('/')[1] for f in bx.file]
for wmo in bx.wmo.unique():
    sub = bx[bx.wmo == wmo]
    ax.plot(sub.longitude, sub.latitude, '-', color='lightgray', transform=transform, label=None, alpha=0.5)
    ax.plot(sub.longitude, sub.latitude, '.', color='black', transform=transform, label=None, markersize=1)

# plot the latest profile with certain markers to mark each variable
vars_of_interest = ['DOXY', 'RAW_DOWNWELLING_PAR', 'CHLA', 'PH_IN_SITU_TOTAL']
marker_types = [mpl.markers.CARETLEFTBASE,mpl.markers.CARETUPBASE,mpl.markers.CARETRIGHTBASE,mpl.markers.CARETDOWNBASE]
colors = ['#2E86C1', '#F1C40F', '#28B463', '#884EA0']
labels = ['Dissolved Oxygen', 'PAR', 'Bio-optics', 'pH']

for v, m, c, l in zip(vars_of_interest, marker_types, colors, labels):
    sub = bx.subset_parameter(v)
    ax.plot(np.nan, np.nan, marker=m, markersize=7.5, markeredgecolor='white', color=c, transform=transform, label=l)
    for wmo in sub.wmo.unique():
        flt = sub[sub.wmo == wmo]
        last_profile = flt[flt.date == flt.date.max()]
        ax.plot(last_profile.longitude, last_profile.latitude, marker=m, markersize=7.5, markeredgecolor='white', color=c, transform=transform, label=None)

ax.legend(fontsize=6, loc=4)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
gl.xlocator = mticker.FixedLocator([-70, -55, -40])
gl.right_labels = False
gl.xformatter = lon_formatter
gl.yformatter = lat_formatter

fig.savefig(Path('../../figures/lab_sea_map_w_variables.png'), dpi=450, bbox_inches='tight')
plt.close(fig)
