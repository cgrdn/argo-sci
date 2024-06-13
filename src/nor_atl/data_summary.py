
from pathlib import Path

import numpy as np
import pandas as pd
import argopandas as argo

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks', palette='colorblind')
import cmocean.cm as cmo

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from netCDF4 import Dataset

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
ix = argo.synthetic_prof.subset_rect(*reshape_extent(na_box)).subset_parameter('DOXY').subset_date('2022-10', '2023-01')
prof = ix.prof
ix = ix.set_index('file')
ix = ix.join(prof[['PLATFORM_NUMBER', 'CYCLE_NUMBER']])
ix['LAST_REPORTED'] = [ix.loc[ix.PLATFORM_NUMBER == p].date.iloc[-1] for p in ix.PLATFORM_NUMBER]
# df = ix.levels[varnames]
ix['deployment_date'] = [argo.float(f).synthetic_prof.date.min() for f in ix.PLATFORM_NUMBER]
ix['deployment_year'] = [d.year for d in ix.deployment_date]
last_cycle = ix.loc[[True if c == ix.loc[ix.PLATFORM_NUMBER == p].CYCLE_NUMBER.max() else False for p in ix.PLATFORM_NUMBER.unique() for c in ix.loc[ix.PLATFORM_NUMBER == p].CYCLE_NUMBER]]

# map figure
fig = plt.figure()

proj = ccrs.PlateCarree(central_longitude=-30)
trans = ccrs.PlateCarree()
ax = fig.add_subplot(projection=trans)

# bathymetry
bath_file = Path('/Users/GordonC/Documents/data/GEBCO/GEBCO_2020.nc')
bath = Dataset(bath_file)
blat = bath['lat'][:]
blon = bath['lon'][:]
elev = bath['elevation'][:]

px = np.logical_and(blon > na_box[0] - 10, blon < na_box[1] + 10)
py = np.logical_and(blat > na_box[2] - 10, blat < na_box[3] + 10)

blon = blon[px]
blat = blat[py]
elev = elev[py,:]
elev = elev[:,px]
elev = -np.ma.masked_array(elev.data, elev > 0)

n = 10
blon = blon[::n]
blat = blat[::n]
elev = elev[::n, :]
elev = elev[:, ::n]

# add bathymetry
im = ax.contourf(
    blon, blat, elev, list(range(0, 3800, 200)),
    transform=trans,
    cmap=cmo.deep,
    vmin=0, extend='max'
)
cb = plt.colorbar(im, ax=ax, orientation='horizontal', label='Depth (m)')

sns.lineplot(data=ix, x='longitude', y='latitude', 
    style='PLATFORM_NUMBER', dashes=False, 
    sort=False, estimator=None,
    ax=ax, transform=trans, legend=False, 
    alpha=0.25, color='grey', zorder=1)
sns.scatterplot(data=last_cycle, x='longitude', y='latitude', hue='deployment_year', palette='Greys', transform=trans, zorder=2, legend=False)
ax.set_extent(na_box)
ax.add_feature(cfeature.GSHHSFeature('high'))

# format map
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
gl = ax.gridlines(crs=trans, draw_labels=True)
gl.right_labels = False
gl.xformatter = lon_formatter
gl.yformatter = lat_formatter

fig.set_size_inches(fig.get_figwidth()*2, fig.get_figheight()*2)
fig.savefig(Path('../../figures/na/map_profiles_2018-2023.png', dpi=350, bbox_inches='tight'))
plt.close(fig)
# plt.show()