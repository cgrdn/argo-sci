
import argopy
argopy.set_options(mode='expert')

import pandas as pd
import shapely
from netCDF4 import Dataset

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
sns.set_theme(style='ticks', context='paper', palette='colorblind')

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# load polygon to select data withing
poly = pd.read_csv('polygon_3300m.csv')
# shapely polygon
polygon = shapely.geometry.Polygon(poly)
# subset region to make polygon searching faster
polygon_box = [poly.longitude.min(), poly.longitude.max(), poly.latitude.min(), poly.latitude.max()]
# load argo index
ix = argopy.IndexFetcher().region(polygon_box).load().to_dataframe()
# points within polygon
ix = ix.loc[[polygon.contains(shapely.geometry.Point(x,y)) for x,y in zip(ix.longitude, ix.latitude)]]
ix['year'] = [d.year for d in ix.date]

# create figure
fig = plt.figure()

# set up map details
projection = ccrs.AzimuthalEquidistant(central_latitude=poly.latitude.mean(), central_longitude=poly.longitude.mean())
transform = ccrs.PlateCarree()
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
extent  = [-66, -42, 50, 64]
# create geo axis
ax = fig.add_subplot(projection=projection)
# format axis, ticks, labels
ax.set_extent(extent)
ax.add_feature(cfeature.GSHHSFeature('intermediate', 
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

# plot profiles
sns.scatterplot(data=ix, x='longitude', y='latitude', hue='year', ax=ax, transform=transform)
plt.show()