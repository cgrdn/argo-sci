
import argopandas as argo

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# grab index for this float
ix = argo.float(6902729).prof

# set up map, make figure
projection = ccrs.LambertConformal(central_latitude=80, central_longitude=15)
transform = ccrs.PlateCarree()
extent  = [-30, 60, 70, 85]
fig = plt.figure()
ax = fig.add_subplot(projection=projection)

# feature on the map
ax.set_extent(extent)
ax.add_feature(cfeature.GSHHSFeature('low', 
    edgecolor='black', facecolor=cfeature.COLORS['land']))
gl = ax.gridlines(crs=ccrs.PlateCarree(), 
    draw_labels=True, color='grey')

# plot the float trajectory
ax.plot(ix.longitude, ix.latitude, '-o', markersize=4, transform=transform)
fig.savefig('../figures/6902729_map.png', bbox_inches='tight', dpi=250)