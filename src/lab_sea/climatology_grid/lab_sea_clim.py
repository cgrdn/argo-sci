#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd

from scipy.interpolate import griddata

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
sns.set(style='ticks', palette='colorblind')
import cmocean.cm as cmo

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import argopandas as argo

# select lab sea region, doxy floats, last 5 years, ascending
lab_sea = [-67, -43, 55, 62.5]
bx = argo.synthetic_prof.subset_parameter('DOXY')\
    .subset_rect(lab_sea[2], lab_sea[0], lab_sea[3], lab_sea[1])\
    .subset_date('2016-01')\
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

prof = bx.prof[['JULD', 'LATITUDE', 'LONGITUDE']]
df = df.join(prof, on=['file', 'N_PROF'])
df['DATE'] = [pd.Timestamp(j -  365.25*20, unit='D') for j in df.JULD]
df['YEAR'] = [t.year for t in df.DATE]
allowed_flags = [b'1', b'2', b'3', b'5']

'''
-------------------------------------------------------------------------------
FIGURE 1 - AVERAGE OXYGEN BY YEAR
-------------------------------------------------------------------------------

'''

surface = df[df.PRES < 100]
fig, axes = plt.subplots(2, 1, sharex=True)
sns.boxplot(x='YEAR', y='DOXY', data=surface[surface.DOXY_QC.isin(allowed_flags)], fliersize=0.5, ax=axes[0])
sns.boxplot(x='YEAR', y='TEMP', data=surface[surface.TEMP_QC.isin(allowed_flags)], fliersize=0.5, ax=axes[1])

axes[0].set_xlabel('')
axes[1].set_xlabel('')
axes[0].set_ylim(top=420)
axes[0].set_ylabel('Diss. Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)')
axes[1].set_ylabel(f'Temperature ({chr(176)}C)')
axes[0].set_title('Surface Values ($P < 100$dbar)', loc='left')

fig.set_size_inches(fig.get_figwidth()/1.5, fig.get_figheight())
fig.savefig(Path('../../figures/lab_sea/2022/mean_surface_oxygen_temperature.png'), bbox_inches='tight', dpi=350)
plt.close(fig)

'''
-------------------------------------------------------------------------------
FIGURE 2 - PCOLOR BY YEAR
-------------------------------------------------------------------------------
'''

# choose projection/transform
projection = ccrs.LambertConformal(central_latitude=55, central_longitude=-55)
transform = ccrs.PlateCarree()

nyear = df.YEAR.unique().shape[0]

fig = plt.figure()
axes = np.array([
    [fig.add_subplot(nyear, 3, 3*i+1) for i in range(nyear)],
    [fig.add_subplot(nyear, 3, 3*i+2) for i in range(nyear)],
    [fig.add_subplot(nyear, 3, 3*i+3, projection=projection) for i in range(nyear)]

]).T
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()

extent  = [-67, -42, 54.5, 63]

df['DAY'] = [t.dayofyear for t in df.DATE]
upper = df[df.PRES < 250]
maps = {}
for axrow, year in zip(axes, np.sort(df.YEAR.unique())):
    yr = upper[upper.YEAR == year]
    for ax, v, cmap, vmin, vmax in zip(axrow, ['TEMP', 'DOXY'], [cmo.thermal, cmo.dense_r], [3, 280], [8, 340]):
        qc = yr[yr[v + '_QC'].isin(allowed_flags)]
        maps[f'{v}-{year}'] = ax.scatter(x=qc.DAY, y=qc.PRES, c=qc[v], s=3, cmap=cmap, vmin=vmin, vmax=vmax)
    
    coords = bx.subset_date(f'{year}-01', f'{year+1}-01')
    axrow[-1].plot(coords.longitude, coords.latitude, 'o', color='k', markersize=1, transform=transform)
    axrow[-1].set_extent(extent)
    axrow[-1].add_feature(cfeature.GSHHSFeature('low', edgecolor='black', facecolor='lightgray'))
    gl = axrow[-1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
    gl.xlocator = mticker.FixedLocator([-65, -45])
    gl.ylocator = mticker.FixedLocator([55, 60])
    gl.left_labels = False
    gl.top_labels = False
    gl.xformatter = lon_formatter
    gl.yformatter = lat_formatter

    axrow[0].set_title(f'{year}', loc='left', fontsize=8, fontweight='bold')
    axrow[-1].set_title('$N_{prof}$ = %d' % coords.shape[0], loc='right', fontsize=8, fontweight='bold')

for ax in axes[:,1]:
    ax.set_yticklabels([])

for ax in axes[:nyear-1, :].flatten():
    ax.set_xticklabels([])

for ax in axes[nyear-1, :2]:
    ax.set_xlabel('Julian Day')

for ax in axes[:, 2]:
    ax.yaxis.tick_right()

for ax in axes[:, :2].flatten():
    ax.set_ylim((200, 0))
    ax.set_xlim((0, 366))

tcax = fig.add_axes([-0.07, 0.545, 0.025, 0.4])
dcax = fig.add_axes([-0.07, 0.075, 0.025, 0.4])
plt.colorbar(maps['TEMP-2021'], cax=tcax, label=f'Temperature ({chr(176)}C)')
plt.colorbar(maps['DOXY-2021'], cax=dcax, label='Diss. Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)')
tcax.yaxis.set_ticks_position('left')
dcax.yaxis.set_ticks_position('left')
tcax.yaxis.set_label_position('left')
dcax.yaxis.set_label_position('left')

fig.set_size_inches(fig.get_figwidth(), fig.get_figheight()*1.75)
fig.tight_layout()
fig.savefig(Path('../../figures/lab_sea/2022/pcolor_and_maps.png'), bbox_inches='tight', dpi=200)

'''
-------------------------------------------------------------------------------
FIGURE 3 - AVERAGE OXYGEN OVER LAST 5 YEARS
-------------------------------------------------------------------------------

'''

clim = df[df.YEAR != 2022]
vec_pres, vec_time = np.arange(0, 2000, 1), np.arange(0, 365, 1)
grid_pres, grid_time = np.meshgrid(vec_pres, vec_time)
dix = (clim.DOXY_QC.isin(allowed_flags)) & (clim.DOXY.notna()) & (clim.PRES.notna()) & (clim.DAY.notna())
tix = (clim.TEMP_QC.isin(allowed_flags)) & (clim.TEMP.notna()) & (clim.PRES.notna()) & (clim.DAY.notna())
clim_doxy = griddata((clim[dix].DAY.values, clim[dix].PRES.values), clim[dix].DOXY.values, (grid_time, grid_pres), method='linear')
clim_temp = griddata((clim[tix].DAY.values, clim[tix].PRES.values), clim[tix].TEMP.values, (grid_time, grid_pres), method='linear')

this_year = df[df.YEAR == 2022]
dix = (this_year.DOXY_QC.isin(allowed_flags)) & (this_year.DOXY.notna()) & (this_year.PRES.notna()) & (this_year.DAY.notna())
tix = (this_year.TEMP_QC.isin(allowed_flags)) & (this_year.TEMP.notna()) & (this_year.PRES.notna()) & (this_year.DAY.notna())
this_year_doxy = griddata((this_year[dix].DAY.values, this_year[dix].PRES.values), this_year[dix].DOXY.values, (grid_time, grid_pres), method='linear')
this_year_temp = griddata((this_year[tix].DAY.values, this_year[tix].PRES.values), this_year[tix].TEMP.values, (grid_time, grid_pres), method='linear')

fig, axes = plt.subplots(2, 3, sharex=True, sharey=True)
maps = 6*['']
maps[0] = axes[0, 0].pcolor(grid_time, grid_pres, clim_temp, cmap=cmo.thermal, vmin=3, vmax=8, shading='auto')
maps[1] = axes[0, 1].pcolor(grid_time, grid_pres, this_year_temp, cmap=cmo.thermal, vmin=3, vmax=8, shading='auto')
maps[2] = axes[0, 2].pcolor(grid_time, grid_pres, this_year_temp - clim_temp, cmap=cmo.balance, vmin=-1, vmax=1, shading='auto')
maps[3] = axes[1, 0].pcolor(grid_time, grid_pres, clim_doxy, cmap=cmo.dense_r, vmin=240, vmax=340, shading='auto')
maps[4] = axes[1, 1].pcolor(grid_time, grid_pres, this_year_doxy, cmap=cmo.dense_r, vmin=240, vmax=340, shading='auto')
maps[5] = axes[1, 2].pcolor(grid_time, grid_pres, this_year_doxy - clim_doxy, cmap=cmo.balance, vmin=-30, vmax=30, shading='auto')

labels = [
    f'T ({chr(176)}C)',
    f'T ({chr(176)}C)',
    f'$\Delta$T ({chr(176)}C)',
    'O$_2$ ($\mathregular{\mu}$mol kg$^{-1}$)', 
    'O$_2$ ($\mathregular{\mu}$mol kg$^{-1}$)', 
    '$\Delta$O$_2$ ($\mathregular{\mu}$mol kg$^{-1}$)']

for m, ax, l in zip(maps, axes.flatten(), labels):
    plt.colorbar(m, ax=ax, label=l)

for ax in axes.flatten():
    ax.set_ylim((200,0))

axes[0,0].set_ylabel('Pessure (dbar)')
axes[1,0].set_ylabel('Pessure (dbar)')
axes[1,0].set_xlabel('Julian Day')
axes[1,1].set_xlabel('Julian Day')
axes[1,2].set_xlabel('Julian Day')

fig.set_size_inches(fig.get_figwidth()*2, fig.get_figheight())
fig.savefig(Path('../../figures/lab_sea/2022/pcolor_clim_and_delta_upper.png'), bbox_inches='tight', dpi=200)

for ax in axes.flatten():
    ax.set_ylim((2000,0))
fig.savefig(Path('../../figures/lab_sea/2022/pcolor_clim_and_delta_full.png'), bbox_inches='tight', dpi=200)

'''
-------------------------------------------------------------------------------
FIGURE 4 - PROFILES
-------------------------------------------------------------------------------

'''
fig, axes = plt.subplots(2, nyear+1, sharey=True)

vec_pres = np.arange(0, 200, 5)

for axrow, year in zip(axes.T, np.sort(df.YEAR.unique())):
    yr = upper[upper.YEAR == year]
    mean_vals = {'TEMP':np.nan*np.ones((vec_pres.shape[0] - 1,)), 'DOXY':np.nan*np.ones((vec_pres.shape[0] - 1,))}
    for j, ax, v, ll, ul, l in zip([0, 1], axrow, ['TEMP', 'DOXY'], [0, 220], [10, 360], [f'T ({chr(176)}C)', 'O$_2$ ($\mathregular{\mu}$mol kg$^{-1}$)']):
        qc = yr[yr[v + '_QC'].isin(allowed_flags)]
        for i in range(vec_pres.shape[0]-1):
            ix = (qc.PRES > vec_pres[i]) & (qc.PRES < vec_pres[i+1])
            mean_vals[v][i] = qc[ix][v].mean()
        for d in qc.DATE.unique():
            prof = qc[qc.DATE == d]
            ax.plot(prof[v], prof.PRES, color='grey', alpha=0.2)
        ax.plot(mean_vals[v], vec_pres[:-1] + 2.5, linewidth=2, color='k')
        axes[j,-1].plot(mean_vals[v], vec_pres[:-1] + 2.5, linewidth=2, label=f'{year}')
        ax.set_xlim((ll, ul))

for ax in axes.flatten():
    ax.set_ylim((200, 0))

for ax in axes[1,:]:
    ax.set_xticks([250, 325])

axes[0,-1].legend(loc=3, fontsize=8, bbox_to_anchor=(1.05, 0.1))
axes[0,0].set_ylabel('Pressure (dbar)')
axes[0,0].set_xlabel(f'T ({chr(176)}C)')
axes[1,0].set_ylabel('Pressure (dbar)')
axes[1,0].set_xlabel('O$_2$ ($\mathregular{\mu}$mol kg$^{-1}$)')
fig.set_size_inches(fig.get_figwidth()*1.5, fig.get_figheight())
fig.tight_layout()
fig.savefig(Path('../../figures/lab_sea/2022/profiles_and_average.png'), bbox_inches='tight', dpi=200)
