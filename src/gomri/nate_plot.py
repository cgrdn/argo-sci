
from pathlib import Path
from netCDF4 import Dataset

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks', palette='colorblind')
import cmocean.cm as cmo

import gsw

data_files = list(Path('../../data/gomri/').glob('*.nc'))
# time for Hurricane Nate + a 1 day on each end
storm_range = [
    pd.Timestamp(year=2017, month=10, day=3), 
    pd.Timestamp(year=2017, month=10, day=9)
]

for fn in [data_files[-3]]:
    nc = Dataset(fn)
    # build dataframe
    df = pd.DataFrame(dict(
        cycles = np.tile(nc['CYCLE_NUMBER'][:], (nc.dimensions['N_LEVELS'].size, 1)).T.flatten(),
        juld = np.tile(nc['JULD'][:], (nc.dimensions['N_LEVELS'].size, 1)).T.flatten(),
        time = [pd.Timestamp(t, unit='D') - pd.Timedelta(days=365.25*20) if t != 999999. else np.nan for t in nc['MTIME'][:].flatten().data],
        date_number = nc['MTIME'][:].flatten(),
        latitude = nc['LATITUDE_GRID'][:].flatten(),
        longitude = nc['LONGITUDE_GRID'][:].flatten(),
        pres = nc['PRES'][:].flatten(),
        psal = nc['PSAL'][:].flatten(),
        temp = nc['TEMP'][:].flatten(),
        doxy = nc['DOXY_ADJUSTED'][:].flatten(),
        chla = nc['CHLA'][:].flatten(),
        bbp  = nc['BBP700'][:].flatten()*1000,
    ))

    df.psal.loc[df.psal < 34] = np.nan
    df = df.loc[(df.time > storm_range[0]) & (df.time < storm_range[1])]
    df = df.loc[df.pres < 120]
    df['pden'] = gsw.pot_rho_t_exact(
        gsw.SA_from_SP(df.psal.values, df.pres.values, df.longitude.values, df.latitude.values),
        df.temp.values, df.pres.values, p_ref=0
    ) - 1000

    vmin = df.juld.min()-0.1
    vmin = None if np.isnan(vmin) else vmin
    vmax = df.juld.max()+0.1
    vmax = None if np.isnan(vmax) else vmax

    fig, axes = plt.subplots(1, 3, sharey=True)

    for v, l, ax in zip(['temp', 'psal', 'pden'], [f'Temperature ({chr(176)}C)', 'Abs. Salinity (g kg$^{-1}$)', 'Pot. Density (kg m$^{-3}$)'], axes):
        g = ax.scatter(df[v], df['pres'], c=df['date_number'], cmap=cmo.thermal, s=10, vmin=vmin, vmax=vmax)
        ax.set_xlabel(l)

    axes[0].set_ylabel('Pressure (dbar)')
    ax.set_ylim((100,0))

    dates = df.juld.unique()
    date_ts = [pd.Timestamp(t, unit='D') - pd.Timedelta(days=365.25*20) for t in dates]
    date_labels = [t.strftime('%d %b, %Y %H:%M') for t in date_ts]

    N = 8
    rect = [0.9125, 0.1075, 0.0325, 0.775]
    cax = fig.add_axes(rect)
    cax2 = fig.add_axes(rect)
    cax2.set_axis_off()
    cb = plt.colorbar(g, ax=ax, cax=cax)
    cax2.plot(len(dates)*[1.05], dates, '_', color='k', zorder=99999)
    cax2.set_xlim((1.05, 1.15))
    cax2.set_ylim((vmin, vmax))
    cb.set_ticks(dates[::N])
    cb.ax.set_yticklabels(date_labels[::N])

    fig.set_size_inches(fig.get_figwidth(), fig.get_figheight()*0.75)
    fig.savefig(f'../../figures/gomri/T-S-rho/profiles_{fn.stem.replace("_Sprof", ".png")}', bbox_inches='tight', dpi=300)

    fig, axes = plt.subplots(1, 3, sharey=True)

    for v, l, ax in zip(['doxy', 'temp', 'pden'], ['Diss. Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)', f'Temperature ({chr(176)}C)', 'Pot. Density (kg m$^{-3}$)'], axes):
        g = ax.scatter(df[v], df['pres'], c=df['date_number'], cmap=cmo.thermal, s=10, vmin=vmin, vmax=vmax)
        ax.set_xlabel(l)

    axes[0].set_ylabel('Pressure (dbar)')
    ax.set_ylim((100,0))

    rect = [0.9125, 0.1075, 0.0325, 0.775]
    cax = fig.add_axes(rect)
    cax2 = fig.add_axes(rect)
    cax2.set_axis_off()
    cb = plt.colorbar(g, ax=ax, cax=cax)
    cax2.plot(len(dates)*[1.05], dates, '_', color='k', zorder=99999)
    cax2.set_xlim((1.05, 1.15))
    cax2.set_ylim((vmin, vmax))
    cb.set_ticks(dates[::N])
    cb.ax.set_yticklabels(date_labels[::N])

    fig.set_size_inches(fig.get_figwidth(), fig.get_figheight()*0.75)
    fig.savefig(f'../../figures/gomri/DO-T-rho/profiles_{fn.stem.replace("_Sprof", ".png")}', bbox_inches='tight', dpi=300)

    fig, axes = plt.subplots(2, 3, sharey=True)

    for v, l, ax in zip(['temp', 'psal', 'pden'], [f'Temperature ({chr(176)}C)', 'Abs. Salinity (g kg$^{-1}$)', 'Pot. Density (kg m$^{-3}$)'], axes[0,:]):
        g = ax.scatter(df[v], df['pres'], c=df['date_number'], cmap=cmo.thermal, s=10, vmin=vmin, vmax=vmax)
        ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
        ax.xaxis.set_label_position('top')
        ax.set_xlabel(l)

    axes[0,0].set_ylabel('Pressure (dbar)')
    ax.set_ylim((100,0))

    for v, l, ax in zip(['doxy', 'chla', 'bbp'], ['Diss. Oxygen ($\mathregular{\mu}$mol kg$^{-1}$)', 'Chl $a$ (mg m$^{-3}$)', 'Backscatter (x10$^3$ m$^{-1}$)'], axes[1,:]):
        g = ax.scatter(df[v], df['pres'], c=df['date_number'], cmap=cmo.thermal, s=10, vmin=vmin, vmax=vmax)
        ax.set_xlabel(l)

    axes[1,0].set_ylabel('Pressure (dbar)')
    ax.set_ylim((100,0))

    rect = [0.9125, 0.1075, 0.0325, 0.775]
    cax = fig.add_axes(rect)
    cax2 = fig.add_axes(rect)
    cax2.set_axis_off()
    cb = plt.colorbar(g, ax=ax, cax=cax)
    cax2.plot(len(dates)*[1.05], dates, '_', color='k', zorder=99999)
    cax2.set_xlim((1.05, 1.15))
    cax2.set_ylim((vmin, vmax))
    cb.set_ticks(dates[::N])
    cb.ax.set_yticklabels(date_labels[::N])

    fig.set_size_inches(fig.get_figwidth(), fig.get_figheight()*2*0.75)
    fig.savefig(f'../../figures/gomri/T-S-rho_DO-CHLA-BBP/profiles_{fn.stem.replace("_Sprof", ".png")}', bbox_inches='tight', dpi=300)

    plt.close('all')