
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import argopandas as argo

wmos = [4902596, 4902597]

var_names = ['PRES', 'TEMP', 'PSAL', 'DOXY', 'CHLA', 'BBP700']

for wmo in wmos:
    ix   = argo.float(wmo).synthetic_prof
    up   = ix.subset_direction('asc')
    down = ix.subset_direction('desc')
    
    up['CYCLE']   = [int(f[-6:-3]) for f in up.file]
    down['CYCLE'] = [int(f[-7:-4]) for f in down.file]

    cycles = set(up.CYCLE.unique()).intersection(down.CYCLE.unique())

    for cycle in cycles:
        fig, ax  = plt.subplots()
        up_sub   = up.loc[up.CYCLE == cycle]
        down_sub = down.loc[down.CYCLE == cycle]

        up_data   = up_sub.levels
        down_data = down_sub.levels

        sns.lineplot(data=up_data, x='DOXY', y='PRES', sort=False, estimator=None, ax=ax)
        sns.lineplot(data=down_data, x='DOXY', y='PRES', sort=False, estimator=None, ax=ax)

        ax.set_ylim((200,0))