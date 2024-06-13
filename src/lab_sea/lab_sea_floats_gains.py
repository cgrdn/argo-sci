
import numpy as np
import pandas as pd

import argopandas as argo
import bgcArgoDMQC as bgc

lab_sea = [-67, -43, 55, 62.5]
bx = argo.synthetic_prof.subset_parameter('DOXY')\
    .subset_rect(lab_sea[2], lab_sea[0], lab_sea[3], lab_sea[1])\
    .subset_date('2016-01', '2024-01')\
    .subset_direction('asc')

# this float is bad
bx = bx[~bx.file.str.contains('4901779')]
bx = bx[~bx.file.str.contains('4902581')]

bx['wmo'] = [int(s.split('/')[1]) for s in bx.file]

gdf = pd.read_csv('lab_sea_gains.csv')

wmo_list = [wmo for wmo in gdf.wmo]
gain_list = [g for g in gdf.gain]

for wmo in bx.wmo.unique():

    if wmo not in wmo_list:

        bgc.io.get_argo(wmo, local_path=bgc.io.Path.ARGO_PATH)
        flt = bgc.sprof(wmo)
        gains = flt.calc_gains(ref='WOA')

        wmo_list.append(wmo)
        gain_list.append(np.nanmean(gains))

gdf = pd.DataFrame(dict(wmo=wmo_list, gain=gain_list))
gdf.to_csv('lab_sea_gains.csv')