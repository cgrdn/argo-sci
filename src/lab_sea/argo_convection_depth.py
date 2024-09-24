
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(context='paper', style='ticks', palette='colorblind')

# import argopy
import argopandas as argo
import gsw

lab_sea = [-67, -43, 55, 62.5]

# argopy.set_options(mode='expert')

# ix = argopy.ArgoIndex().load()
# ix.search_lat_lon_tim(lab_sea + ['2024-06-01', '2024-07-01'])

ix = argo.prof.subset_rect(lab_sea[2], lab_sea[0], lab_sea[3], lab_sea[1])\
    .subset_date('2024-05-23')\
    .subset_direction('asc')

data = ix.levels[['PRES', 'PRES_QC', 'TEMP', 'TEMP_QC', 'PSAL', 'PSAL_QC']]
data = data.loc[:,0,:]
prof = ix.prof[['PLATFORM_NUMBER', 'JULD', 'LATITUDE', 'LONGITUDE', 'CYCLE_NUMBER']]
data = data.join(prof, on=('file', 'N_PROF'))

data['SA'] = gsw.SA_from_SP(data['PSAL'], data['PRES'], data['LONGITUDE'], data['LATITUDE'])
data['CT'] = gsw.CT_from_t(data['SA'], data['TEMP'], data['PRES'])
data['sigma2'] = gsw.sigma2(data['SA'], data['CT'])

lsw_layer = 36.85
g = sns.lineplot(data=data, x='sigma2', y='PRES', hue='file', sort=False, legend=False)
g.axes.set_ylim((2050, -50))
# g.axes.set_ylim((37, 35.6))
g.axes.axvline(lsw_layer, color='k')
for fn in ix.file:
    sub = data.loc[fn]
    ix = sub.sigma2 > lsw_layer
    if ix.any():
        lsw_pres = sub.loc[ix].PRES.iloc[0]
        print(lsw_pres)
        g.axes.axhline(lsw_pres, color='k', alpha=0.5)
plt.show()