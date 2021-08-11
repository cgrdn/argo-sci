#!/usr/bin/python

import bgcArgoDMQC as bgc
import argopy
argopy.set_options(mode='expert')

# define "Northwest Atlantic" area
na_box   = [-62, -36, 52, 65]

# get index, subset to only floats with DOXY sensors
bgc_index = bgc.get_index()
bgc_index = bgc_index[~bgc_index.parameters.isna()]
bgc_index = bgc_index[['DOXY' in p for p in bgc_index.parameters]]

# subset index to only floats within our NW Atlantic box
na_index  = bgc_index[bgc_index.longitude > na_box[0]]
na_index  = na_index[na_index.longitude < na_box[1]]
na_index  = na_index[na_index.latitude > na_box[2]]
na_index  = na_index[na_index.latitude < na_box[3]]

# floats from after 2020 (as a starting point to minimize data load)
na_index  = na_index[na_index.date > 2.02e13]

# argopy data fetcher
fetcher = argopy.DataFetcher()
ds = fetcher.float([int(w) for w in na_index.wmo.unique()]).to_xarray().to_dataframe()

for wmo in ds.PLATFORM_NUMBER.unique():
    bgc.io.get_argo(wmo, local_path=bgc.ARGO_PATH)

