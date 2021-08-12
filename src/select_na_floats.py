#!/usr/bin/python

import bgcArgoDMQC as bgc
import argopy
argopy.set_options(mode='expert')
fetcher = argopy.IndexFetcher()

# define "Northwest Atlantic" area
na_box   = [-62, -36, 52, 65]

na_index = fetcher.region(na_box).to_dataframe()
bgc_index = bgc.get_index()
na_bgc_index = bgc_index[bgc_index.wmo.isin(na_index.wmo)]
na_bgc_index = na_bgc_index[~na_bgc_index.parameters.isna()]
na_doxy_index = na_bgc_index[['DOXY' in p for p in na_bgc_index.parameters]]
# floats from after 2020 (as a starting point to minimize data load)
recent_index = na_doxy_index[na_doxy_index.date > 2.02e13]

# argopy data fetcher
fetcher = argopy.DataFetcher()
# ds = fetcher

# for wmo in ds.PLATFORM_NUMBER.unique():
    # bgc.io.get_argo(wmo, local_path=bgc.ARGO_PATH)
