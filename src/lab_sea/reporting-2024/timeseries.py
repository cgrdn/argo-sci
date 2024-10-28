#!/usr/bin/python

from pathlib import Path

import numpy as np
import pandas as pd
import argopy

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
sns.set_theme(style='ticks', palette='colorblind')

import shapely


# load polygon to select data withing
poly = pd.read_csv('polygon_3300m.csv')
# shapely polygon
polygon = shapely.geometry.Polygon(poly)
# subset region to make polygon searching faster
polygon_box = [poly.longitude.min(), poly.longitude.max(), poly.latitude.min(), poly.latitude.max()]
depths = [0, 200]
dates = ['2020-01-01', '2025-01-01']
# load argo index
index = argopy.IndexFetcher().region(polygon_box).load()
ix = index.to_dataframe()
# points within polygon
ix = ix.loc[[polygon.contains(shapely.geometry.Point(x,y)) for x,y in zip(ix.longitude, ix.latitude)]]
ix['year'] = [d.year for d in ix.date]


# get corresponding data to index above
data = argopy.DataFetcher(src='erddap').region(polygon_box + depths + dates).load()
ds = data.to_xarray()
ds.argo.teos10(['SA', 'SIG0'])
df = data.to_dataframe()

variables = ['TEMP', 'SA', 'SIG0']

fig, axes = plt.subplots(3, 1, sharex=True)
for v, ax in zip(variables, axes):
    sns.lineplot(data=df, x='date', y=v)

plt.show()