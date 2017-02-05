# EUregmort

[![Build Status](https://travis-ci.org/klpn/EUregmort.jl.svg?branch=master)](https://travis-ci.org/klpn/EUregmort.jl)

[![Coverage Status](https://coveralls.io/repos/klpn/EUregmort.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/klpn/EUregmort.jl?branch=master)

[![codecov.io](http://codecov.io/github/klpn/EUregmort.jl/coverage.svg?branch=master)](http://codecov.io/github/klpn/EUregmort.jl?branch=master)

This package can be used to analyze cause-specific mortality in NUTS regions in
the European Union, and, more specifically, to explore the relative importance
of different causes of death. It uses the table with crude age-specific death
rates in NUTS 2 regions averaged over 3 years, `hlth_cd_ycdr2`, which is
available via the [Eurostat database](http://ec.europa.eu/eurostat/data/database).
The table, converted to CSV, is included in the packages at `data/hlth_cd_ycdr2.csv`.
For data about the different NUTS regions, the file
`data/NUTS_AT_2013.csv` (included in the [shapefile
archive](http://ec.europa.eu/eurostat/cache/GISCO/geodatafiles/NUTS_2013_03M_SH.zip)
for the default plots, as described below) is used.

To define an array of the NUTS 2 regions in the Nordic countries with data in
`hlth_cd_ycdr2` (Denmark, Finland, Norway and Sweden):

```julia
using EUregmort
nordnuts2 = vcat(map((x)->nuts2ids(x), ["DK"; "FI"; "NO"; "SE"])...)
```

Use this array to plot the correlation between female and male proportion of
mortality from circulatory causes in the Nordic countries (see the
`data/CL_ICD10_20170129_155451.csv`, based on Eurostat metadata, for
information about causes of death):

```julia
caprop_regsexplot(nordnuts2, "TOTAL", "I")
```

Using [cartopy](https://github.com/SciTools/cartopy), it is also possible to plot
maps showing regions with a lower or higher proportion of deaths from a given
cause. To do this, it is necessary to have shapefiles with the different
regions. You can download files with 1:3 million scale as a [ZIP
archive](http://ec.europa.eu/eurostat/cache/GISCO/geodatafiles/NUTS_2013_03M_SH.zip),
which can be unzipped in the `data` directory. However, the shapefiles from
Eurostat use the EPSG:4258 projection, which is not suited for use with
Cartopy. To convert them to EPSG:3034 projection (default settings for the
package), use e.g. [GDAL](http://www.gdal.org), and run, in the
`data/NUTS_2013_03M_SH/data` directory:

```shell
ogr2ogr -f "ESRI Shapefile" -t_srs EPSG:3034 -s_srs EPSG:4258 NUTS_RG_03M_2013_3034.shp
 NUTS_RG_03M_2013.shp
```

To plot a map of female proportion of deaths due to circulatory causes in the
Nordic countries:

```julia
caprop_mapplot(nordnuts2, "F", "TOTAL", "I")
```

By default, death rates from all causes are used as denominator. However, you
can give another cause as the fourth argument. For example, to plot a map of
male deaths due to circulatory causes relative to neoplasms in the Nordic
countries:

```julia
caprop_mapplot(nordnuts2, "M", "TOTAL", "I", ca2 = "C00-D48")
```

You can also plot the death rates themselves by specifying `"pop"` as denominator:

```julia
caprop_mapplot(nordnuts2, "M", "TOTAL", "I", ca2 = "pop")
```

When plotting death rates, it is interesting to compare average rates over
age groups. Use the function `meanrate` to define an alternative dataframe. Age groups
with average rates are prefixed `YM`. To plot a map of regions with higher
average female mortality from neoplasms over the ages from 45 to 64.

```julia
caprop_mapplot(nordnuts2, "F", "YM45-64", "C00-D48", ca2 = "pop",
inframe = meanrate("Y45-49", "Y60-64"))
```
