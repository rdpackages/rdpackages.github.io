# RDDENSITY

The `rddensity` package provides Stata and R implementations of manipulation tests employing local polynomial density estimation methods. This method is useful for falsification of Regression Discontinuity Designs, as well as for testing for self-selection or sorting in other contexts. This implementation provides hypothesis tests and bandwidth selectors for manipulation testing. 

This work was supported in part by the National Science Foundation through grants [SES-1357561](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1357561), [SES-1459931](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459931), [SES-1459967](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459967), [SES-1947662](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947662), [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805), and [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432).

## Queries and Requests

Please email: [rdpackages@googlegroups.com](mailto:rdpackages@googlegroups.com)

## Major Upgrades

This package was first released in Spring 2017, and had one major upgrade in Summer 2020.

- _Summer 2020 new features include_: (i) speed improvements; (ii) improved integration with [`lpdensity`](https://nppackages.github.io/lpdensity/); (iii) mass points in running variable adjustments; (iv) bandwidth selection adjustments for too few mass points in and/or overshooting of the support of the running variable; (v) density discontinuity plots with histogram and/or confidence bands; and (vi) binomial testing near cutoff as complementary discontinuity testing following results in [`rdlocrand`](https://rdpackages.github.io/rdlocrand/) methods (see references there for details).

## Python Implementation

To install/update in Python type:
```
pip install rddensity
```
- Help: [PyPI](https://pypi.org/project/rddensity/), [Documentation](https://github.com/rdpackages/rddensity/blob/master/Python/rddensity/docs/build/latex/rddensity.pdf)
- Replication: [Python script](https://github.com/rdpackages/rddensity/blob/master/Python/rddensity_illustration.py)

## R Implementation

To install/update in R type:
```
install.packages('rddensity')
```

Plotting employs [`lpdensity`](https://nppackages.github.io/lpdensity/), to install/update in R type:
```
install.packages('lpdensity')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/rddensity/rddensity.pdf), [CRAN repository](https://cran.r-project.org/package=rddensity).

- Replication: [R-script](https://github.com/rdpackages/rddensity/blob/master/R/rddensity_illustration.R), [density plot illustration](https://github.com/rdpackages/rddensity/blob/master/R/rddensity_plot_illustration.R), [senate data](https://github.com/rdpackages/rddensity/blob/master/R/rddensity_senate.csv).

## Stata Implementation

To install/update in Stata type:
```
net install rddensity, from(https://github.com/rdpackages/rddensity/blob/master/stata) replace
```

Plotting employs [`lpdensity`](https://nppackages.github.io/lpdensity/), to install/update in Stata type:
```
net install lpdensity, from(https://raw.githubusercontent.com/nppackages/lpdensity/master/stata) replace
```

- Help: [rddensity](https://github.com/rdpackages/rddensity/blob/master/stata/rddensity.pdf), [rdbwdensity](https://github.com/rdpackages/rddensity/blob/master/stata/rdbwdensity.pdf).

- Replication: [do-file](https://github.com/rdpackages/rddensity/blob/master/stata/rddensity_illustration.do), [do-file plot](https://github.com/rdpackages/rddensity/blob/master/stata/rddensity_plot_illustration.do), [data-senate](https://github.com/rdpackages/rddensity/blob/master/stata/rddensity_senate.dta).

## Repository

For source code and related files, visit [`rddensity` repository](https://github.com/rdpackages/rddensity/).


## References

For overviews and introductions, see [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Cattaneo, Jansson and Ma (2018): [Manipulation Testing based on Density Discontinuity](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2018_Stata.pdf).<br>
_Stata Journal_ 18(1): 234-261.

- Cattaneo, Jansson and Ma (2022): [lpdensity: Local Polynomial Density Estimation and Inference](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2022_JSS.pdf).<br>
_Journal of Statistical Software_, forthcoming.

### Technical and Methodological

- Cattaneo, Jansson and Ma (2020): [Simple Local Polynomial Density Estimators](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf).<br>
_Journal of the American Statistical Association_ 115(531): 1449-1455.<br>
[Supplemental appendix](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA--Supplement.pdf).

- Cattaneo, Jansson and Ma (2022): [Local Regression Distribution Estimators](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2022_JoE.pdf).<br>
_Journal of Econometrics_, forthcoming.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2022_JoE--Supplement.pdf).

<br><br>
