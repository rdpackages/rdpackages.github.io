# RDDENSITY

The `rddensity` package provides Stata and R implementations of manipulation tests employing local polynomial density estimation methods. This method is useful for falsification of Regression Discontinuity Designs, as well as for testing for self-selection or sorting in other contexts. This implementation provides hypothesis tests and bandwidth selectors for manipulation testing. 

This work was supported in part by the National Science Foundation through grants [SES-1357561](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1357561), [SES-1459931](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459931), [SES-1459967](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459967), [SES-1947662](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947662), [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805), and [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432).

## Major Upgrades

This package was first released in Spring 2017, and had one major upgrade in Summer 2020.

- _Summer 2020 new features include_: (i) speed improvements; (ii) improved integration with lpdensity; (iii) mass points in running variable adjustments; (iv) bandwidth selection adjustments for too few mass points in and/or overshooting of the support of the running variable; (v) density discontinuity plots with histogram and/or confidence bands; and (vi) binomial testing near cutoff as complementary discontinuity testing following results in rdlocrand methods (see references there for details).

## Stata Implementation

To install/update in Stata type:
```
net install rddensity, from(https://raw.githubusercontent.com/rdpackages/rddensity/master/stata) replace
```

- Help: [rddensity](stata/rddensity.pdf), [rdbwdensity](stata/rdbwdensity.pdf).
- Replication: [do-file](stata/rddensity_illustration.do), [do-file plot](stata/rddensity_plot_illustration.do), [data-senate](stata/rddensity_senate.dta).

## R Implementation

To install/update in R type:
```
install.packages('rddensity')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/rddensity/rddensity.pdf), [CRAN repository](https://cran.r-project.org/package=rddensity).

- Replication: [R-script](R/rddensity_illustration.R), [rdplot illustration](R/rddensity_plot_illustration.R), [senate data](R/rddensity_senate.csv).

## References

For overviews and introductions, see [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Cattaneo, Jansson and Ma (2018): [Manipulation Testing based on Density Discontinuity](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2018_Stata.pdf), _Stata Journal_ 18(1): 234-261.

- Cattaneo, Jansson and Ma (2020): [lpdensity: Local Polynomial Density Estimation and Inference.](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2020_JSS.pdf), working paper.

### Technical and Methodological

- Cattaneo, Jansson and Ma (2020): [Simple Local Polynomial Density Estimators](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA.pdf), _Journal of the American Statistical Association_, forthcoming. [Supplemental appendix](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2020_JASA--Supplement.pdf).

- Cattaneo, Jansson and Ma (2020): [Local Regression Distribution Estimators](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2020_JoE.pdf), working paper. [Supplemental appendix](https://rdpackages.github.io/references/Cattaneo-Jansson-Ma_2020_JoE--Supplement.pdf).

<br><br>
