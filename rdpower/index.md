# RDPOWER

The `rdpower` package provides Stata and R implementations of power, sample size, and minimum detectable effects calculations using robust bias-corrected local polynomial inference methods.

This work was supported by the National Science Foundation through grant [SES-1357561](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1357561).

## Queries and Requests

Please email: [rdpackages@googlegroups.com](mailto:rdpackages@googlegroups.com)

## Major Upgrades

This package was first released in Fall 2016, and had one major upgrade in Fall 2020.

- _Fall 2020 new feature_: command/function `rdmde` for computing minimum detectable effects.

## Python Implementation

Coming soon.

## R Implementation

To install/update in R type:
```
install.packages('rdpower')
```
- Help: [R Manual](https://cran.r-project.org/web/packages/rdpower/rdpower.pdf), [CRAN repository](https://cran.r-project.org/package=rdpower).

- Replication files: [R-script](https://raw.githubusercontent.com/rdpackages/rdpower/master/R/rdpower_illustration.R), [data-senate](https://raw.githubusercontent.com/rdpackages/rdpower/master/R/rdpower_senate.csv).

## Stata Implementation

To install/update in Stata type:
```
net install rdpower, from(https://raw.githubusercontent.com/rdpackages/rdpower/master/stata) replace
```

- Help: [rdpower](https://raw.githubusercontent.com/rdpackages/rdpower/master/stata/rdpower.pdf), [rdsampsi](https://raw.githubusercontent.com/rdpackages/rdpower/master/stata/rdsampsi.pdf), [rdmde](https://raw.githubusercontent.com/rdpackages/rdpower/master/stata/rdmde.pdf).

- Replication: [do-file](https://raw.githubusercontent.com/rdpackages/rdpower/master/stata/rdpower_illustration.do), [data-senate](https://raw.githubusercontent.com/rdpackages/rdpower/master/stata/rdpower_senate.dta).

## Repository

For source code and related files, visit [`rdpower` repository](https://github.com/rdpackages/rdpower/).


## References

For overviews and introductions, see [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Cattaneo, Titiunik and Vazquez-Bare (2019): [Power Calculations for Regression Discontinuity Designs](https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2019_Stata.pdf).<br>
_Stata Journal_ 19(1): 210-245.

### Technical and Methodological

- Calonico, Cattaneo and Titiunik (2014): [Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf).<br>
_Econometrica_ 82(6): 2295-2326.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA--Supplemental.pdf).

- Calonico, Cattaneo, Farrell and Titiunik (2019): [Regression Discontinuity Designs Using Covariates](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf).<br>
_Review of Economics and Statistics_ 101(3): 442-451.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT--Supplement.pdf).

- Calonico, Cattaneo and Farrell (2020): [Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf).<br>
_Econometrics Journal_ 23(2): 192-210.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ--Supplement.pdf).

<br><br>
