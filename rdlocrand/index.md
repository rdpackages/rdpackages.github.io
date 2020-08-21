# RDLOCRAND

The `rdlocrand` package provides Stata and R implementations of statistical inference and graphical procedures for Regression Discontinuity designs employing local randomization methods. It provides point estimators, confidence intervals estimators, binomial manipulation testing, windows selectors, automatic plots, sensitivity analysis, and other related features.

This work was supported in part by the National Science Foundation through grants [SES-1357561](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1357561).

## Stata Implementation

To install/update in Stata type:
```
net install rdlocrand, from(https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata) replace
```

- Help: [rdrandinf](https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata/rdrandinf.pdf), [rdwinselect](https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata/rdwinselect.pdf), [rdsensitivity](https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata/rdsensitivity.pdf), [rdrbounds](https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata/rdrbounds.pdf).

- Replication: [do-file](https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata/rdlocrand_illustration.do), [senate data](https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata/rdlocrand_senate.dta).

## R Implementation

To install/update in R type:
```
install.packages('rdlocrand')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/rdlocrand/rdlocrand.pdf), [CRAN repository](https://cran.r-project.org/package=rdlocrand).

- Replication: [R-script](https://raw.githubusercontent.com/rdpackages/rdlocrand/master/R/rdlocrand_illustration.R), [senate data](R/rdlocrand_senate.csv). [R illustration](https://raw.githubusercontent.com/rdpackages/rdlocrand/master/R/rdlocrand_illustration.pdf).

## Repository

For source code and related files, visit [`rdlocrand` repository](https://github.com/rdpackages/rdlocrand/).


## References

For overviews and introductions, see [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Cattaneo, Titiunik and Vazquez-Bare (2016): [Inference in Regression Discontinuity Designs under Local Randomization](https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2016_Stata.pdf).<br>
_Stata Journal_ 16(2): 331-367.

### Technical and Methodological

- Cattaneo, Frandsen and Titiunik (2015): [Randomization Inference in the Regression Discontinuity Design: An Application to Party Advantages in the U.S. Senate](https://rdpackages.github.io/references/Cattaneo-Frandsen-Titiunik_2015_JCI.pdf).<br>
_Journal of Causal Inference_ 3(1): 1-24.

- Cattaneo, Titiunik and Vazquez-Bare (2017): [Comparing Inference Approaches for RD Designs: A Reexamination of the Effect of Head Start on Child Mortality](https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2017_JPAM.pdf).<br>
_Journal of Policy Analysis and Management_ 36(3): 643-681.

<br><br>



