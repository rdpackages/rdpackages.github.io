###########################################################################
## RD2D R Package
## Authors: Matias D. Cattaneo, Ruiqi (Rae) Yu and Rocio Titiunik 
###########################################################################
### Clear R environment
rm(list=ls(all=TRUE))

### Install R library
#install.packages('rd2d')

### Load RD2D package
library(rd2d)

### Load data base
rd2d_data <- read.csv("rd2d_data.csv")
attach(rd2d_data)

### Summary stats
summary(rd2d_data)


### Bivariate Location Score Estimation and Inference


### Univariate Distance Score Estimation and Inference


### Plot Results

